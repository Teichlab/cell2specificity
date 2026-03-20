"""
Fast approximate TCR matching via indexed seed-and-extend.

Algorithm overview
------------------
Reference preprocessing builds a per-chain index keyed by ``(V gene, length)``
(optionally also ``J gene``). For each key, CDR3 AA sequences are split into
``(mm+1)`` contiguous parts; for each part position an exact-substring hash
table maps substring → list of reference row indices. This exploits the
pigeonhole principle: if two strings of equal length differ by ≤ mm Hamming
mismatches, at least one of the ``(mm+1)`` parts must match exactly.

Querying repeats the split, collects candidate reference indices from matching
part buckets, then verifies candidates by computing Hamming distance ≤ mm over
the full CDR3 AA. Final paired-chain hits are the intersection of alpha and
beta per-query match sets.

Performance
-----------
Runs in <1 s for 1,000+ TCRs against a reference of 100,000+ clones.

Column convention
-----------------
Input DataFrames must have the compact column names produced by
:func:`cell2specificity.tcr_motifs._preprocess.to_matching_frame`:
``vdj_aa``, ``vdj_v``, ``vdj_j`` (beta chain) and
``vj_aa``, ``vj_v``, ``vj_j`` (alpha chain).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

__all__ = [
    "ChainIndex",
    "build_chain_index",
    "match_chain_all_reference_indices",
    "query_to_ref",
    "strip_allele",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def strip_allele(gene: str) -> str:
    """
    Convert e.g. ``'TRBV6-2*01'`` → ``'TRBV6-2'``.

    Handles ``None`` and ``nan`` safely.
    """
    if gene is None:
        return ""
    s = str(gene)
    return "" if s == "nan" else s.split("*", 1)[0]


def _hamming_leq(a: str, b: str, mm: int) -> bool:
    mismatches = 0
    for x, y in zip(a, b):
        if x != y:
            mismatches += 1
            if mismatches > mm:
                return False
    return True


def _split_positions(length: int, parts: int) -> List[Tuple[int, int]]:
    """Divide ``length`` into ``parts`` contiguous, balanced spans."""
    base, rem = divmod(length, parts)
    spans, start = [], 0
    for i in range(parts):
        end = start + base + (1 if i < rem else 0)
        spans.append((start, end))
        start = end
    return spans


# ---------------------------------------------------------------------------
# Index data structure
# ---------------------------------------------------------------------------

@dataclass
class ChainIndex:
    """
    Pre-built index for a single TCR chain.

    Attributes
    ----------
    by_key
        Nested mapping: key → list-of-part-dicts, where each part-dict maps
        CDR3 substring → list of reference row indices.
    include_j
        Whether J gene was included in the lookup key.
    chain
        ``'vdj'`` (beta) or ``'vj'`` (alpha).
    mm
        Maximum Hamming mismatches allowed.
    """
    by_key: Dict[Tuple, List[Dict[str, List[int]]]]
    include_j: bool
    chain: str
    mm: int


# ---------------------------------------------------------------------------
# Index construction
# ---------------------------------------------------------------------------

def build_chain_index(
    df_ref: pd.DataFrame,
    chain: str,
    mm: int = 1,
    include_j: bool = False,
) -> ChainIndex:
    """
    Build a seed-and-extend index for fast CDR3 matching on one chain.

    Parameters
    ----------
    df_ref
        Reference TCR DataFrame in compact column format
        (``vdj_aa``, ``vdj_v``, ``vdj_j``, ``vj_aa``, ``vj_v``, ``vj_j``).
    chain
        ``'vdj'`` for beta chain, ``'vj'`` for alpha chain.
    mm
        Maximum Hamming mismatches. ``mm+1`` parts are used for seeding.
    include_j
        If ``True``, include J gene (allele-stripped) in the lookup key,
        producing slightly more specific matching.

    Returns
    -------
    ChainIndex
    """
    if mm < 0:
        raise ValueError("mm must be >= 0")
    parts = mm + 1

    aa  = df_ref[f"{chain}_aa"].astype(str).to_numpy()
    v   = df_ref[f"{chain}_v"].map(strip_allele).astype(str).to_numpy()
    j   = df_ref[f"{chain}_j"].map(strip_allele).astype(str).to_numpy() if include_j else None

    by_key: Dict[Tuple, List[Dict[str, List[int]]]] = {}

    for ridx, (s, vi) in enumerate(zip(aa, v)):
        if not s or s == "nan":
            continue
        L = len(s)
        key = (vi, j[ridx], L) if include_j else (vi, L)  # type: ignore[index]
        if key not in by_key:
            by_key[key] = [dict() for _ in range(parts)]

        for p, (a0, a1) in enumerate(_split_positions(L, parts)):
            sub = s[a0:a1]
            bucket = by_key[key][p]
            if sub in bucket:
                bucket[sub].append(ridx)
            else:
                bucket[sub] = [ridx]

    return ChainIndex(by_key=by_key, include_j=include_j, chain=chain, mm=mm)


# ---------------------------------------------------------------------------
# Querying
# ---------------------------------------------------------------------------

def match_chain_all_reference_indices(
    df_ref: pd.DataFrame,
    df_query: pd.DataFrame,
    index: ChainIndex,
    max_candidates: int = 10_000,
) -> List[List[int]]:
    """
    For each query row, return all reference row indices that match.

    Match criteria:
      - Same V gene (allele-stripped)
      - Same CDR3 length
      - Same J gene (allele-stripped) if ``index.include_j``
      - Hamming distance on CDR3 ≤ ``index.mm``

    Parameters
    ----------
    df_ref
        Reference DataFrame (same used to build ``index``).
    df_query
        Query DataFrame in compact column format.
    index
        Pre-built :class:`ChainIndex` for this chain.
    max_candidates
        Safety cap on candidate size per query (prevents slow queries against
        very common seeds).

    Returns
    -------
    list of lists
        One inner list of reference row indices per query row.
    """
    chain = index.chain
    mm    = index.mm
    parts = mm + 1

    r_aa = df_ref[f"{chain}_aa"].astype(str).to_numpy()
    q_aa = df_query[f"{chain}_aa"].astype(str).to_numpy()
    q_v  = df_query[f"{chain}_v"].map(strip_allele).astype(str).to_numpy()
    q_j  = df_query[f"{chain}_j"].map(strip_allele).astype(str).to_numpy() if index.include_j else None

    out: List[List[int]] = []

    for qidx, (s, vi) in enumerate(zip(q_aa, q_v)):
        if not s or s == "nan":
            out.append([])
            continue
        L = len(s)
        key = (vi, q_j[qidx], L) if index.include_j else (vi, L)  # type: ignore[index]
        part_dicts = index.by_key.get(key)
        if part_dicts is None:
            out.append([])
            continue

        candidates: Set[int] = set()
        for p, (a0, a1) in enumerate(_split_positions(L, parts)):
            hits = part_dicts[p].get(s[a0:a1])
            if hits:
                candidates.update(hits)
                if len(candidates) > max_candidates:
                    break

        matches: List[int] = [
            ridx for ridx in candidates
            if len(r_aa[ridx]) == L and _hamming_leq(s, r_aa[ridx], mm)
        ]
        out.append(matches)

    return out


def query_to_ref(
    df_ref: pd.DataFrame,
    df_query: pd.DataFrame,
    idx_alpha: ChainIndex,
    idx_beta: ChainIndex,
    motif_col: str = "motif",
) -> pd.DataFrame:
    """
    Map query TCRs to a reference atlas and annotate with matching motifs.

    Computes paired-chain hits (intersection of alpha and beta matches) and
    annotates ``df_query`` with:

    * ``ref_hits`` — list of reference row indices for each query TCR
    * ``n_ref_hits`` — count of matching reference TCRs
    * ``ref_motifs`` — unique motif IDs from matching reference TCRs
    * ``n_ref_motifs`` — number of distinct motifs

    Parameters
    ----------
    df_ref
        Reference TCR DataFrame with a ``motif`` column (or ``motif_col``).
    df_query
        Query TCR DataFrame in compact column format.
    idx_alpha
        Pre-built alpha-chain :class:`ChainIndex`.
    idx_beta
        Pre-built beta-chain :class:`ChainIndex`.
    motif_col
        Column in ``df_ref`` holding motif IDs.

    Returns
    -------
    pd.DataFrame
        A copy of ``df_query`` with the four annotation columns added.
    """
    beta_hits  = match_chain_all_reference_indices(df_ref, df_query, idx_beta)
    alpha_hits = match_chain_all_reference_indices(df_ref, df_query, idx_alpha)

    paired = [
        list(set(b) & set(a)) if b and a else []
        for b, a in zip(beta_hits, alpha_hits)
    ]

    df_query = df_query.copy()
    df_query["ref_hits"]     = paired
    df_query["n_ref_hits"]   = [len(h) for h in paired]
    df_query["ref_motifs"]   = [
        list(set(df_ref.iloc[h][motif_col].tolist())) if h else []
        for h in paired
    ]
    df_query["n_ref_motifs"] = [len(m) for m in df_query["ref_motifs"]]

    return df_query
