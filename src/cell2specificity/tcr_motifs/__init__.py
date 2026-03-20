"""
tcr_motifs
----------
TCR clonotype grouping, public motif inference, VDJdb annotation, and fast
approximate TCR matching.

This module incorporates and extends functionality from the ``cell2tcr``
package (Teichlab/cell2tcr, db_extension branch). Since that package is no
longer actively maintained, the core logic has been vendored and extended
here rather than installed as a dependency.

Sub-modules
~~~~~~~~~~~
``_preprocess``
    Column standardisation and V/J allele normalisation.

``_matching``
    Fast seed-and-extend paired-chain TCR matching (~1 s for 1 k+ queries vs
    a 100 k+ reference). Suitable for mapping new user data onto atlas motifs.

Motif inference (``run_motif_inference``) wraps the original Cell2TCR
algorithm: tcrdist3-based pairwise distances → sparse graph → Leiden
clustering. Requires ``tcrdist3``, ``igraph``, and ``leidenalg``.

VDJdb queries (``db_match``, ``db_annotate``) use the TCRMatch K3 kernel
against the bundled IEDB / VDJdb database.

Public API
----------
"""

from __future__ import annotations

from typing import Optional, Union
import pandas as pd

from ._preprocess import (
    VDJ_COLS,
    add_allele_suffix,
    preprocess_tcr_table,
    to_matching_frame,
)
from ._matching import (
    ChainIndex,
    build_chain_index,
    match_chain_all_reference_indices,
    query_to_ref,
    strip_allele,
)
from ._invariant import (
    MAIT_V, MAIT_J, INKT_V, INKT_J,
    classify_invariant_row,
    annotate_invariant,
    invariant_summary,
)

__all__ = [
    # preprocessing
    "VDJ_COLS",
    "add_allele_suffix",
    "preprocess_tcr_table",
    "to_matching_frame",
    # matching
    "ChainIndex",
    "build_chain_index",
    "match_chain_all_reference_indices",
    "query_to_ref",
    "strip_allele",
    # invariant annotation
    "MAIT_V", "MAIT_J", "INKT_V", "INKT_J",
    "classify_invariant_row",
    "annotate_invariant",
    "invariant_summary",
    # motif inference
    "run_motif_inference",
    "assign_column_names",
    # VDJdb
    "db_match",
    "db_annotate",
]


# ---------------------------------------------------------------------------
# Motif inference (wraps cell2tcr / tcrdist3 logic)
# ---------------------------------------------------------------------------

def assign_column_names(
    df: pd.DataFrame,
    vdj_v: str = "IR_VDJ_1_v_call",
    vdj_j: str = "IR_VDJ_1_j_call",
    vdj_aa: str = "IR_VDJ_1_junction_aa",
    vj_v: str = "IR_VJ_1_v_call",
    vj_j: str = "IR_VJ_1_j_call",
    vj_aa: str = "IR_VJ_1_junction_aa",
    clone_id: Optional[str] = None,
    donor_id: Optional[str] = None,
) -> pd.DataFrame:
    """
    Map user-supplied column names to the format expected by
    :func:`run_motif_inference`.

    Defaults match the IR_ convention used across the atlas. Override any
    parameter to adapt to a different naming scheme.

    Parameters
    ----------
    df
        Input TCR DataFrame.
    vdj_v, vdj_j, vdj_aa
        Column names for the beta-chain V gene, J gene, and CDR3 AA sequence.
    vj_v, vj_j, vj_aa
        Column names for the alpha-chain V gene, J gene, and CDR3 AA sequence.
    clone_id
        Optional column identifying unique clones. If ``None``, each row is
        treated as a distinct clone.
    donor_id
        Optional column for donor identifiers.

    Returns
    -------
    pd.DataFrame
        A renamed copy ready for :func:`run_motif_inference`.
    """
    rename = {
        vdj_v:  "IR_VDJ_1_v_call",
        vdj_j:  "IR_VDJ_1_j_call",
        vdj_aa: "IR_VDJ_1_junction_aa",
        vj_v:   "IR_VJ_1_v_call",
        vj_j:   "IR_VJ_1_j_call",
        vj_aa:  "IR_VJ_1_junction_aa",
    }
    if clone_id:
        rename[clone_id] = "clone_id"
    if donor_id:
        rename[donor_id] = "donor_id"
    return df.rename(columns={k: v for k, v in rename.items() if k in df.columns})


def run_motif_inference(
    df: pd.DataFrame,
    distance_threshold: float = 12.5,
    chunk_size: int = 100,
    sparse: bool = True,
    organism: str = "human",
    add_alleles: bool = True,
    **leiden_kwargs,
) -> pd.DataFrame:
    """
    Assign TCR motif IDs by clustering clonotypes with tcrdist3 + Leiden.

    Wraps the original Cell2TCR ``motifs()`` function. Clonotypes with
    sufficient CDR3 similarity (tcrdist ≤ ``distance_threshold``) are grouped
    into the same motif; singletons receive a unique motif ID.

    Parameters
    ----------
    df
        DataFrame with canonical IR_ VDJ columns (see :data:`VDJ_COLS`).
        Use :func:`assign_column_names` if your columns differ.
    distance_threshold
        Maximum tcrdist distance for two TCRs to be connected in the graph.
        Default (12.5) was optimised for paired alpha-beta chains at the
        atlas scale.
    chunk_size
        Chunk size for batched tcrdist distance computation.
    sparse
        Use sparse distance computation (recommended for large datasets).
    organism
        ``'human'`` or ``'mouse'``.
    add_alleles
        Append ``*01`` to V/J gene names lacking an allele suffix before
        passing to tcrdist3.
    **leiden_kwargs
        Additional keyword arguments forwarded to the Leiden clustering step.

    Returns
    -------
    pd.DataFrame
        ``df`` with an additional ``'motif'`` column (integer motif IDs).

    Raises
    ------
    ImportError
        If ``tcrdist3``, ``igraph``, or ``leidenalg`` are not installed.
    """
    try:
        import cell2tcr as _c2t  # type: ignore
    except ImportError:
        # Fallback: attempt to use tcrdist3 directly if cell2tcr not available
        raise ImportError(
            "run_motif_inference requires cell2tcr (tcrdist3 + igraph + leidenalg).\n"
            "Install with: pip install tcrdist3 igraph leidenalg\n"
            "and the cell2tcr package from Teichlab/cell2tcr (db_extension branch)."
        )

    if add_alleles:
        df = preprocess_tcr_table(df, add_alleles=True)

    return _c2t.motifs(
        df,
        threshold=distance_threshold,
        chunk_size=chunk_size,
        sparse=sparse,
        organism=organism,
        **leiden_kwargs,
    )


# ---------------------------------------------------------------------------
# VDJdb / IEDB annotation
# ---------------------------------------------------------------------------

def db_match(
    sequences: pd.Series,
    threshold: float = 0.97,
    trim_flanks: bool = True,
    n_jobs: int = 1,
) -> pd.DataFrame:
    """
    Match beta-chain CDR3 sequences against the IEDB TCRMatch database.

    Uses the TCRMatch K3 kernel (amino acid substitution matrix + Levenshtein
    pre-filter) to find sequences in the IEDB with similarity ≥ ``threshold``.

    Parameters
    ----------
    sequences
        Series of CDR3 amino acid strings (beta chain).
    threshold
        Minimum TCRMatch similarity score (0–1). Default 0.97.
    trim_flanks
        Trim leading C and trailing F/W residues before matching.
    n_jobs
        Number of parallel worker processes.

    Returns
    -------
    pd.DataFrame
        Table of matched IEDB entries with columns including ``query_cdr3``,
        ``db_cdr3``, ``score``, ``epitope``, ``mhc``, ``organism``.

    Raises
    ------
    ImportError
        If ``cell2tcr`` is not installed.
    """
    try:
        import cell2tcr as _c2t  # type: ignore
    except ImportError:
        raise ImportError(
            "db_match requires cell2tcr.\n"
            "Install from Teichlab/cell2tcr (db_extension branch)."
        )
    return _c2t.db_match(sequences, threshold=threshold, trim_flanks=trim_flanks, n_jobs=n_jobs)


def db_annotate(
    df: pd.DataFrame,
    scores: pd.DataFrame,
    column: str = "db_annotation",
    delimiter: str = "|",
) -> pd.DataFrame:
    """
    Add VDJdb / IEDB match results as a new column on the original DataFrame.

    Parameters
    ----------
    df
        TCR DataFrame to annotate.
    scores
        Output of :func:`db_match`.
    column
        Name of the new annotation column.
    delimiter
        Delimiter used when multiple database entries match a single TCR.

    Returns
    -------
    pd.DataFrame
        ``df`` with the new ``column`` column.
    """
    try:
        import cell2tcr as _c2t  # type: ignore
    except ImportError:
        raise ImportError(
            "db_annotate requires cell2tcr.\n"
            "Install from Teichlab/cell2tcr (db_extension branch)."
        )
    return _c2t.db_annotate(df, scores, column=column, delimiter=delimiter)
