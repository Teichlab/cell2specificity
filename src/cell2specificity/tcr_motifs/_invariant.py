"""
Annotation of invariant T cell receptors (MAIT and iNKT).

MAIT cells are defined by co-expression of TRAV1-2 paired with TRAJ33,
TRAJ20, or TRAJ12.  iNKT cells are defined by TRAV10 paired with TRAJ18.
All other cells are labelled 'Conventional'.

Classification is applied row-wise and is independent of any clustering or
transcriptomic data, relying solely on the V/J gene calls in the six standard
IR_ VDJ columns.

Reference
---------
Dratva et al. (2026) — Methods: "Assigning MAIT and iNKT Cell Receptors"
"""

from __future__ import annotations

from typing import Union

import pandas as pd

__all__ = [
    "MAIT_V",
    "MAIT_J",
    "INKT_V",
    "INKT_J",
    "classify_invariant_row",
    "annotate_invariant",
]

# ---------------------------------------------------------------------------
# Gene constants — easy to override if studying other organisms / alleles
# ---------------------------------------------------------------------------

#: Alpha-chain V gene prefix defining MAIT cells.
MAIT_V: str = "TRAV1-2"

#: Alpha-chain J genes defining MAIT cells (any of these).
MAIT_J: tuple[str, ...] = ("TRAJ33", "TRAJ20", "TRAJ12")

#: Alpha-chain V gene prefix defining iNKT cells.
INKT_V: str = "TRAV10"

#: Alpha-chain J gene defining iNKT cells.
INKT_J: str = "TRAJ18"

# Label constants
_MAIT  = "MAIT"
_INKT  = "iNKT"
_CONV  = "Conventional"


# ---------------------------------------------------------------------------
# Row-level classifier
# ---------------------------------------------------------------------------

def classify_invariant_row(
    vj_v: Union[str, float],
    vj_j: Union[str, float],
    mait_v: str = MAIT_V,
    mait_j: tuple = MAIT_J,
    inkt_v: str = INKT_V,
    inkt_j: str = INKT_J,
) -> str:
    """
    Classify a single cell as MAIT, iNKT, or Conventional.

    Parameters
    ----------
    vj_v
        Alpha-chain V gene call (e.g. ``'TRAV1-2*01'``).
    vj_j
        Alpha-chain J gene call (e.g. ``'TRAJ33*01'``).
    mait_v, mait_j, inkt_v, inkt_j
        Gene definitions. Override to study other organisms.

    Returns
    -------
    str
        ``'MAIT'``, ``'iNKT'``, or ``'Conventional'``.
    """
    # Treat missing values as Conventional
    v = str(vj_v) if pd.notna(vj_v) else ""
    j = str(vj_j) if pd.notna(vj_j) else ""

    if mait_v in v:
        if any(jgene in j for jgene in mait_j):
            return _MAIT
    elif inkt_v in v:
        if inkt_j in j:
            return _INKT

    return _CONV


# ---------------------------------------------------------------------------
# DataFrame-level annotation
# ---------------------------------------------------------------------------

def annotate_invariant(
    df: pd.DataFrame,
    vj_v_col: str = "IR_VJ_1_v_call",
    vj_j_col: str = "IR_VJ_1_j_call",
    out_col: str = "invariant",
    mait_v: str = MAIT_V,
    mait_j: tuple = MAIT_J,
    inkt_v: str = INKT_V,
    inkt_j: str = INKT_J,
) -> pd.DataFrame:
    """
    Annotate a TCR DataFrame with invariant T cell identity.

    Adds a new column (default ``'invariant'``) with values
    ``'MAIT'``, ``'iNKT'``, or ``'Conventional'`` for each cell.

    Parameters
    ----------
    df
        TCR DataFrame with at minimum the alpha-chain V and J gene columns.
    vj_v_col
        Column holding the alpha-chain V gene call.
        Default: ``'IR_VJ_1_v_call'``.
    vj_j_col
        Column holding the alpha-chain J gene call.
        Default: ``'IR_VJ_1_j_call'``.
    out_col
        Name of the output annotation column. Default: ``'invariant'``.
    mait_v
        V gene prefix for MAIT classification. Default: ``'TRAV1-2'``.
    mait_j
        Tuple of J gene substrings for MAIT classification.
        Default: ``('TRAJ33', 'TRAJ20', 'TRAJ12')``.
    inkt_v
        V gene prefix for iNKT classification. Default: ``'TRAV10'``.
    inkt_j
        J gene substring for iNKT classification. Default: ``'TRAJ18'``.

    Returns
    -------
    pd.DataFrame
        A copy of ``df`` with the ``out_col`` column added.

    Examples
    --------
    >>> from cell2specificity.tcr_motifs import annotate_invariant
    >>> df = annotate_invariant(df)
    >>> df['invariant'].value_counts()
    Conventional    180423
    MAIT              8712
    iNKT               423
    Name: invariant, dtype: int64
    """
    missing = [c for c in [vj_v_col, vj_j_col] if c not in df.columns]
    if missing:
        raise ValueError(
            f"Column(s) not found in DataFrame: {missing}. "
            f"Use vj_v_col / vj_j_col parameters to specify the correct names."
        )

    df = df.copy()
    df[out_col] = df.apply(
        lambda row: classify_invariant_row(
            row[vj_v_col], row[vj_j_col],
            mait_v=mait_v, mait_j=mait_j,
            inkt_v=inkt_v, inkt_j=inkt_j,
        ),
        axis=1,
    )
    return df


def invariant_summary(
    df: pd.DataFrame,
    invariant_col: str = "invariant",
    groupby: Union[str, None] = None,
    normalize: bool = True,
) -> pd.DataFrame:
    """
    Summarise the proportion of MAIT, iNKT, and Conventional T cells.

    Parameters
    ----------
    df
        DataFrame with an ``invariant_col`` column (from :func:`annotate_invariant`).
    invariant_col
        Column holding invariant labels.
    groupby
        Optional metadata column to stratify by (e.g. ``'pathogen'``,
        ``'tissue'``, ``'donor_id'``).
    normalize
        Return proportions (``True``) or raw counts (``False``).

    Returns
    -------
    pd.DataFrame
        Counts or proportions of each invariant category, optionally
        stratified by ``groupby``.

    Examples
    --------
    >>> summary = invariant_summary(df, groupby='pathogen')
    """
    if invariant_col not in df.columns:
        raise ValueError(
            f"Column '{invariant_col}' not found. Run annotate_invariant() first."
        )

    if groupby is not None:
        if groupby not in df.columns:
            raise ValueError(f"groupby column '{groupby}' not found in DataFrame.")
        counts = pd.crosstab(df[groupby], df[invariant_col])
    else:
        counts = df[invariant_col].value_counts().to_frame("count").T

    if normalize:
        return counts.div(counts.sum(axis=1), axis=0)
    return counts
