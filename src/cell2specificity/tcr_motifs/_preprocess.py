"""
Preprocessing helpers for TCR input data.

Standardises the IR_ column convention used throughout the atlas and
validates V/J gene nomenclature before motif inference or matching.
"""

from __future__ import annotations

import pandas as pd

__all__ = ["VDJ_COLS", "add_allele_suffix", "preprocess_tcr_table", "to_matching_frame"]

#: Canonical six VDJ column names used across this package.
VDJ_COLS = [
    "IR_VDJ_1_v_call",
    "IR_VDJ_1_j_call",
    "IR_VDJ_1_junction_aa",
    "IR_VJ_1_v_call",
    "IR_VJ_1_j_call",
    "IR_VJ_1_junction_aa",
]


def add_allele_suffix(val) -> str:
    """
    Append ``*01`` to a gene name that has no allele suffix.

    Handles ``None``, ``float('nan')``, and empty strings gracefully.
    """
    s = str(val).strip()
    if s.lower() in ("nan", "none", ""):
        return val
    return val if "*" in s else f"{s}*01"


def preprocess_tcr_table(
    df: pd.DataFrame,
    add_alleles: bool = True,
    drop_incomplete: bool = False,
) -> pd.DataFrame:
    """
    Validate and standardise a TCR table for downstream use.

    Parameters
    ----------
    df
        DataFrame that must contain the six canonical VDJ columns
        (see :data:`VDJ_COLS`).
    add_alleles
        Append ``*01`` to V/J gene names that lack an allele suffix.
    drop_incomplete
        Drop rows where any VDJ column is null.

    Returns
    -------
    pd.DataFrame
        A copy of ``df`` with standardised VDJ columns.
    """
    missing = [c for c in VDJ_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"DataFrame is missing required VDJ columns: {missing}")

    df = df.copy()

    if add_alleles:
        gene_cols = [c for c in VDJ_COLS if c.endswith("_call")]
        for col in gene_cols:
            df[col] = df[col].apply(add_allele_suffix)

    if drop_incomplete:
        df = df.dropna(subset=VDJ_COLS)

    return df


def to_matching_frame(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename IR_ columns to the compact ``vdj_*`` / ``vj_*`` convention
    required by :mod:`cell2specificity.tcr_motifs._matching`.

    Parameters
    ----------
    df
        DataFrame with canonical IR_ VDJ columns.

    Returns
    -------
    pd.DataFrame
        A view with added ``vdj_aa``, ``vdj_v``, ``vdj_j``,
        ``vj_aa``, ``vj_v``, ``vj_j`` columns.
    """
    rename_map = {
        "IR_VDJ_1_junction_aa": "vdj_aa",
        "IR_VDJ_1_v_call":      "vdj_v",
        "IR_VDJ_1_j_call":      "vdj_j",
        "IR_VJ_1_junction_aa":  "vj_aa",
        "IR_VJ_1_v_call":       "vj_v",
        "IR_VJ_1_j_call":       "vj_j",
    }
    return df.rename(columns=rename_map)
