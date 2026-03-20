"""
Load the reference motif metadata tables bundled with the package.

Two tables are shipped:

``disease_associated_motifs_hla.csv``
    One row per TCR motif. Key column: ``predicted_pathogen``.
    Used for pathogen exposure inference.

``df_motifs_with_hla.csv``
    One row per TCR motif with HLA restriction annotations.
    Key columns: ``MHC_I_restricted``, ``MHC_I_restricted_allele``,
    ``is_CD8``, ``n_donors``.
    Used for HLA genotype inference.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd

_DATA_DIR = Path(__file__).parent / "data"

_PATHOGEN_FILE = _DATA_DIR / "disease_associated_motifs_hla.csv"
_HLA_FILE      = _DATA_DIR / "df_motifs_with_hla.csv"


@lru_cache(maxsize=1)
def load_pathogen_motifs() -> pd.DataFrame:
    """
    Load the pathogen-association table.

    Returns
    -------
    pd.DataFrame
        Columns include ``motif``, ``predicted_pathogen`` and per-pathogen
        association metadata.
    """
    _require(_PATHOGEN_FILE, "disease_associated_motifs_hla.csv")
    return pd.read_csv(_PATHOGEN_FILE)


@lru_cache(maxsize=1)
def load_hla_motifs(
    min_donors: int = 4,
    cd8_only: bool = True,
    exclude_multiple: bool = True,
) -> pd.DataFrame:
    """
    Load and filter the HLA-restriction table.

    Parameters
    ----------
    min_donors
        Minimum number of donors a motif must span to be included.
    cd8_only
        Keep only motifs with ``is_CD8 > 0.5``.
    exclude_multiple
        Exclude motifs whose restriction is labelled ``'Multiple'``.

    Returns
    -------
    pd.DataFrame
        Filtered HLA motif table.
    """
    _require(_HLA_FILE, "df_motifs_with_hla.csv")
    df = pd.read_csv(_HLA_FILE)
    if cd8_only:
        df = df[df["is_CD8"] > 0.5]
    df = df[df["MHC_I_restricted"] == True]  # noqa: E712
    df = df[df["n_donors"] >= min_donors]
    if exclude_multiple:
        df = df[df["MHC_I_restricted_allele"] != "Multiple"]
    return df


def _require(path: Path, name: str) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Bundled data file '{name}' not found at {path}.\n"
            "Place the CSV files in src/cell2specificity/specificity/data/."
        )
