"""
Motif-to-metadata mapping and clinical inference.

Core design
-----------
All prediction tasks share a single underlying step:

    donor × motif presence matrix  →  donor × label score matrix

``map_motifs_to_metadata()`` performs this step given any motif→label mapping
dictionary. Specialised wrappers (``score_pathogen_exposure``,
``score_hla``) populate that mapping from the bundled reference tables.

HLA inference
~~~~~~~~~~~~~
For each HLA allele, a Youden's-J-optimal threshold is learned from donors
with ground-truth HLA calls. Donors without ground truth are then scored
against those thresholds. Where no ground-truth cohort is available the
default threshold of 1 (≥ 1 matching motif = positive call) is applied.

Pathogen exposure inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~
A donor is called as exposed to a pathogen if their repertoire contains
≥ ``threshold`` distinct motifs associated with that pathogen. This simple
count approach achieves p < 10⁻¹⁷ for SARS-CoV-2 (chi-squared) in the
atlas cohort.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from sklearn.metrics import roc_curve

__all__ = [
    "build_donor_motif_matrix",
    "map_motifs_to_metadata",
    "score_pathogen_exposure",
    "predict_pathogen_exposure",
    "score_hla",
    "predict_hla_type",
    "evaluate_hla_prediction",
]


# ---------------------------------------------------------------------------
# Shared core
# ---------------------------------------------------------------------------

def build_donor_motif_matrix(
    clone_df: pd.DataFrame,
    donor_col: str = "donor_id",
    motif_col: str = "motif",
) -> pd.DataFrame:
    """
    Build a binary donor × motif presence/absence matrix.

    Parameters
    ----------
    clone_df
        DataFrame with one row per cell/clone. Must contain ``donor_col``
        and ``motif_col``.
    donor_col
        Column identifying donors.
    motif_col
        Column holding motif IDs. Rows with ``NaN`` motif are ignored.

    Returns
    -------
    pd.DataFrame
        Boolean (int 0/1) matrix, index = donors, columns = motif IDs.
    """
    df = clone_df.dropna(subset=[motif_col])
    return pd.crosstab(df[donor_col], df[motif_col]).gt(0).astype(int)


def map_motifs_to_metadata(
    donor_motif_matrix: pd.DataFrame,
    motif_to_label: Dict,
) -> pd.DataFrame:
    """
    Aggregate motif presence scores per donor into per-label scores.

    For each label, the score is the number of label-associated motifs
    present in a donor's repertoire.

    Parameters
    ----------
    donor_motif_matrix
        Binary donor × motif matrix (output of :func:`build_donor_motif_matrix`
        or any 0/1 DataFrame with motif IDs as columns).
    motif_to_label
        Mapping from motif ID → label string (pathogen name, HLA allele, etc.).
        Only motifs present in ``donor_motif_matrix.columns`` are used.

    Returns
    -------
    pd.DataFrame
        Donor × label integer count matrix.
    """
    label_to_motifs: Dict[str, List] = defaultdict(list)
    for motif, label in motif_to_label.items():
        if motif in donor_motif_matrix.columns:
            label_to_motifs[label].append(motif)

    scores = pd.DataFrame(index=donor_motif_matrix.index)
    for label, motifs in label_to_motifs.items():
        scores[label] = donor_motif_matrix[motifs].sum(axis=1)
    return scores


# ---------------------------------------------------------------------------
# Pathogen exposure
# ---------------------------------------------------------------------------

def score_pathogen_exposure(
    donor_motif_matrix: pd.DataFrame,
    motifs_disease: Optional[pd.DataFrame] = None,
    pathogen_col: str = "predicted_pathogen",
    motif_col: str = "motif",
) -> pd.DataFrame:
    """
    Compute per-donor pathogen exposure scores from motif counts.

    Parameters
    ----------
    donor_motif_matrix
        Binary donor × motif matrix.
    motifs_disease
        Pathogen-association table. If ``None``, the bundled table is loaded
        automatically via :func:`~._data.load_pathogen_motifs`.
    pathogen_col
        Column in ``motifs_disease`` with pathogen labels.
    motif_col
        Column in ``motifs_disease`` with motif IDs.

    Returns
    -------
    pd.DataFrame
        Donor × pathogen integer count matrix.
    """
    if motifs_disease is None:
        from ._data import load_pathogen_motifs
        motifs_disease = load_pathogen_motifs()

    motif_to_patho = (
        motifs_disease.dropna(subset=[pathogen_col])
        .set_index(motif_col)[pathogen_col]
        .to_dict()
    )
    return map_motifs_to_metadata(donor_motif_matrix, motif_to_patho)


def predict_pathogen_exposure(
    donor_motif_matrix: pd.DataFrame,
    motifs_disease: Optional[pd.DataFrame] = None,
    threshold: int = 2,
    pathogens: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Call binary pathogen exposure status per donor.

    A donor is considered exposed to a pathogen if they carry ≥ ``threshold``
    distinct motifs associated with that pathogen.

    Parameters
    ----------
    donor_motif_matrix
        Binary donor × motif matrix.
    motifs_disease
        Pathogen-association table (loaded automatically if ``None``).
    threshold
        Minimum number of pathogen-associated motifs to call a donor exposed.
        Threshold of 2 (the "double-hit" rule) provides high precision for
        SARS-CoV-2 (p=8×10⁻¹⁷, chi-squared, atlas cohort).
    pathogens
        Subset of pathogen labels to return. If ``None``, all pathogens in
        the reference table are included.

    Returns
    -------
    pd.DataFrame
        Boolean donor × pathogen exposure matrix.
    """
    scores = score_pathogen_exposure(donor_motif_matrix, motifs_disease)
    if pathogens is not None:
        scores = scores[[p for p in pathogens if p in scores.columns]]
    return scores.ge(threshold)


# ---------------------------------------------------------------------------
# HLA inference
# ---------------------------------------------------------------------------

def score_hla(
    donor_motif_matrix: pd.DataFrame,
    df_motifs_hla: Optional[pd.DataFrame] = None,
    min_donors: int = 4,
    cd8_only: bool = True,
) -> pd.DataFrame:
    """
    Compute per-donor HLA scores from MHC-I restricted motif counts.

    Parameters
    ----------
    donor_motif_matrix
        Binary donor × motif matrix.
    df_motifs_hla
        HLA motif table (loaded and filtered automatically if ``None``).
    min_donors
        Minimum donor count filter (passed to loader if ``df_motifs_hla`` is
        ``None``).
    cd8_only
        CD8-only filter (passed to loader if ``df_motifs_hla`` is ``None``).

    Returns
    -------
    pd.DataFrame
        Donor × HLA allele integer count matrix.
    """
    if df_motifs_hla is None:
        from ._data import load_hla_motifs
        df_motifs_hla = load_hla_motifs(min_donors=min_donors, cd8_only=cd8_only)

    motif_to_hla = (
        df_motifs_hla.set_index("motif")["MHC_I_restricted_allele"].to_dict()
    )
    return map_motifs_to_metadata(donor_motif_matrix, motif_to_hla)


def _youden_threshold(y_true: np.ndarray, scores: np.ndarray) -> float:
    """Return the Youden's J optimal threshold on training data."""
    fpr, tpr, thresholds = roc_curve(y_true, scores)
    thresholds = np.asarray(thresholds)
    j = tpr - fpr
    finite = np.isfinite(thresholds)
    if not finite.any():
        return 1.0
    return float(thresholds[finite][np.argmax(j[finite])])


def predict_hla_type(
    donor_motif_matrix: pd.DataFrame,
    df_motifs_hla: Optional[pd.DataFrame] = None,
    hla_ground_truth: Optional[pd.DataFrame] = None,
    default_threshold: float = 0.5,
    min_positives_for_training: int = 5,
) -> pd.DataFrame:
    """
    Predict HLA genotype per donor from TCR motif composition.

    For each HLA allele, a Youden's-J threshold is learned from the subset of
    donors in ``donor_motif_matrix`` who also appear in ``hla_ground_truth``.
    The threshold is then applied to all remaining donors.

    When ``hla_ground_truth`` is ``None`` (or too few training donors are
    available), ``default_threshold`` (≥ 1 matching motif = positive call) is
    used instead.

    Parameters
    ----------
    donor_motif_matrix
        Binary donor × motif matrix including all donors.
    df_motifs_hla
        HLA motif table (loaded automatically if ``None``).
    hla_ground_truth
        Donor × HLA allele binary ground-truth matrix (e.g. from ArcasHLA).
        Index must be donor IDs. If provided, per-allele thresholds are
        learned from donors present in both this table and
        ``donor_motif_matrix``.
    default_threshold
        Count threshold applied to alleles where training data are
        insufficient.
    min_positives_for_training
        Minimum number of positive training donors required to learn a
        Youden's threshold. Alleles with fewer positives fall back to
        ``default_threshold``.

    Returns
    -------
    pd.DataFrame
        Boolean donor × HLA allele prediction matrix.
    """
    hla_scores = score_hla(donor_motif_matrix, df_motifs_hla)

    thresholds: Dict[str, float] = {}

    if hla_ground_truth is not None:
        common_donors = hla_ground_truth.index.intersection(hla_scores.index)
        common_hla    = hla_ground_truth.columns.intersection(hla_scores.columns)

        for allele in common_hla:
            y = hla_ground_truth.loc[common_donors, allele]
            s = hla_scores.loc[common_donors, allele]
            if y.sum() >= min_positives_for_training and (len(y) - y.sum()) >= min_positives_for_training:
                thresholds[allele] = _youden_threshold(y.to_numpy(), s.to_numpy())
            else:
                thresholds[allele] = default_threshold

    predictions = pd.DataFrame(False, index=hla_scores.index, columns=hla_scores.columns)
    for allele in hla_scores.columns:
        thr = thresholds.get(allele, default_threshold)
        predictions[allele] = hla_scores[allele] >= thr

    return predictions


# ---------------------------------------------------------------------------
# Evaluation helpers
# ---------------------------------------------------------------------------

def evaluate_hla_prediction(
    hla_scores: pd.DataFrame,
    hla_ground_truth: pd.DataFrame,
    n_splits: int = 5,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Evaluate HLA inference via stratified K-fold cross-validation.

    Uses Youden's-J threshold selection within each training fold.

    Parameters
    ----------
    hla_scores
        Donor × HLA allele score matrix (output of :func:`score_hla`).
    hla_ground_truth
        Ground-truth binary HLA matrix aligned to the same donors/alleles.
    n_splits
        Number of cross-validation folds.
    random_state
        Random seed for fold splitting.

    Returns
    -------
    pd.DataFrame
        Per-allele cross-validated metrics: Mean_AUC, Std_AUC,
        Mean_BalAcc, Mean_F1, Mean_Precision, Mean_Recall, N_total,
        Pos_total.
    """
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import (
        roc_auc_score, balanced_accuracy_score, f1_score,
        precision_score, recall_score,
    )

    common_donors = hla_ground_truth.index.intersection(hla_scores.index)
    common_hla    = hla_ground_truth.columns.intersection(hla_scores.columns)

    X = hla_scores.loc[common_donors, common_hla]
    Y = hla_ground_truth.loc[common_donors, common_hla]

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    results = []

    for allele in common_hla:
        y = Y[allele]
        if y.sum() < n_splits or (len(y) - y.sum()) < n_splits:
            continue

        fold_metrics = []
        indices = np.arange(len(y))

        for train_raw, test_raw in skf.split(indices, y):
            train_idx = y.index[train_raw]
            test_idx  = y.index[test_raw]
            y_train, y_test = y.loc[train_idx], y.loc[test_idx]
            s_train, s_test = X.loc[train_idx, allele], X.loc[test_idx, allele]

            thr = _youden_threshold(y_train.to_numpy(), s_train.to_numpy())
            y_pred = (s_test >= thr).astype(int)

            fold_metrics.append({
                "AUC":       roc_auc_score(y_test, s_test) if y_test.nunique() > 1 else np.nan,
                "BalAcc":    balanced_accuracy_score(y_test, y_pred),
                "F1":        f1_score(y_test, y_pred, zero_division=0),
                "Precision": precision_score(y_test, y_pred, zero_division=0),
                "Recall":    recall_score(y_test, y_pred, zero_division=0),
            })

        fm = pd.DataFrame(fold_metrics)
        results.append({
            "HLA":            allele,
            "Mean_AUC":       fm["AUC"].mean(),
            "Std_AUC":        fm["AUC"].std(),
            "Mean_BalAcc":    fm["BalAcc"].mean(),
            "Mean_F1":        fm["F1"].mean(),
            "Mean_Precision": fm["Precision"].mean(),
            "Mean_Recall":    fm["Recall"].mean(),
            "N_total":        len(y),
            "Pos_total":      int(y.sum()),
        })

    return (
        pd.DataFrame(results)
        .set_index("HLA")
        .sort_values("Mean_AUC", ascending=False)
    )
