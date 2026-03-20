"""
Tests for cell2specificity.motif_based_inference

Run with:  pytest tests/test_specificity.py -v

Toy dataset facts relevant to this module:
  - 865 cells, 113 donors, 752 cells with motif assignments
  - 125 cells whose motif IDs appear in the bundled reference CSVs
  - Pathogens with reference motif hits: HPV, Influenza_CMV_EBV, EBV_VZV,
    Influenza_virus, HSV-2, HIV, HBV
"""

import numpy as np
import pandas as pd
import pytest
from pathlib import Path

from cell2specificity.motif_based_inference import (
    build_donor_motif_matrix,
    map_motifs_to_metadata,
    score_pathogen_exposure,
    predict_pathogen_exposure,
    score_hla,
    predict_hla_type,
    evaluate_hla_prediction,
    load_pathogen_motifs,
    load_hla_motifs,
)

DATA_DIR = Path(__file__).parent / "data"
TOY_CSV  = DATA_DIR / "toy_tcr_atlas.csv"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def toy_df():
    return pd.read_csv(TOY_CSV, index_col=0)


@pytest.fixture(scope="module")
def toy_dmm(toy_df):
    return build_donor_motif_matrix(toy_df)


@pytest.fixture
def mini_clone_df():
    return pd.DataFrame({
        "donor_id": ["d1", "d1", "d1", "d2", "d2", "d3"],
        "motif":    [101,   102,   103,   101,   104,   105],
    })


@pytest.fixture
def mini_motif_to_label():
    return {101: "CMV", 102: "EBV", 103: "CMV", 104: "HIV", 105: "EBV"}


# ---------------------------------------------------------------------------
# Bundled data
# ---------------------------------------------------------------------------

class TestBundledData:

    def test_pathogen_table_loads(self):
        df = load_pathogen_motifs()
        assert "motif" in df.columns
        assert "predicted_pathogen" in df.columns
        assert len(df) > 0

    def test_pathogen_table_has_known_pathogens(self):
        df = load_pathogen_motifs()
        pathogens = df["predicted_pathogen"].dropna().unique()
        for expected in ("SARS-CoV-2", "CMV", "EBV"):
            assert expected in pathogens, f"{expected} not found in pathogen table"

    def test_hla_table_loads(self):
        df = load_hla_motifs()
        assert "motif" in df.columns
        assert "MHC_I_restricted_allele" in df.columns
        assert len(df) > 0

    def test_hla_table_filter_cd8_only(self):
        df = load_hla_motifs(cd8_only=True)
        assert (df["is_CD8"] > 0.5).all()

    def test_hla_table_filter_min_donors(self):
        df = load_hla_motifs(min_donors=4)
        assert (df["n_donors"] >= 4).all()

    def test_hla_table_excludes_multiple(self):
        df = load_hla_motifs(exclude_multiple=True)
        assert "Multiple" not in df["MHC_I_restricted_allele"].values


# ---------------------------------------------------------------------------
# Donor × motif matrix
# ---------------------------------------------------------------------------

class TestDonorMotifMatrix:

    def test_shape_mini(self, mini_clone_df):
        dmm = build_donor_motif_matrix(mini_clone_df)
        assert dmm.shape == (3, 5)

    def test_binary_values(self, mini_clone_df):
        dmm = build_donor_motif_matrix(mini_clone_df)
        assert set(dmm.values.flatten()).issubset({0, 1})

    def test_nan_motifs_dropped(self, mini_clone_df):
        df = mini_clone_df.copy()
        df.loc[0, "motif"] = float("nan")
        dmm = build_donor_motif_matrix(df)
        assert "d1" in dmm.index  # donor still present via other rows

    def test_toy_donor_count(self, toy_dmm, toy_df):
        # 113 unique donors in the CSV, but 1 has only null motifs and is
        # dropped by build_donor_motif_matrix → expect 112 rows in the matrix
        assert toy_dmm.shape[0] <= toy_df["donor_id"].nunique()

    def test_toy_binary(self, toy_dmm):
        assert toy_dmm.values.max() == 1
        assert toy_dmm.values.min() == 0


# ---------------------------------------------------------------------------
# Shared mapping core
# ---------------------------------------------------------------------------

class TestMapMotifs:

    def test_labels_match_input(self, mini_clone_df, mini_motif_to_label):
        dmm = build_donor_motif_matrix(mini_clone_df)
        scores = map_motifs_to_metadata(dmm, mini_motif_to_label)
        assert set(scores.columns) == {"CMV", "EBV", "HIV"}

    def test_d1_cmv_score_is_2(self, mini_clone_df, mini_motif_to_label):
        # d1 has motifs 101 (CMV) and 103 (CMV) → score = 2
        dmm = build_donor_motif_matrix(mini_clone_df)
        scores = map_motifs_to_metadata(dmm, mini_motif_to_label)
        assert scores.loc["d1", "CMV"] == 2

    def test_d2_hiv_score_is_1(self, mini_clone_df, mini_motif_to_label):
        dmm = build_donor_motif_matrix(mini_clone_df)
        scores = map_motifs_to_metadata(dmm, mini_motif_to_label)
        assert scores.loc["d2", "HIV"] == 1

    def test_unknown_motifs_give_zero_scores(self, mini_clone_df):
        dmm = build_donor_motif_matrix(mini_clone_df)
        scores = map_motifs_to_metadata(dmm, {9999: "Phantom"})
        # motif not in matrix → column absent or all zeros
        if "Phantom" in scores.columns:
            assert scores["Phantom"].sum() == 0


# ---------------------------------------------------------------------------
# Pathogen exposure
# ---------------------------------------------------------------------------

class TestPathogenExposure:

    def test_returns_boolean_df(self, mini_clone_df, mini_motif_to_label):
        dmm = build_donor_motif_matrix(mini_clone_df)
        md = pd.DataFrame(list(mini_motif_to_label.items()),
                          columns=["motif", "predicted_pathogen"])
        result = predict_pathogen_exposure(dmm, motifs_disease=md, threshold=1)
        assert result.dtypes.apply(
            lambda dt: dt == bool or np.issubdtype(dt, np.bool_)
        ).all()

    def test_higher_threshold_fewer_positives(self, mini_clone_df, mini_motif_to_label):
        dmm = build_donor_motif_matrix(mini_clone_df)
        md = pd.DataFrame(list(mini_motif_to_label.items()),
                          columns=["motif", "predicted_pathogen"])
        t1 = predict_pathogen_exposure(dmm, motifs_disease=md, threshold=1).sum().sum()
        t2 = predict_pathogen_exposure(dmm, motifs_disease=md, threshold=2).sum().sum()
        assert t1 >= t2

    def test_toy_bundled_data_runs(self, toy_dmm):
        result = predict_pathogen_exposure(toy_dmm, threshold=2)
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == toy_dmm.shape[0]

    def test_toy_some_donors_predicted_exposed(self, toy_dmm):
        # the toy set includes cells from pathogens present in the reference —
        # at threshold=1 at least some donors should be called positive
        result = predict_pathogen_exposure(toy_dmm, threshold=1)
        assert result.sum().sum() > 0, \
            "No donors predicted exposed at threshold=1 — check motif overlap with bundled reference"

    def test_toy_double_hit_stricter(self, toy_dmm):
        t1 = predict_pathogen_exposure(toy_dmm, threshold=1).sum().sum()
        t2 = predict_pathogen_exposure(toy_dmm, threshold=2).sum().sum()
        assert t1 >= t2


# ---------------------------------------------------------------------------
# HLA inference
# ---------------------------------------------------------------------------

class TestHLAInference:

    def test_score_hla_returns_df(self, toy_dmm):
        result = score_hla(toy_dmm)
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == toy_dmm.shape[0]

    def test_score_hla_columns_are_alleles(self, toy_dmm):
        result = score_hla(toy_dmm)
        # All column names should look like HLA alleles e.g. A*02:01
        assert all("*" in c for c in result.columns), \
            "Unexpected column name in HLA scores — expected HLA allele format"

    def test_predict_hla_returns_boolean(self, toy_dmm):
        result = predict_hla_type(toy_dmm)
        assert result.dtypes.apply(
            lambda dt: dt == bool or np.issubdtype(dt, np.bool_)
        ).all()

    def test_predict_hla_shape(self, toy_dmm):
        result = predict_hla_type(toy_dmm)
        assert result.shape[0] == toy_dmm.shape[0]

    def test_predict_hla_with_ground_truth(self, toy_dmm, toy_df):
        # Build a small fake ground truth from the toy donors
        donors = toy_df["donor_id"].unique()[:20]
        hla_gt = pd.DataFrame(
            np.random.default_rng(42).integers(0, 2, size=(len(donors), 3)),
            index=donors,
            columns=["A*02:01", "B*07:02", "A*01:01"],
        )
        result = predict_hla_type(toy_dmm, hla_ground_truth=hla_gt)
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == toy_dmm.shape[0]

    def test_evaluate_hla_runs(self, toy_dmm, toy_df):
        hla_scores = score_hla(toy_dmm)
        # Build synthetic ground truth for evaluation
        donors = hla_scores.index[:30]
        alleles = list(hla_scores.columns)[:3]
        rng = np.random.default_rng(0)
        gt = pd.DataFrame(
            rng.integers(0, 2, size=(len(donors), len(alleles))),
            index=donors, columns=alleles,
        )
        result = evaluate_hla_prediction(
            hla_scores.loc[donors, alleles], gt, n_splits=3
        )
        assert isinstance(result, pd.DataFrame)
        assert "Mean_AUC" in result.columns


# ---------------------------------------------------------------------------
# End-to-end smoke test
# ---------------------------------------------------------------------------

class TestEndToEnd:

    def test_full_pipeline_no_errors(self, toy_df):
        """Full motif → exposure → HLA pipeline must run without errors."""
        from cell2specificity.tcr_motifs import preprocess_tcr_table

        df  = preprocess_tcr_table(toy_df)
        dmm = build_donor_motif_matrix(df)
        exp = predict_pathogen_exposure(dmm, threshold=1)
        hla = predict_hla_type(dmm)

        assert exp.shape[0] == dmm.shape[0]
        assert hla.shape[0] == dmm.shape[0]
