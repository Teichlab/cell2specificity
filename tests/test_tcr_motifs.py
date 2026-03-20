"""
Tests for cell2specificity.tcr_motifs

Run with:  pytest tests/test_tcr_motifs.py -v

Toy dataset (tests/data/toy_tcr_atlas.csv):
  865 cells | 113 donors | 11 pathogens
  50 MAIT | 20 iNKT | 795 Conventional
"""

import pandas as pd
import pytest
from pathlib import Path

from cell2specificity.tcr_motifs import (
    VDJ_COLS,
    preprocess_tcr_table,
    annotate_invariant,
    invariant_summary,
    classify_invariant_row,
    to_matching_frame,
    build_chain_index,
    query_to_ref,
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
def toy_preprocessed(toy_df):
    return preprocess_tcr_table(toy_df)


@pytest.fixture
def minimal_df():
    """Three-row DataFrame covering MAIT, iNKT, and Conventional."""
    return pd.DataFrame({
        "IR_VDJ_1_v_call":      ["TRBV4-2",    "TRBV20-1",   "TRBV6-2"],
        "IR_VDJ_1_j_call":      ["TRBJ2-3",    "TRBJ2-7",    "TRBJ1-2"],
        "IR_VDJ_1_junction_aa": ["CASSLGNTIYF","CSARDLGTEAFF","CASSPSGNTIYF"],
        "IR_VJ_1_v_call":       ["TRAV1-2*01", "TRAV10*01",  "TRAV12-2"],
        "IR_VJ_1_j_call":       ["TRAJ33*01",  "TRAJ18*01",  "TRAJ9"],
        "IR_VJ_1_junction_aa":  ["CAVMDSNYQLIW","CVVGDRGSYLNQLIW","CATSDA"],
        "donor_id":             ["d1",          "d1",          "d2"],
        "motif":                [101,            102,           103],
    })


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------

class TestPreprocess:

    def test_returns_copy(self, minimal_df):
        result = preprocess_tcr_table(minimal_df)
        assert result is not minimal_df

    def test_all_vdj_cols_present(self, minimal_df):
        result = preprocess_tcr_table(minimal_df)
        for col in VDJ_COLS:
            assert col in result.columns

    def test_allele_suffix_added_to_bare_gene(self, minimal_df):
        # TRBV4-2 has no * — should become TRBV4-2*01
        result = preprocess_tcr_table(minimal_df, add_alleles=True)
        assert all("*" in v for v in result["IR_VDJ_1_v_call"].dropna())

    def test_existing_allele_not_doubled(self, minimal_df):
        # TRAV1-2*01 already has *01 — must stay as-is
        result = preprocess_tcr_table(minimal_df, add_alleles=True)
        assert result.loc[0, "IR_VJ_1_v_call"] == "TRAV1-2*01"

    def test_missing_cols_raises(self):
        with pytest.raises(ValueError, match="missing required VDJ columns"):
            preprocess_tcr_table(pd.DataFrame({"a": [1]}))

    def test_toy_data_survives_preprocess(self, toy_df):
        result = preprocess_tcr_table(toy_df)
        assert result.shape[0] == 865
        assert all(c in result.columns for c in VDJ_COLS)


# ---------------------------------------------------------------------------
# Invariant annotation — unit level
# ---------------------------------------------------------------------------

class TestClassifyInvariantRow:

    def test_mait_traj33(self):
        assert classify_invariant_row("TRAV1-2*01", "TRAJ33*01") == "MAIT"

    def test_mait_traj20(self):
        assert classify_invariant_row("TRAV1-2*01", "TRAJ20*01") == "MAIT"

    def test_mait_traj12(self):
        assert classify_invariant_row("TRAV1-2*01", "TRAJ12*01") == "MAIT"

    def test_inkt(self):
        assert classify_invariant_row("TRAV10*01", "TRAJ18*01") == "iNKT"

    def test_conventional(self):
        assert classify_invariant_row("TRAV12-2*01", "TRAJ9*01") == "Conventional"

    def test_trav1_2_wrong_j_is_conventional(self):
        # TRAV1-2 but wrong J → not MAIT
        assert classify_invariant_row("TRAV1-2*01", "TRAJ9*01") == "Conventional"

    def test_trav10_wrong_j_is_conventional(self):
        assert classify_invariant_row("TRAV10*01", "TRAJ9*01") == "Conventional"

    def test_nan_v_is_conventional(self):
        assert classify_invariant_row(float("nan"), "TRAJ33") == "Conventional"

    def test_none_v_is_conventional(self):
        assert classify_invariant_row(None, "TRAJ33") == "Conventional"


# ---------------------------------------------------------------------------
# Invariant annotation — DataFrame level
# ---------------------------------------------------------------------------

class TestAnnotateInvariant:

    def test_adds_column(self, minimal_df):
        result = annotate_invariant(minimal_df)
        assert "invariant" in result.columns

    def test_correct_labels(self, minimal_df):
        result = annotate_invariant(minimal_df)
        assert result.loc[0, "invariant"] == "MAIT"
        assert result.loc[1, "invariant"] == "iNKT"
        assert result.loc[2, "invariant"] == "Conventional"

    def test_returns_copy(self, minimal_df):
        annotate_invariant(minimal_df)
        assert "invariant" not in minimal_df.columns

    def test_custom_out_col(self, minimal_df):
        result = annotate_invariant(minimal_df, out_col="tcr_type")
        assert "tcr_type" in result.columns

    def test_missing_column_raises(self, minimal_df):
        with pytest.raises(ValueError, match="not found"):
            annotate_invariant(minimal_df, vj_v_col="nonexistent")

    # --- toy data integration ---

    def test_toy_mait_count(self, toy_preprocessed):
        result = annotate_invariant(toy_preprocessed)
        assert result["invariant"].eq("MAIT").sum() == 50

    def test_toy_inkt_count(self, toy_preprocessed):
        result = annotate_invariant(toy_preprocessed)
        assert result["invariant"].eq("iNKT").sum() == 20

    def test_toy_conventional_count(self, toy_preprocessed):
        result = annotate_invariant(toy_preprocessed)
        assert result["invariant"].eq("Conventional").sum() == 795

    def test_toy_total_unchanged(self, toy_preprocessed):
        result = annotate_invariant(toy_preprocessed)
        assert len(result) == 865


# ---------------------------------------------------------------------------
# Invariant summary
# ---------------------------------------------------------------------------

class TestInvariantSummary:

    def test_proportions_sum_to_one(self, minimal_df):
        import numpy as np
        df = annotate_invariant(minimal_df)
        summary = invariant_summary(df, normalize=True)
        assert np.isclose(summary.values.sum(), 1.0)

    def test_raw_counts_correct(self, minimal_df):
        df = annotate_invariant(minimal_df)
        summary = invariant_summary(df, normalize=False)
        assert summary["MAIT"].sum() == 1
        assert summary["iNKT"].sum() == 1
        assert summary["Conventional"].sum() == 1

    def test_groupby_pathogen(self, toy_preprocessed):
        df = annotate_invariant(toy_preprocessed)
        summary = invariant_summary(df, groupby="pathogen", normalize=True)
        assert "SARS-CoV-2" in summary.index
        # proportions per group must sum to 1
        import numpy as np
        assert np.allclose(summary.sum(axis=1), 1.0)

    def test_missing_invariant_col_raises(self, toy_preprocessed):
        with pytest.raises(ValueError, match="Run annotate_invariant"):
            invariant_summary(toy_preprocessed)


# ---------------------------------------------------------------------------
# TCR matching
# ---------------------------------------------------------------------------

class TestMatching:

    def test_self_match_zero_mm(self, minimal_df):
        df_m = to_matching_frame(preprocess_tcr_table(minimal_df))
        idx_b = build_chain_index(df_m, chain="vdj", mm=0)
        idx_a = build_chain_index(df_m, chain="vj",  mm=0)
        result = query_to_ref(df_m, df_m.copy(), idx_a, idx_b)
        # every cell must find itself
        assert (result["n_ref_hits"] >= 1).all()

    def test_result_annotation_columns(self, minimal_df):
        df_m = to_matching_frame(preprocess_tcr_table(minimal_df))
        idx_b = build_chain_index(df_m, chain="vdj", mm=0)
        idx_a = build_chain_index(df_m, chain="vj",  mm=0)
        result = query_to_ref(df_m, df_m.copy(), idx_a, idx_b)
        for col in ["ref_hits", "n_ref_hits", "ref_motifs", "n_ref_motifs"]:
            assert col in result.columns

    def test_ref_motifs_are_lists(self, minimal_df):
        df_m = to_matching_frame(preprocess_tcr_table(minimal_df))
        idx_b = build_chain_index(df_m, chain="vdj", mm=0)
        idx_a = build_chain_index(df_m, chain="vj",  mm=0)
        result = query_to_ref(df_m, df_m.copy(), idx_a, idx_b)
        assert all(isinstance(m, list) for m in result["ref_motifs"])

    def test_no_false_match_different_v(self, minimal_df):
        df_m = to_matching_frame(preprocess_tcr_table(minimal_df))
        idx_b = build_chain_index(df_m, chain="vdj", mm=0)
        idx_a = build_chain_index(df_m, chain="vj",  mm=0)
        # query with a completely different V gene → should get 0 hits
        query = df_m.copy()
        query["vj_v"] = "TRAV99-1"
        query["vdj_v"] = "TRBV99-1"
        result = query_to_ref(df_m, query, idx_a, idx_b)
        assert (result["n_ref_hits"] == 0).all()

    def test_toy_data_matching(self, toy_preprocessed):
        df_m  = to_matching_frame(toy_preprocessed)
        ref   = df_m.dropna(subset=["motif"]).head(400)
        query = df_m.dropna(subset=["motif"]).tail(100)
        idx_b = build_chain_index(ref, chain="vdj", mm=1)
        idx_a = build_chain_index(ref, chain="vj",  mm=1)
        result = query_to_ref(ref, query.copy(), idx_a, idx_b)
        # not every query will match (different donors / pathogens), but shape must hold
        assert len(result) == 100
        assert result["n_ref_hits"].dtype.kind in ("i", "u")  # integer type
