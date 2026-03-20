# Tutorial: cell2specificity

This tutorial walks through the three core workflows using the test dataset
shipped with the package under `tests/data/`, except for the initial cell state annotation, where you should substitute your own expression dataset.


## Installation

```bash
git clone https://github.com/lisadratva/cell2specificity.git
cd cell2specificity
pip install -e ".[dev]"
```

For TCR motif inference (requires tcrdist3):
```bash
pip install -e ".[motifs]"
```


## 0. Cell type annotation with CellTypist

```python
import scanpy as sc
from cell2specificity.annotation import annotate, list_bundled_models

# See available bundled models
print(list_bundled_models())
# {
#   'paninfection_level2':     'paninfection_annotation_level_2.pkl',
#   'paninfection_CD4_level3': 'paninfection_CD4_annotation_level_3.pkl',
#   'paninfection_CD8_level3': 'paninfection_CD8_annotation_level_3.pkl',
# }

# Load an AnnData object (log1p-normalised counts in adata.X)
adata = sc.read_h5ad("your_data.h5ad")

# Annotate with the pan-infection level-2 model
predictions = annotate(adata, model="paninfection_level2")

# For CD8 T cells only:
predictions_cd8 = annotate(adata, model="paninfection_CD8_level3",
                            majority_voting=True)

# Embed results back into AnnData
adata = predictions.to_adata()
print(adata.obs[["predicted_labels", "conf_score"]].head())
```

You can also use any built-in CellTypist model by name:
```python
predictions = annotate(adata, model="Immune_All_Low.pkl")
```

## 1. Loading test data

The package ships a small subset of the pan-infection T cell atlas for
testing and exploration. It contains ~500 cells with paired scRNA-TCR
sequences across several donors.

```python
import pandas as pd

df = pd.read_csv("tests/data/toy_tcr_atlas.csv", index_col=0)
print(df.shape)        # e.g. (500, 20)
print(df.columns.tolist())
```

Expected columns include the six canonical VDJ fields:
`IR_VDJ_1_v_call`, `IR_VDJ_1_j_call`, `IR_VDJ_1_junction_aa`,
`IR_VJ_1_v_call`, `IR_VJ_1_j_call`, `IR_VJ_1_junction_aa`,
plus `donor_id`, `motif`, `pathogen`, and cell annotation columns.

## 2. TCR preprocessing and invariant annotation

```python
from cell2specificity.tcr_motifs import (
    preprocess_tcr_table,
    annotate_invariant,
    invariant_summary,
)

# Standardise V/J gene allele suffixes
df = preprocess_tcr_table(df)

# Annotate MAIT, iNKT, and Conventional T cells
df = annotate_invariant(df)
print(df["invariant"].value_counts())
# Conventional    ...
# MAIT            ...
# iNKT            ...

# Summarise by pathogen
summary = invariant_summary(df, groupby="pathogen", normalize=True)
print(summary)
```

## 3. Fast TCR matching to atlas motifs

Map new (query) TCR sequences to the reference atlas and retrieve the
motif IDs of the closest matching clones.

```python
from cell2specificity.tcr_motifs import (
    to_matching_frame,
    build_chain_index,
    query_to_ref,
)

# Convert column names to the compact vdj_* / vj_* format
df_m = to_matching_frame(df)

# Split into a reference set (cells with known motifs) and a query set
df_ref   = df_m.dropna(subset=["motif"])
df_query = df_m.sample(50, random_state=42)

# Build indices once — reuse for many queries
idx_beta  = build_chain_index(df_ref, chain="vdj", mm=1)
idx_alpha = build_chain_index(df_ref, chain="vj",  mm=1)

# Map query TCRs to reference and annotate with motif IDs
df_query = query_to_ref(df_ref, df_query, idx_alpha, idx_beta)

print(df_query[["n_ref_hits", "n_ref_motifs", "ref_motifs"]].head())
```

## 4. Pathogen exposure inference

```python
from cell2specificity.motif_based_inference import (
    build_donor_motif_matrix,
    predict_pathogen_exposure,
)

# Binary donor × motif presence matrix
dmm = build_donor_motif_matrix(df)

# Predict exposure using the bundled pathogen-motif reference
# threshold=2 → "double-hit" rule (high-precision call)
exposure = predict_pathogen_exposure(dmm, threshold=2)
print(exposure)
# rows = donors, columns = pathogens, values = True/False
```


## 5. HLA genotype inference

```python
from cell2specificity.motif_based_inference import (
    score_hla,
    predict_hla_type,
)

# Score donors by HLA-restricted motif counts
hla_scores = score_hla(dmm)
print(hla_scores.head())

# Predict HLA type — no ground truth needed; uses default threshold ≥ 1
hla_calls = predict_hla_type(dmm)
print(hla_calls.head())

# With ground-truth HLA for threshold learning:
# hla_gt = pd.read_csv("hla_ground_truth.csv", index_col=0)
# hla_calls = predict_hla_type(dmm, hla_ground_truth=hla_gt)
```


## 7. Running the test suite

```bash
pytest tests/ -v
```

All tests use the test data in `tests/data/` and require no external data
or compute resources.
