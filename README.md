# cell2specificity

A computational toolkit for systematic antigen-specificity inference from single-cell TCR and transcriptomic data.

**cell2specificity** integrates TCR clonotype analysis, HLA genotype inference, epitope prediction, and structure-informed modelling of TCR-peptide-HLA binding. It is the companion software to:

> Dratva et al. (2026) *Single-cell analysis of human T cells across infections unlocks systematic antigen-specificity inference.*

---

## Overview

Starting from scRNA+TCR-seq data, the toolkit enables:

- **Invariant T cell annotation** — classify MAIT and iNKT cells from V/J gene usage
- **TCR motif discovery** — group clonotypes into shared-specificity clusters using Cell2TCR
- **Fast TCR matching** — map new repertoires onto atlas motifs in <1 s via indexed seed-and-extend
- **Pathogen exposure inference** — predict donor infection history from TCR motif composition
- **HLA genotype inference** — impute HLA alleles from MHC-I-restricted public TCR motifs
- **Cell type annotation** — predict T cell states using bundled CellTypist models trained on the atlas
- **Structural modelling** — run TCRdock and classify TCR-pHLA binding (AUC 0.910) *(collaborator module)*

---

## Installation

```bash
git clone https://github.com/needle-bio/cell2specificity.git
cd cell2specificity
pip install -e ".[dev]"
```

For TCR motif inference (requires tcrdist3):
```bash
pip install -e ".[motifs]"
```

Python ≥ 3.10 required.

---

## Quickstart

```python
import pandas as pd
from cell2specificity.tcr_motifs import preprocess_tcr_table, annotate_invariant
from cell2specificity.specificity import build_donor_motif_matrix, predict_pathogen_exposure, predict_hla_type
from cell2specificity.annotation import annotate

# 1. Preprocess VDJ table and annotate invariant T cells
df = preprocess_tcr_table(pd.read_csv("my_tcr_data.csv"))
df = annotate_invariant(df)

# 2. Build donor × motif matrix and run inference
dmm      = build_donor_motif_matrix(df)
exposure = predict_pathogen_exposure(dmm, threshold=2)  # double-hit rule
hla      = predict_hla_type(dmm)

# 3. Annotate cell states with bundled CellTypist models
predictions = annotate(adata, model="paninfection_level2")
adata = predictions.to_adata()
```

**→ See the full step-by-step walkthrough in [docs/tutorial.md](docs/tutorial.md)**

---

## Bundled models and reference data

Three CellTypist models trained on the pan-infection T cell atlas are shipped
with the package:

| Alias | Description |
|---|---|
| `paninfection_level2` | Broad T cell lineages (CD4, CD8, MAIT, iNKT, γδ, NKT) |
| `paninfection_CD4_level3` | Fine-grained CD4 T cell subtypes (29 states) |
| `paninfection_CD8_level3` | Fine-grained CD8 T cell subtypes (12 states) |

Two reference tables for clinical inference are bundled under
`src/cell2specificity/specificity/data/`:

- `disease_associated_motifs_hla.csv` — motif → predicted pathogen
- `df_motifs_with_hla.csv` — motif → MHC-I restricted HLA allele + metadata

---

## Package structure

```
src/cell2specificity/
├── tcr_motifs/     # Preprocessing, invariant annotation, seed-and-extend
│   │               # matching, Cell2TCR motif inference, VDJdb queries
│   ├── _preprocess.py
│   ├── _invariant.py
│   └── _matching.py
├── specificity/    # Pathogen exposure + HLA inference from motif composition
│   ├── _predict.py
│   └── data/       # Bundled reference CSVs
├── annotation/     # CellTypist wrapper + bundled pan-infection models
│   └── models/     # .pkl model files
├── structural/     # TCRdock + random forest binding classifier [collaborator]
├── hla/            # HLA genotype inference from raw scRNA-seq reads
├── epitope/        # NetMHCpan wrapper and genome-wide peptide scanning
└── utils/          # Shared helpers
tests/
├── data/           # Toy atlas subset for self-contained testing
└── test_*.py
docs/
└── tutorial.md     # Full worked tutorial
```

---

## Running tests

```bash
pytest tests/ -v
```

All tests run against the toy dataset in `tests/data/` — no external data or
compute required.

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). New modules follow the `src/` layout
and must include tests and docstrings. The `structural/` module is reserved for
a collaborator — see its docstring for the intended API.

---

## Citation

If you use this toolkit, please cite:

> Dratva LM et al. (2026) *Single-cell analysis of human T cells across infections unlocks systematic antigen-specificity inference.*

## License

Apache 2.0 — see [LICENSE](LICENSE).
Copyright 2026 Lisa M Dratva (lmd76@cam.ac.uk)
