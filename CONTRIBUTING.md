# Contributing to cell2specificity

Thanks for contributing. This document covers the conventions for adding new modules, submitting code, and working with collaborators.

---

## Repository layout

```
cell2specificity/
├── src/cell2specificity/   # All source code lives here (src layout)
│   └── <module>/
│       ├── __init__.py     # Module docstring + public API
│       └── *.py
├── tests/                  # Mirrors src layout; one test file per module
├── docs/                   # Sphinx documentation
├── pyproject.toml          # Single source of truth for build + deps
└── CONTRIBUTING.md
```

## Adding a new module

1. Create `src/cell2specificity/<your_module>/` with an `__init__.py`.
2. Export public functions from `__init__.py`; keep internals in private submodules (`_helpers.py`, etc.).
3. Add your module to the imports in `src/cell2specificity/__init__.py`.
4. Add any new dependencies to `pyproject.toml` under `[project.dependencies]`.
5. Write tests in `tests/test_<your_module>.py`.

## Code style

- Formatting and linting: `ruff` (config in `pyproject.toml`).
- Type hints on all public functions.
- NumPy-style docstrings.

Run before committing:
```bash
ruff check src/ tests/
ruff format src/ tests/
pytest
```

## Branching and PRs

- Work on a feature branch: `git checkout -b feat/<module-name>`
- Keep PRs focused — one module or fix per PR.
- All PRs require at least one review before merge to `main`.
- `main` is the stable branch; do not push directly.

## Commit messages

Use the conventional commits style:
```
feat(tcr_motifs): add public motif filtering function
fix(hla): handle missing class II alleles gracefully
docs(structural): update classifier feature description
```
