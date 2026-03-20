"""
structural
----------
Structure-informed classification of TCR–peptide–HLA (TCR-pHLA) binding.

.. note::
   **This module is reserved for a collaborator contribution.**
   The public API skeleton below documents the intended interface so that
   downstream code can import stubs; implementations should be added as
   submodules and exposed here.

Overview
~~~~~~~~
The pipeline described in the manuscript (Dratva et al. 2026) proceeds in
three stages:

1. **Structural modelling** — Submit candidate TCR-pHLA complexes to TCRdock,
   which returns predicted structures with PAE and pLDDT confidence scores.

2. **Feature extraction** — Parse TCRdock output to compute per-complex
   structural features: overall PAE, per-complex PAE, per-complex pLDDT,
   and contact matrices between TCR-beta/alpha chains and the peptide/HLA.

3. **Binding classification** — A random forest classifier trained on
   BEAM-validated binders and in-silico controls scores each complex.
   The final model (selected from 1,000+ parallel random forests via a
   decision-tree meta-selector) achieves mean AUC 0.910, accuracy 0.856.

Intended public API
~~~~~~~~~~~~~~~~~~~
::

    from cell2specificity.structural import (
        run_tcrdock,
        extract_structural_features,
        predict_binding,
        train_rf_classifier,   # optional: re-train on new labelled data
    )

    candidates = run_tcrdock(tcr_peptide_hla_table)
    features   = extract_structural_features(candidates)
    scores     = predict_binding(features)            # uses bundled classifier
    hits       = scores[scores['binding_score'] >= 0.5]

Notes for the contributor
~~~~~~~~~~~~~~~~~~~~~~~~~
* Add submodules (e.g. ``_tcrdock.py``, ``_features.py``, ``_classifier.py``)
  and expose their public symbols in this ``__init__.py``.
* The bundled random forest classifier (``.joblib`` or ``.pkl``) should be
  placed in ``src/cell2specificity/structural/models/`` and loaded lazily.
* External dependencies (TCRdock, any structural libraries) should be listed
  in ``pyproject.toml`` under ``[project.optional-dependencies]`` with a
  suitable extras key, e.g. ``pip install cell2specificity[structural]``.
* Add tests in ``tests/test_structural.py`` following the conventions in
  ``CONTRIBUTING.md``.
"""

# Stub — to be implemented by contributor.
# Remove this block and replace with real imports once modules are added.

__all__: list = []


def run_tcrdock(*args, **kwargs):
    """Submit TCR-pHLA complexes to TCRdock for structural modelling. **(Not yet implemented.)**"""
    raise NotImplementedError("structural.run_tcrdock is not yet implemented.")


def extract_structural_features(*args, **kwargs):
    """Parse TCRdock output into a feature matrix. **(Not yet implemented.)**"""
    raise NotImplementedError("structural.extract_structural_features is not yet implemented.")


def predict_binding(*args, **kwargs):
    """Score candidate complexes with the random forest classifier. **(Not yet implemented.)**"""
    raise NotImplementedError("structural.predict_binding is not yet implemented.")


def train_rf_classifier(*args, **kwargs):
    """(Re-)train the random forest on new labelled TCR-pHLA data. **(Not yet implemented.)**"""
    raise NotImplementedError("structural.train_rf_classifier is not yet implemented.")
