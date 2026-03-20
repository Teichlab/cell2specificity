"""
annotation
----------
Cell type annotation of T cells using CellTypist.

Thin wrapper around the `celltypist` package that resolves the `model`
argument to either:
  - a built-in CellTypist model (any string accepted by celltypist itself, e.g.
    ``'Immune_All_Low.pkl'``), or
  - one of the three models shipped with this package, referenced by their
    short aliases:

    =====================================  ===================================
    Alias                                  Description
    =====================================  ===================================
    ``'paninfection_level2'``              Pan-infection atlas, annotation
                                           level 2 (broad lineages)
    ``'paninfection_CD4_level3'``          CD4 T cell fine-grained subtypes
    ``'paninfection_CD8_level3'``          CD8 T cell fine-grained subtypes
    =====================================  ===================================

  - an explicit path to any ``.pkl`` CellTypist model file on disk.

All other keyword arguments are forwarded directly to
``celltypist.annotate()``, so the full CellTypist API remains accessible.

Examples
--------
>>> from cell2specificity.annotation import annotate
>>> predictions = annotate(adata, model='paninfection_level2',
...                        majority_voting=True)
>>> predictions = annotate(adata, model='Immune_All_Low.pkl')
"""

from __future__ import annotations

import importlib.resources
from pathlib import Path
from typing import Union

import anndata
import celltypist
try:
    from celltypist import AnnotationResult
except ImportError:
    # celltypist >= 1.6 exposes the result class under a different path
    from celltypist.classifier import AnnotationResult  # type: ignore[no-redef]

__all__ = ["annotate", "list_bundled_models", "BUNDLED_MODELS"]

# ---------------------------------------------------------------------------
# Bundled model registry
# ---------------------------------------------------------------------------

#: Short alias → filename for models shipped with the package.
BUNDLED_MODELS: dict[str, str] = {
    "paninfection_level2":    "paninfection_annotation_level_2.pkl",
    "paninfection_CD4_level3": "paninfection_CD4_annotation_level_3.pkl",
    "paninfection_CD8_level3": "paninfection_CD8_annotation_level_3.pkl",
}

_MODELS_DIR = Path(__file__).parent / "models"


def _resolve_model(model: Union[str, Path, "celltypist.models.Model"]):
    """
    Return a celltypist-compatible model argument.

    Priority order:
      1. Already a ``celltypist.models.Model`` instance — pass through.
      2. Known bundled alias — load from the package ``models/`` directory.
      3. Existing file path — load directly.
      4. Anything else — pass through to celltypist (built-in model name).
    """
    if isinstance(model, celltypist.models.Model):
        return model

    if isinstance(model, str) and model in BUNDLED_MODELS:
        pkl_path = _MODELS_DIR / BUNDLED_MODELS[model]
        if not pkl_path.exists():
            raise FileNotFoundError(
                f"Bundled model '{model}' not found at {pkl_path}.\n"
                "Make sure the model .pkl files have been added to "
                "src/cell2specificity/annotation/models/."
            )
        return celltypist.models.Model.load(str(pkl_path))

    path = Path(model) if isinstance(model, str) else model
    if path.exists():
        return celltypist.models.Model.load(str(path))

    # Fall through to celltypist's own model resolution (built-in names)
    return model


def annotate(
    input_data: Union[anndata.AnnData, "str", "Path"],
    model: Union[str, Path, "celltypist.models.Model"] = "paninfection_level2",
    *,
    majority_voting: bool = False,
    mode: str = "best match",
    p_thres: float = 0.5,
    transpose_input: bool = False,
    over_clustering: Union[str, None] = None,
    **kwargs,
) -> AnnotationResult:
    """
    Predict T cell types using CellTypist.

    Parameters
    ----------
    input_data
        An AnnData object (cells × genes, normalised log1p counts), or a path
        to a ``.h5ad`` file. Gene expression values must be log1p-normalised
        counts in ``adata.X`` (or ``adata.raw.X``).
    model
        Which model to use. Accepts:

        * A bundled alias: ``'paninfection_level2'``,
          ``'paninfection_CD4_level3'``, ``'paninfection_CD8_level3'``.
        * A path to any ``.pkl`` CellTypist model file.
        * Any built-in CellTypist model name (e.g. ``'Immune_All_Low.pkl'``).
        * A pre-loaded ``celltypist.models.Model`` instance.

    majority_voting
        Whether to refine predictions by majority voting within local
        over-clustering groups. Recommended when ``over_clustering`` is set.
    mode
        ``'best match'`` (default) or ``'prob match'``. See CellTypist docs.
    p_thres
        Probability threshold for ``'prob match'`` mode.
    transpose_input
        Set ``True`` if input is genes × cells.
    over_clustering
        ``adata.obs`` key for pre-computed cluster labels to use for majority
        voting. Ignored when ``majority_voting=False``.
    **kwargs
        Any additional keyword arguments forwarded to ``celltypist.annotate()``.

    Returns
    -------
    celltypist.AnnotationResult
        Contains per-cell predicted labels, decision-function scores, and
        (optionally) majority-voted labels. Call ``.to_adata()`` to embed
        results into the AnnData object.

    Examples
    --------
    Annotate with the pan-infection level-2 model:

    >>> predictions = annotate(adata, model='paninfection_level2',
    ...                        majority_voting=True)
    >>> adata = predictions.to_adata()

    Annotate CD8 T cells with the fine-grained model:

    >>> predictions = annotate(adata_cd8, model='paninfection_CD8_level3')

    Use a standard CellTypist built-in model:

    >>> predictions = annotate(adata, model='Immune_All_Low.pkl')
    """
    resolved = _resolve_model(model)
    return celltypist.annotate(
        input_data,
        model=resolved,
        majority_voting=majority_voting,
        mode=mode,
        p_thres=p_thres,
        transpose_input=transpose_input,
        over_clustering=over_clustering,
        **kwargs,
    )


def list_bundled_models() -> dict[str, str]:
    """
    Return the registry of models bundled with this package.

    Returns
    -------
    dict
        ``{alias: filename}`` for each bundled model.
    """
    return dict(BUNDLED_MODELS)
