"""
cell2specificity
================
A computational toolkit for systematic antigen-specificity inference from
single-cell TCR and transcriptomic data.

Integrates TCR clonotype analysis, HLA genotype inference, epitope prediction,
and structural modelling of TCR-peptide-HLA binding.

Reference
---------
Dratva et al. (2026) "Single-cell analysis of human T cells across infections
unlocks systematic antigen-specificity inference."
"""

__version__ = "0.1.0"
__author__ = "Lisa M Dratva"

from cell2specificity import tcr_motifs, hla, epitope, structural, specificity, annotation, utils

__all__ = [
    "tcr_motifs",
    "hla",
    "epitope",
    "structural",
    "specificity",
    "annotation",
    "utils",
]
