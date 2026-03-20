"""
hla
---
HLA genotype inference from raw single-cell RNA-sequencing data.

Infers six-digit class I (HLA-A, -B, -C) and class II (HLA-DPA1/DPB1,
HLA-DQA1/DQB1, HLA-DRB1) haplotypes directly from scRNA-seq reads.
Supports downstream per-motif HLA restriction analysis and infection-history
reconstruction.

Key functions (to be implemented)
----------------------------------
- infer_hla_from_scrna()       : Run HLA typing from BAM / fastq input
- load_hla_table()             : Parse and validate a donor HLA table
- quantify_hla_restriction()  : Compute HLA sharing within TCR motifs (chi-squared)
- bootstrap_hla_overlap()     : Compare within-motif vs random donor HLA overlap
"""
