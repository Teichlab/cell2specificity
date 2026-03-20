"""
epitope
-------
Epitope prediction and genome-wide peptide scanning.

Uses NetMHCpan to nominate candidate peptides from pathogen proteomes that
are predicted to bind HLA alleles enriched in a given TCR motif's donor set.
Generates ranked TCR-peptide-HLA candidate pairs for downstream structural
modelling.

Key functions (to be implemented)
----------------------------------
- scan_proteome()           : Enumerate candidate peptides from a FASTA proteome
- run_netmhcpan()           : Wrapper for NetMHCpan binding predictions
- filter_binders()          : Apply affinity / rank thresholds
- rank_candidates()         : Prioritise peptides per motif / HLA allele
"""
