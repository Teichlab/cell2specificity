"""
specificity
-----------
Clinical metadata inference from TCR motif composition.

This module enables two key inference tasks on new donor TCR repertoires:

**Pathogen exposure prediction**
    Given a donor's set of TCR motifs, count how many match motifs
    associated with each pathogen in the reference atlas. A "double-hit"
    threshold (≥ 2 motifs) provides high-precision calls validated at
    p = 8×10⁻¹⁷ for SARS-CoV-2.

**HLA genotype inference**
    MHC-I-restricted TCR motifs serve as proxies for HLA alleles. For each
    allele, a Youden's-J-optimal threshold learned from donors with known
    HLA enables imputation across donors lacking genetic data.

Bundled reference data
~~~~~~~~~~~~~~~~~~~~~~
Two CSV tables are shipped with the package (place in
``src/cell2specificity/specificity/data/``):

* ``disease_associated_motifs_hla.csv`` — motif → predicted pathogen
* ``df_motifs_with_hla.csv``            — motif → MHC-I restricted HLA allele

Typical workflow
~~~~~~~~~~~~~~~~
::

    from cell2specificity.tcr_motifs import (
        preprocess_tcr_table, to_matching_frame,
        build_chain_index, query_to_ref,
    )
    from cell2specificity.specificity import (
        build_donor_motif_matrix,
        predict_pathogen_exposure,
        predict_hla_type,
    )

    # 1. Map new user VDJ data onto atlas motifs
    df_ref  = pd.read_csv('atlas_clone_df_with_motifs.csv')
    df_user = preprocess_tcr_table(pd.read_csv('my_tcr_data.csv'))
    df_ref_m  = to_matching_frame(df_ref)
    df_user_m = to_matching_frame(df_user)

    idx_alpha = build_chain_index(df_ref_m, chain='vj',  mm=1)
    idx_beta  = build_chain_index(df_ref_m, chain='vdj', mm=1)
    df_user_m = query_to_ref(df_ref_m, df_user_m, idx_alpha, idx_beta)

    # 2. Build presence/absence matrix
    dmm = build_donor_motif_matrix(df_user_m)

    # 3. Predict
    exposure  = predict_pathogen_exposure(dmm, threshold=2)
    hla_calls = predict_hla_type(dmm)
"""

from ._predict import (
    build_donor_motif_matrix,
    map_motifs_to_metadata,
    score_pathogen_exposure,
    predict_pathogen_exposure,
    score_hla,
    predict_hla_type,
    evaluate_hla_prediction,
)
from ._data import load_pathogen_motifs, load_hla_motifs

__all__ = [
    # core mapping
    "build_donor_motif_matrix",
    "map_motifs_to_metadata",
    # pathogen
    "score_pathogen_exposure",
    "predict_pathogen_exposure",
    # HLA
    "score_hla",
    "predict_hla_type",
    "evaluate_hla_prediction",
    # data loaders
    "load_pathogen_motifs",
    "load_hla_motifs",
]
