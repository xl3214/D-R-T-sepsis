# D-R-T-sepsis
## Multi-omic Interrogation of Host Resistance, Tolerance, and Damage in Sepsis and Severe Infection Reveals Conserved Immune Programs with Pathogen-Dependent Heterogeneity 

# D/R/T domain scores: analysis code

Code for the Damage / Resistance / Tolerance (D/R/T) domain-score analyses. All scripts are run from the project root (the folder containing `data/`), in the order below.

| Step | Script | What it does | Main outputs |
|---|---|---|---|
| 1 | `01_DRT_scores.R` | Derives the multi-omic PC1 weights and the transcriptomic projection models in RESERVE-U-2-TOR (`DERIVE = TRUE`), checks them against the published values, and scores every cohort. Also runs the TOR-axis checks in RESERVE-U-1-EBB and the projection performance (TOR out-of-fold, EBB external). | `DRT_scores_output/<cohort>/`, `derivation/`, `transfer_checks/` |
| 2 | `02_Figure1_multiomic_scores.R` | Figure 1 (panels A–E) from the step 1 outputs. | `DRT_scores_output/figures/` |
| 3 | `clinical/03_RESERVE_U_2_TOR_DRT_Clinical.R` | RESERVE-U-2-TOR clinical analyses (joint 30-day mortality and MEWS models, Firth sensitivity, infection → score regressions, projection correlations and supportive models). | `DRT_scores_output/Clinical analyses/RESERVE_U_2_TOR/` |
| 3a | `clinical/03a_RESERVE_U_2_TOR_Table_Characteristics.R` | RESERVE-U-2-TOR participant characteristics table. | same folder (.docx) |
| 4 | `clinical/04_RESERVE_U_1_EBB_DRT_Clinical.R` | Same analyses in RESERVE-U-1-EBB. | `.../RESERVE_U_1_EBB/` |
| 5 | `clinical/05_RESERVE_U_Pooled_DRT_Clinical.R` | Post hoc pooled RESERVE-U analyses, infection × domain interaction models, exploratory urine LAM analysis. | `.../Pooled/` |
| 5a | `clinical/05a_RESERVE_U_Figure3.R` | Figure 3 (source at the end of step 5, or run after it). | `.../Pooled/` |
| 6 | `clinical/06_ACESO_Ghana_DRT_Clinical.R` | ACESO-Ghana: sepsis vs healthy, joint MEWS and 28-day mortality models, longitudinal mixed models, 72-hour landmark analysis. | `.../Ghana/` |
| 7 | `clinical/07_GSE185263_DRT_Clinical.R` | GSE185263: sepsis vs healthy, joint in-hospital mortality and median SOFA models, recruitment-site and ED/ICU analyses. | `.../GSE185263/` |
| — | `08_DRT_clinical_template.R` | A configurable version of the clinical framework for new cohorts (see below). | user-defined |

## Scoring a new cohort

1. Add the cohort to `COHORTS` in `01_DRT_scores.R` (RNA-seq counts, metadata, optionally Olink), set `DERIVE <- FALSE`, and run it. The published locked weights and coefficients are built into the script, so no RESERVE-U data are needed. Set `reference` to the samples that define the z-score scale (for example, baseline sepsis samples in a longitudinal cohort).
2. Edit Section 1 of `08_DRT_clinical_template.R` (file paths, ID columns, covariates, case/control definition, outcomes, exposures, longitudinal settings), switch on the modules you need, and run it.

The template supports four modules: (A) case vs control Wilcoxon tests; (B) mutually adjusted D/R/T outcome models (logistic with Firth sensitivity, proportional-odds, median quantile, or linear regression, with VIFs); (C) exposure → score linear regressions with BH correction and exploratory exposure × domain interactions; (D) longitudinal mixed-effects models with a time × outcome-group interaction, time-specific contrasts, and a landmark change analysis.

## Notes

- Paths in the clinical scripts are set at the top of each script, relative to the project root.
- RESERVE-U-2-TOR projection scores used in all analyses are the out-of-fold predictions from the first cross-validation repeat (`OOF_METHOD = "first_repeat"`). The range of out-of-fold r across the 10 repeats is in `derivation/surrogate_crossvalidated_performance.csv`.
- Internal identifiers keep the original name "surrogate" (e.g. method `rna_surrogate_lasso14_1se`, file `surrogate_coefficients.csv`); in the manuscript these are the transcriptomic projections.
- Software versions are written to `sessionInfo.txt` in each output folder.
