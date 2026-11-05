# =============================================================================
# RESERVE-U-2-TOR (Tororo) participant characteristics table
# Script: 03a_RESERVE_U_2_TOR_Table_Characteristics.R
# =============================================================================
#
# INPUT
#   Clinical data, then TOR_DRT_scores_wide_perSD.csv (paths set below).
#
# OUTPUT
#   DRT_complete_clinical_Table_N_274.docx
# =============================================================================

# Clear R environment
rm(list = ls())

# Load libraries (retained from the original script)
library(boot)
library(broom)
library(dplyr)
library(FSA)
library(fmsb)
library(ggeffects)
library(ggbeeswarm)
library(ggplot2)
library(ggprism)
library(ggpubr)
library(lmtest)
library(MASS)
library(pROC)
library(purrr)
library(rstatix)
library(rms)
library(RVAideMemoire)
library(sandwich)
library(scales)
library(stringr)
library(tibble)
library(tidyr)
library(stats)
library(gmodels)
library(skimr)
library(gtsummary)
library(gt)

# Import the data: first the clinical data, then TOR_DRT_scores_wide_perSD.csv
# ------------------------------------------------------------
# INPUT / OUTPUT PATHS (relative to the project root; edit if needed).
# Scores are produced by publication_code/01_DRT_scores.R.
# ------------------------------------------------------------
clinical_file <- "data/tororo_all_298_clinical_data.csv"
score_file    <- "DRT_scores_output/TOR/TOR_DRT_scores_wide_perSD.csv"
clinical_data <- read.csv(clinical_file, header = TRUE, check.names = FALSE)
factor_scores <- read.csv(score_file, header = TRUE, check.names = FALSE)

# IMPORTANT: Restrict to clinical participants who appear in the multi-omic
# score file (sample_id in score file = pid in clinical file).
# The scores are used ONLY to select participants; no scores appear in Table 1.
pc1_columns <- c(
  "multiomic_PC1__Damage",
  "multiomic_PC1__Resistance",
  "multiomic_PC1__Tolerance"
)
stopifnot("pid" %in% names(clinical_data))
stopifnot(all(c("sample_id", pc1_columns) %in% names(factor_scores)))

scored_pids <- factor_scores %>%
  filter(if_all(all_of(pc1_columns), ~ !is.na(.))) %>%
  transmute(pid = trimws(as.character(sample_id))) %>%
  distinct()

clinical_data <- clinical_data %>%
  mutate(pid = trimws(as.character(pid)))

if (anyDuplicated(clinical_data$pid)) {
  stop("Clinical data have duplicate pid values; check before creating Table 1.")
}
if (anyNA(scored_pids$pid) || any(scored_pids$pid == "")) {
  stop("Missing or blank sample_id in the scored participant list.")
}
if (!all(scored_pids$pid %in% clinical_data$pid)) {
  stop("Some scored sample_id values do not match clinical pid values.")
}

reserve_tor <- clinical_data %>%
  semi_join(scored_pids, by = "pid")

cat("Clinical participants:", nrow(clinical_data), "\n")
cat("Multi-omic PC1 participants:", nrow(scored_pids), "\n")
cat("Clinical participants included in Table 1:", nrow(reserve_tor), "\n")
stopifnot(nrow(reserve_tor) == nrow(scored_pids))

# Convert 'pid' to fixed row identifier (original script behavior)
rownames(reserve_tor) <- reserve_tor$pid
reserve_tor$pid <- NULL

# Compute Shock Index and ROX Index
reserve_tor$shock_index <- reserve_tor$heart_rate_3 / reserve_tor$sbp_3

# Define FiO2 lookup based on flow rate
fio2_lookup <- c(
  `1` = 0.24, `2` = 0.28, `3` = 0.32, `4` = 0.36, `5` = 0.40,
  `6` = 0.44, `7` = 0.50, `8` = 0.60, `9` = 0.90, `10` = 0.95
)

# Calculate ROX index using o2_sat_3 and supp_o2_hosp
reserve_tor <- reserve_tor %>%
  mutate(
    fio2 = case_when(
      supp_o2_hosp == 0 ~ 0.21,
      supp_o2_hosp == 1 & !is.na(supp_02_hosp_amt) ~
        unname(fio2_lookup[as.character(round(supp_02_hosp_amt))]),
      supp_o2_hosp == 1 & is.na(supp_02_hosp_amt) ~ 0.40,
      TRUE ~ NA_real_
    ),
    rox_index = (o2_sat_3 / fio2) / resp_rate_3
  )

# Table of clinical/micro characteristics
# Recode discharge KPS only in patients discharged alive or transferred
reserve_tor[reserve_tor$hosp_outcome == 2, "dc_perf_16_and_above"] <- NA

# Original variable list: use actual clinical column 'sex' instead of 'gender';
# remove Factor1 (multi-omic scores are for eligibility only).
reserve_tor_table_df <- reserve_tor[, c(
  "age",
  "sex",
  "muac_imputed",
  "illnessduration_to_enroll",
  "admit_perf_16_and_above",
  "assist_oob_walk",
  "cough",
  "sob",
  "headache",
  "dysuria",
  "diarrhea",
  "vomiting",
  "rash",
  "tempmax_imputed",
  "heart_rate_3",
  "resp_rate_3",
  "sbp_3",
  "o2_sat_3",
  "ams_avpu",
  "lactate_result_imputed",
  "lactate4plus",
  "qsofa_1p_avpu",
  "qsofa_2p_avpu",
  "mews_score",
  "uva_score_avpu",
  "shock_index",
  "rox_index",
  "hiv_rdt_result",
  "hiv_stage_3_4_all",
  "art_prior_all",
  "cd4_imputed",
  "suppressedviralload_imputed",
  "malaria_rdt_result",
  "microtbdx",
  "urine_lam_result",
  "cryto_ag_result",
  "flu_pcr_uvri_result",
  "covid_pcr_uvri_result",
  "supp_o2_hosp",
  "ivf_hosp",
  "antimalarials_hosp",
  "abx_hosp",
  "anti_tb_hosp",
  "blood_transfusion_hosp",
  "hosp_outcome",
  "dc_perf_16_and_above",
  "patient_status_30d",
  "death30d",
  "kps_30d_16plus",
  "patient_status_60d",
  "death60d",
  "kps_60d_16plus"
)]

# Same gtsummary table as original, without Factor1
# Keep categorical denominators as in the original specification.
theme_gtsummary_compact()

table_1 <- reserve_tor_table_df %>%
  tbl_summary(
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n}/{N} ({p}%)"
    ),
    type = list(
      age ~ "continuous",
      sex ~ "categorical",
      muac_imputed ~ "continuous",
      illnessduration_to_enroll ~ "continuous",
      cough ~ "categorical",
      sob ~ "categorical",
      headache ~ "categorical",
      dysuria ~ "categorical",
      diarrhea ~ "categorical",
      vomiting ~ "categorical",
      rash ~ "categorical",
      admit_perf_16_and_above ~ "continuous",
      assist_oob_walk ~ "dichotomous",
      lactate_result_imputed ~ "continuous",
      lactate4plus ~ "dichotomous",
      qsofa_1p_avpu ~ "dichotomous",
      qsofa_2p_avpu ~ "dichotomous",
      mews_score ~ "continuous",
      uva_score_avpu ~ "continuous",
      shock_index ~ "continuous",
      rox_index ~ "continuous",
      hiv_rdt_result ~ "dichotomous",
      art_prior_all ~ "dichotomous",
      cd4_imputed ~ "continuous",
      suppressedviralload_imputed ~ "dichotomous",
      malaria_rdt_result ~ "dichotomous",
      microtbdx ~ "dichotomous",
      urine_lam_result ~ "dichotomous",
      cryto_ag_result ~ "dichotomous",
      flu_pcr_uvri_result ~ "dichotomous",
      covid_pcr_uvri_result ~ "dichotomous",
      supp_o2_hosp ~ "dichotomous",
      ivf_hosp ~ "dichotomous",
      antimalarials_hosp ~ "dichotomous",
      abx_hosp ~ "dichotomous",
      anti_tb_hosp ~ "dichotomous",
      blood_transfusion_hosp ~ "dichotomous",
      hosp_outcome ~ "categorical",
      dc_perf_16_and_above ~ "continuous",
      patient_status_30d ~ "categorical",
      death30d ~ "dichotomous",
      kps_30d_16plus ~ "continuous",
      patient_status_60d ~ "categorical",
      death60d ~ "dichotomous",
      kps_60d_16plus ~ "continuous"
    ),
    digits = list(
      all_continuous() ~ 0,
      all_dichotomous() ~ c(0, 0, 1),
      all_categorical() ~ c(0, 0, 1),
      tempmax_imputed ~ c(1, 1),
      lactate_result_imputed ~ c(1, 1),
      shock_index ~ c(1, 1),
      rox_index ~ c(2, 2)
    ),
    missing_text = "(Missing)",
    label = list(
      age ~ "Age",
      sex ~ "Sex",
      muac_imputed ~ "Mid-upper arm circumference, cm",
      illnessduration_to_enroll ~ "Illness duration prior to enrollment, days",
      admit_perf_16_and_above ~ "KPS at hospital admission",
      assist_oob_walk ~ "Unable to ambulate without assistance",
      cough ~ "Cough",
      sob ~ "Shortness of breath",
      headache ~ "Headache",
      dysuria ~ "Dysuria",
      diarrhea ~ "Diarrhea",
      vomiting ~ "Vomiting",
      rash ~ "Skin rash or discoloration",
      tempmax_imputed ~ "Temperature, \u00b0C",
      heart_rate_3 ~ "Heart rate, beats/min",
      resp_rate_3 ~ "Respiratory rate, breaths/min",
      sbp_3 ~ "Systolic blood pressure, mmHg",
      o2_sat_3 ~ "Oxygen saturation, %",
      ams_avpu ~ "Altered mental status",
      lactate_result_imputed ~ "Whole-blood lactate, mmol/L",
      lactate4plus ~ "Whole-blood lactate \u2265 4 mmol/L",
      qsofa_1p_avpu ~ "qSOFA \u2265 1",
      qsofa_2p_avpu ~ "qSOFA \u2265 2",
      mews_score ~ "Modified Early Warning Score",
      uva_score_avpu ~ "Universal Vital Assessment score",
      shock_index ~ "Shock index",
      rox_index ~ "ROX index",
      hiv_rdt_result ~ "PLWH",
      hiv_stage_3_4_all ~ "WHO HIV clinical stage 3 or 4",
      art_prior_all ~ "Receiving ART prior to hospitalization",
      cd4_imputed ~ "CD4 count, cells/mm3",
      suppressedviralload_imputed ~ "HIV-1 viral suppression",
      malaria_rdt_result ~ "Malaria RDT positive",
      microtbdx ~ "Microbiological TB positive",
      urine_lam_result ~ "Urine TB-LAM positive if PLWH",
      cryto_ag_result ~ "Cryptococcal ag. positive if PLWH",
      flu_pcr_uvri_result ~ "Influenza PCR positive",
      covid_pcr_uvri_result ~ "SARS-CoV-2 PCR positive",
      supp_o2_hosp ~ "Received oxygen therapy",
      ivf_hosp ~ "Received intravenous fluids",
      antimalarials_hosp ~ "Received anti-malarial agent(s)",
      abx_hosp ~ "Received anti-bacterial agent(s)",
      anti_tb_hosp ~ "Received anti-TB agent(s)",
      blood_transfusion_hosp ~ "Received red blood cell transfusion",
      hosp_outcome ~ "Hospital outcome",
      dc_perf_16_and_above ~ "KPS at discharge or transfer",
      patient_status_30d ~ "Vital status at 30 days post-discharge",
      death30d ~ "Death at 30 days post-discharge",
      kps_30d_16plus ~ "KPS at 30 days post-discharge",
      patient_status_60d ~ "Vital status at 60 days post-discharge",
      death60d ~ "Death at 60 days post-discharge",
      kps_60d_16plus ~ "KPS at 60 days post-discharge"
    )
  ) %>%
  add_n()

out_dir <- "DRT_scores_output/Clinical analyses/RESERVE_U_2_TOR"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(out_dir, "DRT_complete_clinical_Table_N_274.docx")

gt::gtsave(as_gt(table_1), filename = output_file)
