# =============================================================================
# RESERVE-U-1-EBB (Entebbe) clinical analyses of the D/R/T domain PC1 scores
# Script: 04_RESERVE_U_1_EBB_DRT_Clinical.R
# =============================================================================
#
# SCORES
#   Original multi-omic PC1: Resistance_PC1 / Tolerance_PC1 / Damage_PC1.
#   Merges EBB wide domain scores (sample_id or pid) to clinical data on pid.
#   Domain-score cohort filter: cohort == "EBB".
#   Uses the Entebbe clinical variable schema from the existing
#   RESERVE-U-1-EBB clinical-analysis workflow.
#
# ANALYSES
#   One script, all 3 domains together (faceted within each plot).
#   Outcome: 30-day mortality (death30d); transcriptomic projection methods
#   include all matched patients.
#   MEWS modeled with a mutually adjusted joint domain proportional-odds model,
#   matching the harmonized RESERVE-U domain-PC1 workflow.
#
# OUTPUT
#   CSV files and plot PDFs are exported to the output directory specified below.
#   Plots are displayed in R/RStudio and exported as landscape PDF files.
#
# DOMAIN COLOURS
#   Resistance = blue, Tolerance = green, Damage = orange
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggprism)
  library(ggpubr)
  library(ggbeeswarm)
  library(lmtest)
  library(sandwich)
  library(scales)
  library(broom)
})

# ------------------------------------------------------------
# 1. Import clinical data and domain PC1 scores
# ------------------------------------------------------------

# ------------------------------------------------------------
# INPUT / OUTPUT PATHS (relative to the project root; edit if needed).
# Scores are produced by publication_code/01_DRT_scores.R.
# ------------------------------------------------------------
clinical_file <- "data/entebbe_all_301_clinical_data.csv"
domain_file   <- "DRT_scores_output/EBB/EBB_DRT_scores_wide_perSD.csv"

clinical_data <- read.csv(
  clinical_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

domain_scores <- read.csv(
  domain_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Standardize gender -> sex if needed.
if ("gender" %in% names(clinical_data) && !"sex" %in% names(clinical_data)) {
  clinical_data <- clinical_data %>% dplyr::rename(sex = gender)
}

# Harmonize the two Entebbe score-file schemas before ANY downstream analysis.
# New score export: sample_id, multiomic_PC1__{Resistance,Tolerance,Damage}.
# Older score export: pid, multiomic_option2_PC1__{Resistance,Tolerance,Damage}.
# The locked external transcriptomic projection columns are unchanged. No scores are refitted.
if (!"pid" %in% names(clinical_data)) {
  stop("Entebbe clinical data must contain a 'pid' column.")
}
if (!"sample_id" %in% names(domain_scores) && !"pid" %in% names(domain_scores)) {
  stop("Entebbe scores must contain either 'sample_id' (new) or 'pid' (old).")
}

clinical_data$pid <- trimws(as.character(clinical_data$pid))
if ("sample_id" %in% names(domain_scores)) {
  score_pid <- trimws(as.character(domain_scores$sample_id))
  if ("pid" %in% names(domain_scores) &&
      !identical(score_pid, trimws(as.character(domain_scores$pid)))) {
    stop("The score file has conflicting 'sample_id' and 'pid' identifiers.")
  }
  domain_scores$pid <- score_pid
  cat("Using Entebbe score-file identifier: sample_id -> pid\n")
} else {
  domain_scores$pid <- trimws(as.character(domain_scores$pid))
  cat("Using Entebbe score-file identifier: pid\n")
}
if (anyNA(clinical_data$pid) || any(!nzchar(clinical_data$pid)) ||
    anyNA(domain_scores$pid) || any(!nzchar(domain_scores$pid))) {
  stop("Missing or blank participant IDs in the clinical or score file.")
}
if (anyDuplicated(clinical_data$pid)) stop("Duplicate pid values in clinical file.")
if (anyDuplicated(domain_scores$pid)) stop("Duplicate participant IDs in Entebbe score file.")

for (d in c("Resistance", "Tolerance", "Damage")) {
  new_col <- paste0("multiomic_PC1__", d)
  old_col <- paste0("multiomic_option2_PC1__", d)
  to <- paste0(d, "_PC1")
  if (!new_col %in% names(domain_scores) && !old_col %in% names(domain_scores)) {
    stop("Missing multi-omic PC1 score for ", d,
         ". Expected either ", new_col, " or ", old_col)
  }
  if (new_col %in% names(domain_scores) && old_col %in% names(domain_scores)) {
    a <- suppressWarnings(as.numeric(domain_scores[[new_col]]))
    b <- suppressWarnings(as.numeric(domain_scores[[old_col]]))
    if (!isTRUE(all.equal(a, b, tolerance = 1e-8, check.attributes = FALSE))) {
      stop("Conflicting old/new multi-omic PC1 values for ", d)
    }
  }
  from <- if (new_col %in% names(domain_scores)) new_col else old_col
  # Keep legacy aliases for the unchanged downstream correlations and plots.
  domain_scores[[old_col]] <- suppressWarnings(as.numeric(domain_scores[[from]]))
  domain_scores[[to]] <- domain_scores[[old_col]]
  cat("Mapped Entebbe ", from, " -> ", to, "\n", sep = "")
}

required_rna_cols <- as.vector(outer(
  c("rna_surrogate_lasso14_1se__"),
  c("Resistance", "Tolerance", "Damage"), paste0
))
required_score_cols <- c(
  "cohort", "pid", paste0(c("Resistance", "Tolerance", "Damage"), "_PC1"),
  paste0("multiomic_option2_PC1__", c("Resistance", "Tolerance", "Damage")),
  required_rna_cols
)
missing_score_cols <- setdiff(required_score_cols, names(domain_scores))
if (length(missing_score_cols) > 0) {
  stop("Missing columns in the EBB wide domain score file: ",
       paste(missing_score_cols, collapse = ", "))
}
# Only the declared multi-omic PC1 and transcriptomic projection measurements enter downstream data.
domain_scores <- domain_scores[, unique(required_score_cols), drop = FALSE]

domain_scores_ebb <- domain_scores %>%
  dplyr::mutate(cohort = trimws(as.character(cohort))) %>%
  dplyr::filter(cohort == "EBB") %>%
  dplyr::select(-cohort)
if (nrow(domain_scores_ebb) == 0) stop("No cohort == 'EBB' rows found in score file.")
if (anyDuplicated(domain_scores_ebb$pid)) stop("Duplicate pid values in EBB score file.")

# EBB has some participants with transcriptomic projection scores but no paired multi-omic PC1.
# Keep all matched patients for transcriptomic projection analyses; restrict only the original
# multi-omic clinical figures to participants with all 3 PC1 scores.
reserve_ebb_all <- clinical_data %>%
  dplyr::inner_join(domain_scores_ebb, by = "pid")

cat("\nClinical participants: ", nrow(clinical_data), "\n", sep = "")
cat("EBB score participants: ", nrow(domain_scores_ebb), "\n", sep = "")
cat("Clinical/score matched participants: ", nrow(reserve_ebb_all), "\n", sep = "")
cat("Participants with paired multi-omic PC1: ",
    sum(stats::complete.cases(reserve_ebb_all[, c("Resistance_PC1", "Tolerance_PC1", "Damage_PC1")])),
    "\n", sep = "")
if (nrow(reserve_ebb_all) == 0) stop("No Entebbe clinical participants matched scores by pid.")

# Require variables used by the replicated Tororo clinical workflow.
required_clinical_vars <- c(
  "age",
  "sex",
  "illnessduration_enroll",
  "MEWS_score",
  "death30d",
  "hivrdtresult",
  "malariardtresult",
  "microtbdx"
)

missing_clinical_vars <- setdiff(required_clinical_vars, names(reserve_ebb_all))
if (length(missing_clinical_vars) > 0) {
  stop(
    "Entebbe clinical data are missing variables required by this workflow: ",
    paste(missing_clinical_vars, collapse = ", ")
  )
}

# Coerce variables using the Entebbe schema.
# The supplied Entebbe template codes sex as 0 = Male, 1 = Female.
reserve_ebb_all <- reserve_ebb_all %>%
  dplyr::mutate(
    sex = dplyr::case_when(
      as.character(sex) %in% c("0", "Male", "male", "M") ~ "Male",
      as.character(sex) %in% c("1", "Female", "female", "F") ~ "Female",
      TRUE ~ as.character(sex)
    ),
    sex = factor(sex, levels = c("Male", "Female")),
    age = suppressWarnings(as.numeric(age)),
    illnessduration_enroll = suppressWarnings(as.numeric(illnessduration_enroll)),
    MEWS_score = suppressWarnings(as.numeric(as.character(MEWS_score))),
    death30d = suppressWarnings(as.numeric(as.character(death30d))),
    hivrdtresult = factor(as.character(hivrdtresult), levels = c("0", "1")),
    malariardtresult = factor(as.character(malariardtresult), levels = c("0", "1")),
    microtbdx = factor(as.character(microtbdx), levels = c("0", "1"))
  )

# The multi-omic figures/models use participants with complete PC1 scores.
# Correlation and transcriptomic projection comparison sections instead use reserve_ebb_all.
reserve_ebb <- reserve_ebb_all %>%
  dplyr::filter(!is.na(Resistance_PC1), !is.na(Tolerance_PC1), !is.na(Damage_PC1))
if (nrow(reserve_ebb) == 0) stop("No participants with complete multi-omic PC1 scores.")

# Output directory (forward slashes work correctly on Windows in R).
out_dir <- "DRT_scores_output/Clinical analyses/RESERVE_U_1_EBB"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(out_dir)) {
  stop("Output directory could not be created or accessed: ", out_dir)
}


# Save plots to the same output directory used for CSV files.
# Each plot is still printed to the R/RStudio graphics device.

save_plot <- function(plot, filename) {
  # Journal figures: vector PDF only, exactly 8.59 × 6.54 in (landscape).
  plot <- plot + ggplot2::labs(caption = NULL) +
    ggplot2::theme(plot.caption = ggplot2::element_blank())
  ggplot2::ggsave(
    filename = file.path(out_dir, filename), plot = plot,
    device = "pdf", width = 8.59, height = 6.54, units = "in", bg = "white")
}


# Prevent intermittent RStudio interactive grid errors from interrupting exports.
safe_print_plot <- function(plot, label = "plot") {
  tryCatch({
    print(plot)
    invisible(TRUE)
  }, error = function(e) {
    warning("Interactive display failed for ", label, ": ",
            conditionMessage(e),
            ". Export will still be attempted; consider graphics.off().",
            call. = FALSE)
    invisible(FALSE)
  })
}

write.csv(
  reserve_ebb[, !grepl("mean.?z", names(reserve_ebb), ignore.case = TRUE), drop = FALSE],
  file.path(out_dir, "Entebbe_merged_clinical_domain_scores.csv"),
  row.names = FALSE
)

write.csv(
  reserve_ebb_all[, !grepl("mean.?z", names(reserve_ebb_all), ignore.case = TRUE), drop = FALSE],
  file.path(out_dir, "RESERVE_U_1_EBB_all_matched_DRT_clinical_scores.csv"),
  row.names = FALSE
)

# Simple merge-QC table.
merge_qc <- tibble::tibble(
  metric = c(
    "Clinical rows",
    "EBB domain-score rows",
    "Clinical/score matched rows (all transcriptomic projection methods)",
    "Merged rows with all three PC1 scores",
    "Clinical pids without EBB PC1 scores",
    "EBB PC1-score pids without clinical row"
  ),
  n = c(
    nrow(clinical_data),
    nrow(domain_scores_ebb),
    nrow(reserve_ebb_all),
    nrow(reserve_ebb),
    sum(!clinical_data$pid %in% domain_scores_ebb$pid),
    sum(!domain_scores_ebb$pid %in% clinical_data$pid)
  )
)

print(merge_qc)

write.csv(
  merge_qc,
  file.path(out_dir, "Entebbe_domain_PC1_merge_QC.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Long format for all 3 domains
# ------------------------------------------------------------

domain_order <- c("Resistance_PC1", "Tolerance_PC1", "Damage_PC1")
domain_labels <- c(
  Resistance_PC1 = "Resistance",
  Tolerance_PC1  = "Tolerance",
  Damage_PC1     = "Damage"
)

domain_colors <- c(
  Resistance = "#2A6EBB",
  Tolerance  = "#2CA25F",
  Damage     = "#E6550D"
)

reserve_long <- reserve_ebb %>%
  pivot_longer(
    cols = all_of(domain_order),
    names_to = "Domain",
    values_to = "Domain_PC1"
  ) %>%
  mutate(
    Domain = factor(Domain, levels = domain_order, labels = unname(domain_labels))
  )

# ============================================================
# ADDITIONAL ANALYSES: multi-omic / transcriptomic projection correlations and associations
# Transcriptomic projection scores are not substituted for the primary multi-omic PC1 scores.
# The primary clinical outcome models and figures use multi-omic PC1.
# ============================================================

# ------------------------------------------------------------
# A. Correlations between score construction approaches
# ------------------------------------------------------------
# For RESERVE-U-1-EBB, LASSO correlations use locked predictions from
# the model developed in RESERVE-U-2-TOR (no refitting in EBB).
correlation_methods <- c(
  "Multi-omic PC1" = "multiomic_option2_PC1__",
  "Transcriptomic projection (locked external)" = "rna_surrogate_lasso14_1se__"
)

comparison_pairs <- tibble::tribble(
  ~method_1, ~method_2, ~comparison,
  "Multi-omic PC1", "Transcriptomic projection (locked external)", "Multi-omic PC1 vs. transcriptomic projection"
)

correlation_results <- purrr::map_dfr(seq_len(nrow(comparison_pairs)), function(i) {
  purrr::map_dfr(c("Damage", "Resistance", "Tolerance"), function(d) {
    x_col <- paste0(correlation_methods[[comparison_pairs$method_1[i]]], d)
    y_col <- paste0(correlation_methods[[comparison_pairs$method_2[i]]], d)
    paired <- stats::complete.cases(reserve_ebb_all[, c(x_col, y_col), drop = FALSE])
    tibble::tibble(
      comparison = comparison_pairs$comparison[i],
      domain = d,
      method_1 = comparison_pairs$method_1[i],
      method_2 = comparison_pairs$method_2[i],
      n = sum(paired),
      pearson_r = if (sum(paired) > 2) {
        stats::cor(reserve_ebb_all[[x_col]][paired], reserve_ebb_all[[y_col]][paired])
      } else NA_real_,
      spearman_r = if (sum(paired) > 2) {
        stats::cor(reserve_ebb_all[[x_col]][paired], reserve_ebb_all[[y_col]][paired], method = "spearman")
      } else NA_real_
    )
  })
})

utils::write.csv(
  correlation_results,
  file.path(out_dir, "RESERVE_U_1_EBB_DRT_scoring_correlations_long.csv"),
  row.names = FALSE
)

correlation_table <- correlation_results %>%
  dplyr::select(comparison, domain, pearson_r) %>%
  tidyr::pivot_wider(names_from = domain, values_from = pearson_r) %>%
  dplyr::mutate(dplyr::across(c(Damage, Resistance, Tolerance), ~ round(.x, 3)))

utils::write.csv(
  correlation_table,
  file.path(out_dir, "RESERVE_U_1_EBB_DRT_scoring_correlation_matrix.csv"),
  row.names = FALSE
)
cat("\nPairwise Pearson correlations (transcriptomic projection = locked external):\n")
print(as.data.frame(correlation_table), row.names = FALSE)

# A heatmap formatted like the supplied domain comparison table:
# comparisons as rows, domains as columns, Pearson r printed in tiles.
# Focus the sequential blue color scale on the observed range (about 0.55-0.96),
# instead of 0-1. The higher-contrast palette makes
# modest versus very strong correlations visually distinguishable.
# The numeric correlation labels remain the exact computed values.
plot_correlation_data <- correlation_results %>%
  dplyr::mutate(
    comparison = factor(comparison, levels = rev(comparison_pairs$comparison)),
    domain = factor(domain, levels = c("Damage", "Resistance", "Tolerance")),
    r_label = ifelse(is.na(pearson_r), "NA", sprintf("%.3f", pearson_r)),
    label_color = dplyr::if_else(!is.na(pearson_r) & pearson_r >= 0.82,
                                 "white", "black")
  )

p_drt_correlation_matrix <- ggplot2::ggplot(
  plot_correlation_data, ggplot2::aes(x = domain, y = comparison, fill = pearson_r)
) +
  ggplot2::geom_tile(color = "white", linewidth = 1.3) +
  ggplot2::geom_text(
    ggplot2::aes(label = r_label, color = label_color),
    size = 5, fontface = "bold"
  ) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_fill_gradientn(
    colours = c("#DEEBF7", "#9ECAE1", "#4292C6", "#2171B5", "#084594"),
    values = scales::rescale(c(0.50, 0.62, 0.74, 0.86, 1.00),
                             from = c(0.50, 1.00)),
    limits = c(0.50, 1.00),
    breaks = c(0.50, 0.60, 0.70, 0.80, 0.90, 1.00),
    labels = c("0.50", "0.60", "0.70", "0.80", "0.90", "1.00"),
    oob = scales::squish,
    na.value = "gray90",
    name = "Pearson r",
    guide = ggplot2::guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(9, "cm"),
      barheight = grid::unit(0.4, "cm")
    )
  ) +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::labs(
    title = "Correlations between domain scoring approaches",
    subtitle = "RESERVE-U-1-EBB | Pairwise-complete Pearson correlations",
    x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(face = "bold", color = "black", size = 12),
    axis.text.y = ggplot2::element_text(color = "black", size = 11),
    plot.title = ggplot2::element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 11, face = "bold"),
    legend.text = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::margin(10, 16, 12, 12)
  )
safe_print_plot(p_drt_correlation_matrix, "p_drt_correlation_matrix")


ggplot2::ggsave(
  filename = file.path(out_dir, "RESERVE_U_1_EBB_DRT_scoring_correlation_matrix.pdf"),
  plot = p_drt_correlation_matrix, bg = "white",
    width = 8.59, height = 6.54, units = "in", device = "pdf")

# ------------------------------------------------------------
# B. Transcriptomic projection: JOINT domain sensitivity models (no separate-domain models)
# ------------------------------------------------------------
# Entebbe: locked external predictions from the Tororo-trained transcriptomic projection.
# Compare with primary multi-omic PC1 on the PC1-available analysis cohort.
# Standardize each transcriptomic projection domain once on the transcriptome-complete matched cohort BEFORE
# restricting to outcome-specific complete cases. Odds ratios are per 1 SD.
# All three transcriptomic projection domains enter every severity/mortality model TOGETHER.

rna_domains <- c("Resistance", "Tolerance", "Damage")
rna_prefix <- "rna_surrogate_lasso14_1se__"
rna_raw_columns <- paste0(rna_prefix, rna_domains)
rna_missing <- setdiff(rna_raw_columns, names(reserve_ebb))
if (length(rna_missing)) stop("Missing transcriptomic projection score columns: ",
                             paste(rna_missing, collapse = ", "))

rna_joint_data <- reserve_ebb
for (d in rna_domains) {
  raw <- suppressWarnings(as.numeric(rna_joint_data[[paste0(rna_prefix, d)]]))
  if (all(!is.finite(raw)) || stats::sd(raw, na.rm = TRUE) <= 0)
    stop("Transcriptomic projection domain missing or invariant: ", d)
  rna_joint_data[[paste0(d, "_rna_z")]] <-
    (raw - mean(raw, na.rm = TRUE)) / stats::sd(raw, na.rm = TRUE)
}

rna_terms <- paste0(rna_domains, "_rna_z")
rna_covariates <- c("age", "sex", "illnessduration_enroll", "hivrdtresult", "malariardtresult", "microtbdx")
rna_predictors <- c(rna_terms, rna_covariates)

extract_rna_joint <- function(fit, outcome_name, analysis_data) {
  b <- stats::coef(fit)[rna_terms]
  se <- sqrt(diag(stats::vcov(fit)))[rna_terms]
  if (any(!is.finite(b)) || any(!is.finite(se)) || any(se <= 0))
    stop("Transcriptomic projection joint model has nonestimable domain coefficients: ", outcome_name)
  pp <- 2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
  tibble::tibble(
    method = "Transcriptomic projection (locked external)",
    outcome = outcome_name,
    domain = rna_domains,
    n = stats::nobs(fit),
    deaths = if (outcome_name == "30-day mortality") sum(analysis_data$death30d == 1L) else NA_integer_,
    log_OR = unname(b), SE = unname(se),
    OR = exp(unname(b)), lower95 = exp(unname(b - 1.96 * se)),
    upper95 = exp(unname(b + 1.96 * se)),
    p = unname(pp), q_BH_3domains = stats::p.adjust(pp, method = "BH")
  )
}

rna_mortality_data <- rna_joint_data[
  stats::complete.cases(rna_joint_data[, c("death30d", rna_predictors), drop = FALSE]),
  , drop = FALSE]
rna_mortality_data <- droplevels(rna_mortality_data)
if (!all(rna_mortality_data$death30d %in% c(0, 1)) ||
    length(unique(rna_mortality_data$death30d)) != 2L)
  stop("Transcriptomic projection mortality requires both 30-day outcome classes")
rna_mortality_formula <- stats::reformulate(rna_predictors, response = "death30d")
rna_mortality_fit <- stats::glm(rna_mortality_formula, data = rna_mortality_data,
                                family = stats::binomial())

rna_mews_data <- rna_joint_data[
  stats::complete.cases(rna_joint_data[, c("MEWS_score", rna_predictors), drop = FALSE]),
  , drop = FALSE]
rna_mews_data <- droplevels(rna_mews_data)
rna_mews_data$RNA_MEWS_ordinal <- ordered(
  rna_mews_data$MEWS_score, levels = sort(unique(rna_mews_data$MEWS_score)))
if (nlevels(rna_mews_data$RNA_MEWS_ordinal) < 3L)
  stop("Transcriptomic projection joint MEWS requires >=3 ordinal categories")
rna_mews_formula <- stats::reformulate(rna_predictors, response = "RNA_MEWS_ordinal")
rna_mews_fit <- MASS::polr(rna_mews_formula, data = rna_mews_data,
                            Hess = TRUE, method = "logistic")

rna_joint_estimates <- dplyr::bind_rows(
  extract_rna_joint(rna_mortality_fit, "30-day mortality", rna_mortality_data),
  extract_rna_joint(rna_mews_fit, "Ordinal MEWS", rna_mews_data)
)
utils::write.csv(rna_joint_estimates,
  file.path(out_dir, "RESERVE_U_1_EBB_RNA_LASSO_JOINT_DRT_30d_mortality_MEWS.csv"),
  row.names = FALSE)

rna_vif <- function(dd, label) {
  purrr::map_dfr(rna_terms, function(term) {
    aux <- stats::lm(stats::reformulate(setdiff(rna_predictors, term),
                                        response = term), data = dd)
    r2 <- summary(aux)$r.squared
    tibble::tibble(outcome = label, domain = sub("_rna_z$", "", term),
                   n = stats::nobs(aux), VIF = 1 / (1 - r2))
  })
}
rna_joint_vif <- dplyr::bind_rows(
  rna_vif(rna_mortality_data, "30-day mortality"),
  rna_vif(rna_mews_data, "Ordinal MEWS")
)
utils::write.csv(rna_joint_vif,
  file.path(out_dir, "RESERVE_U_1_EBB_RNA_LASSO_JOINT_DRT_VIF.csv"),
  row.names = FALSE)
cat("\nTranscriptomic projection joint mortality / MEWS (per 1 SD):\n")
print(rna_joint_estimates)


# ============================================================
# PRIMARY ANALYSES: mutually adjusted multi-omic domain scores
# - Joint 30-day mortality logistic regression
# - Joint MEWS proportional-odds regression
# - Firth penalized 30-day mortality sensitivity analysis
# - Same six clinical covariates in both RESERVE-U cohorts
# - Domain scores standardized within cohort (1-SD scale)
# ============================================================

primary_domain_names <- c("Resistance", "Tolerance", "Damage")
primary_domain_cols <- paste0(primary_domain_names, "_PC1")
primary_covariates <- c("age", "sex", "illnessduration_enroll", "hivrdtresult", "malariardtresult", "microtbdx")
primary_domain_terms <- paste0(primary_domain_names, "_z")
primary_predictors <- c(primary_domain_terms, primary_covariates)
primary_mortality_formula <- stats::reformulate(primary_predictors, response = "death30d")
primary_mews_formula <- stats::reformulate(primary_predictors, response = "MEWS_ordinal")

# Standardization is based on the original multi-omic analysis population;
# the identical domain scaling is used for the mortality and MEWS models.
standardize_domain <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  sx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx <= 0) stop("A multi-omic domain has zero/undefined SD")
  (x - mean(x, na.rm = TRUE)) / sx
}

primary_data <- reserve_ebb %>%
  dplyr::mutate(
    Resistance_z = standardize_domain(Resistance_PC1),
    Tolerance_z = standardize_domain(Tolerance_PC1),
    Damage_z = standardize_domain(Damage_PC1),
    death30d = suppressWarnings(as.numeric(as.character(death30d))),
    MEWS_raw = suppressWarnings(as.numeric(as.character(MEWS_score)))
  )

mortality_joint_data <- primary_data %>%
  dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(c("death30d", primary_predictors))))) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(c("hivrdtresult", "malariardtresult", "microtbdx")), droplevels))
if (!all(unique(mortality_joint_data$death30d) %in% c(0, 1)) ||
    dplyr::n_distinct(mortality_joint_data$death30d) != 2L) {
  stop("The 30-day mortality variable must be binary, with both outcome groups present.")
}

mews_joint_data <- primary_data %>%
  dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(c("MEWS_raw", primary_predictors))))) %>%
  dplyr::mutate(
    MEWS_ordinal = ordered(MEWS_raw, levels = sort(unique(MEWS_raw))),
    dplyr::across(dplyr::all_of(c("hivrdtresult", "malariardtresult", "microtbdx")), droplevels)
  )
if (nlevels(mews_joint_data$MEWS_ordinal) < 3L) {
  stop("Insufficient distinct MEWS levels for a joint proportional-odds model")
}

primary_logistic_fit <- stats::glm(primary_mortality_formula,
                                   data = mortality_joint_data,
                                   family = stats::binomial())
primary_mews_fit <- MASS::polr(primary_mews_formula,
                               data = mews_joint_data, Hess = TRUE,
                               method = "logistic")

extract_joint_or <- function(fit, model_name, dat) {
  cf <- stats::coef(fit)
  vv <- stats::vcov(fit)
  ss <- sqrt(diag(vv))
  tibble::tibble(
    model = model_name,
    domain = primary_domain_names,
    term = primary_domain_terms,
    n = stats::nobs(fit),
    deaths = if ("death30d" %in% names(dat)) sum(dat$death30d == 1) else NA_integer_,
    log_OR = unname(cf[primary_domain_terms]),
    SE = unname(ss[primary_domain_terms])
  ) %>%
    dplyr::mutate(
      OR = exp(log_OR),
      lower95 = exp(log_OR - 1.96 * SE),
      upper95 = exp(log_OR + 1.96 * SE),
      p = 2 * stats::pnorm(abs(log_OR / SE), lower.tail = FALSE),
      q_BH = stats::p.adjust(p, method = "BH"),
      estimate_95CI = sprintf("%.2f (%.2f–%.2f)", OR, lower95, upper95)
    )
}

# MEWS data contain no mortality endpoint for reporting.
primary_mortality_table <- extract_joint_or(
  primary_logistic_fit, "Primary joint 30-day mortality", mortality_joint_data)
primary_mews_table <- extract_joint_or(
  primary_mews_fit, "Primary joint ordinal MEWS", mews_joint_data %>% dplyr::select(-dplyr::any_of("death30d")))
primary_joint_results <- dplyr::bind_rows(primary_mortality_table, primary_mews_table)
utils::write.csv(primary_joint_results,
  file.path(out_dir, "RESERVE_U_1_EBB_PRIMARY_joint_DRT_30d_Mortality_and_MEWS.csv"),
  row.names = FALSE)

# Report VIF for each molecular domain, conditional on other domains and
# the six clinical covariates. Distinct diagnostic samples for each outcome.
primary_domain_vif <- function(dat, model_name) {
  purrr::map_dfr(primary_domain_terms, function(term) {
    auxiliary <- stats::reformulate(setdiff(primary_predictors, term), response = term)
    fit <- stats::lm(auxiliary, data = dat)
    tibble::tibble(model = model_name, domain = sub("_z$", "", term),
                   n = stats::nobs(fit), VIF = 1/(1-summary(fit)$r.squared))
  })
}
primary_vifs <- dplyr::bind_rows(
  primary_domain_vif(mortality_joint_data, "Joint 30-day mortality"),
  primary_domain_vif(mews_joint_data, "Joint ordinal MEWS"))
utils::write.csv(primary_vifs,
  file.path(out_dir, "RESERVE_U_1_EBB_PRIMARY_joint_DRT_VIF.csv"), row.names = FALSE)

# Firth is a sensitivity, NOT a replacement for the primary logistic model.
# Install once in R if necessary: install.packages("logistf")
if (requireNamespace("logistf", quietly = TRUE)) {
  firth_fit <- tryCatch(
    logistf::logistf(primary_mortality_formula, data = mortality_joint_data,
                     firth = TRUE, pl = TRUE),
    error = function(e) { warning("Firth logistic regression failed: ", conditionMessage(e)); NULL }
  )
} else {
  firth_fit <- NULL
  warning("Package 'logistf' not installed. Skipping Firth sensitivity; install.packages('logistf') to enable.")
}

if (!is.null(firth_fit)) {
  firth_table <- tibble::tibble(
    domain = primary_domain_names, term = primary_domain_terms,
    n = nrow(mortality_joint_data), deaths = sum(mortality_joint_data$death30d == 1),
    OR = exp(unname(firth_fit$coefficients[primary_domain_terms])),
    lower95 = exp(unname(firth_fit$ci.lower[primary_domain_terms])),
    upper95 = exp(unname(firth_fit$ci.upper[primary_domain_terms])),
    p = unname(firth_fit$prob[primary_domain_terms]),
    note = "Firth logistic regression; profile penalized-likelihood CI"
  ) %>% dplyr::mutate(q_BH = stats::p.adjust(p, method = "BH"))
} else {
  firth_table <- tibble::tibble(
    domain = primary_domain_names, term = primary_domain_terms,
    n = nrow(mortality_joint_data), deaths = sum(mortality_joint_data$death30d == 1),
    OR = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
    p = NA_real_, note = "Firth fit unavailable; see console warning", q_BH = NA_real_
  )
}
utils::write.csv(firth_table,
  file.path(out_dir, "RESERVE_U_1_EBB_Firth_joint_DRT_30d_mortality_sensitivity.csv"),
  row.names = FALSE)

# Forest plots display mutually adjusted ORs, unlike the older single-domain
# figures retained below for descriptive comparison.
make_primary_joint_forest <- function(tbl, plot_title, x_axis_label) {
  plot_df <- tbl %>% dplyr::mutate(
    domain = factor(domain, levels = rev(primary_domain_names)),
    annotation = sprintf("%.2f (%.2f–%.2f); p=%.3g", OR, lower95, upper95, p)
  )
  ggplot2::ggplot(plot_df, ggplot2::aes(x = OR, y = domain, color = domain)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower95, xmax = upper95),
                          orientation = "y", width = 0.14, linewidth = 0.85) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::geom_text(ggplot2::aes(x = upper95 * 1.07, label = annotation),
                       hjust = 0, size = 3.5, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = domain_colors) +
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.08, 0.55))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = x_axis_label, y = NULL,
                  title = plot_title) +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(legend.position = "none",
                   plot.margin = ggplot2::margin(8, 160, 8, 8),
                   plot.title = ggplot2::element_text(face = "bold"),
                   axis.text.y = ggplot2::element_text(face = "bold"))
}

p_primary_joint_mortality <- make_primary_joint_forest(
  primary_mortality_table, "RESERVE-U-1-EBB: primary joint domain and 30-day mortality",
  "Adjusted odds ratio for 30-day mortality per 1-SD increase")
p_primary_joint_mews <- make_primary_joint_forest(
  primary_mews_table, "RESERVE-U-1-EBB: primary joint domain and MEWS",
  "Adjusted odds ratio for MEWS per 1-SD increase")
safe_print_plot(p_primary_joint_mortality, "p_primary_joint_mortality")
save_plot(p_primary_joint_mortality, "RESERVE_U_1_EBB_PRIMARY_joint_DRT_30d_mortality_forest.pdf")
safe_print_plot(p_primary_joint_mews, "p_primary_joint_mews")
save_plot(p_primary_joint_mews, "RESERVE_U_1_EBB_PRIMARY_joint_DRT_MEWS_forest.pdf")

# Conditional predicted mortality from the SAME joint logistic model;
# varying one domain at a time with the other two at their cohort means (0 SD).
# Clinical covariates held at cohort median/reference values.
make_primary_mortality_curves <- function(fit, dd) {
  template <- dd[rep(1L, 120L), , drop = FALSE]
  for (v in primary_domain_terms) template[[v]] <- 0
  for (v in primary_covariates) {
    if (is.factor(dd[[v]])) {
      lev <- if ("0" %in% levels(dd[[v]])) "0" else levels(dd[[v]])[1]
      template[[v]] <- factor(lev, levels = levels(dd[[v]]))
    } else {
      template[[v]] <- stats::median(dd[[v]], na.rm = TRUE)
    }
  }
  curves <- purrr::map_dfr(primary_domain_names, function(domain) {
    score_name <- paste0(domain, "_z")
    nd <- template
    nd[[score_name]] <- seq(stats::quantile(dd[[score_name]], .05),
                            stats::quantile(dd[[score_name]], .95), length.out = nrow(nd))
    pred <- stats::predict(fit, newdata = nd, type = "link", se.fit = TRUE)
    tibble::tibble(domain = domain, x = nd[[score_name]],
                   fit = stats::plogis(pred$fit),
                   lower = stats::plogis(pred$fit - 1.96 * pred$se.fit),
                   upper = stats::plogis(pred$fit + 1.96 * pred$se.fit))
  })
  # The annotation estimates come directly from this same joint model.
  # Keep Wald 95% CIs/p-values consistent with the primary joint forest/CSV.
  terms <- paste0(primary_domain_names, "_z")
  coef_joint <- stats::coef(fit)
  se_joint <- sqrt(diag(stats::vcov(fit)))
  labels <- tibble::tibble(
    domain = primary_domain_names,
    x = -Inf, y = Inf,
    annotation = sprintf(
      "OR = %.2f (%.2f–%.2f)\np = %.3g",
      exp(coef_joint[terms]),
      exp(coef_joint[terms] - 1.96 * se_joint[terms]),
      exp(coef_joint[terms] + 1.96 * se_joint[terms]),
      2 * stats::pnorm(abs(coef_joint[terms] / se_joint[terms]), lower.tail = FALSE)
    )
  )
  ggplot2::ggplot(curves, ggplot2::aes(x=x, y=fit, color=domain, fill=domain)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper), alpha=0.16, color=NA) +
    ggplot2::geom_line(linewidth=1.1) +
    ggplot2::geom_label(
      data=labels,
      ggplot2::aes(x=x, y=y, label=annotation, color=domain),
      inherit.aes=FALSE, hjust=-0.03, vjust=1.10,
      size=3.5, fontface="bold", fill="white",
      label.size=0, label.padding=grid::unit(0.14, "lines"),
      show.legend=FALSE
    ) +
    ggplot2::facet_wrap(~domain, nrow=1, scales="free_x") +
    ggplot2::scale_color_manual(values=domain_colors) +
    ggplot2::scale_fill_manual(values=domain_colors) +
    ggplot2::labs(x="Domain score (within-cohort SD)", y="Adjusted 30-day mortality probability",
                  title="RESERVE-U-1-EBB: 30-day mortality, joint domain model") +
    ggprism::theme_prism(base_size=13) +
    ggplot2::theme(legend.position="none", strip.text=ggplot2::element_text(face="bold"))
}
p_primary_joint_mortality_curve <- make_primary_mortality_curves(
  primary_logistic_fit, mortality_joint_data)
safe_print_plot(p_primary_joint_mortality_curve, "p_primary_joint_mortality_curve")
save_plot(p_primary_joint_mortality_curve,
  "RESERVE_U_1_EBB_PRIMARY_joint_DRT_30d_mortality_probabilities.pdf")

cat("\nPRIMARY JOINT domain MODELS:\n")
print(primary_joint_results)
cat("\nFIRTH MORTALITY SENSITIVITY:\n")
print(firth_table)
cat("\nVIF diagnostics:\n")
print(primary_vifs)

# ------------------------------------------------------------
# 4. Model map / covariate transparency
# ------------------------------------------------------------

model_map <- tibble::tribble(
  ~analysis_block, ~model_type, ~covariates_used,
  "PRIMARY joint domain 30-day mortality", "Binary logistic regression", "All 3 PC1 scores (within-cohort SD), age, sex, illness duration, HIV, malaria, TB",
  "PRIMARY joint domain MEWS", "Proportional-odds ordinal regression", "All 3 PC1 scores (within-cohort SD), age, sex, illness duration, HIV, malaria, TB",
  "Firth joint 30-day mortality", "Firth penalized logistic regression", "Same predictors/covariates as primary 30-day model",
  "Continuous PC1 plots", "Robust linear regression", "none",
  "Sex PC1 plot", "Wilcoxon test", "none",
  "Transcriptomic projection joint 30-day mortality and MEWS (supportive)", "Joint binary/ordinal logistic", "All 3 transcriptomic projection scores (per SD), age, sex, illness duration, HIV, malaria, TB",
  "Domain-specific pathogen-PC1 forest plots", "Linear regression", "age, sex, illnessduration_enroll, other major pathogen groups"
)

write.csv(
  model_map,
  file.path(out_dir, "Entebbe_domain_PC1_model_covariate_map.csv"),
  row.names = FALSE
)

print(model_map)

# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------

safe_domain_text_label <- function(est, lo, hi, p) {
  paste0(
    "OR = ", signif(est, 3),
    " (", signif(lo, 3), " - ", signif(hi, 3), ")",
    "\np = ", signif(p, 3)
  )
}

make_continuous_plot <- function(data, xvar, xlab, ylab = "Domain PC1 score",
                                 title = NULL, x_limits = NULL,
                                 adjust_covariates = NULL,
                                 caption_text = "Unadjusted robust linear regression") {
  dat <- data %>% filter(!is.na(.data[[xvar]]), !is.na(Domain_PC1))
  dat <- dat %>% mutate(Domain = factor(Domain, levels = levels(reserve_long$Domain)))
  
  if (!is.null(adjust_covariates) && length(adjust_covariates) > 0) {
    model_formula <- reformulate(c(xvar, adjust_covariates), response = "Domain_PC1")
    beta_prefix <- "Adjusted β"
  } else {
    model_formula <- reformulate(xvar, response = "Domain_PC1")
    beta_prefix <- "β"
  }
  
  lab_df <- dat %>%
    group_by(Domain) %>%
    group_modify(~ {
      dd <- .x
      if (nrow(dd) < 5) {
        return(tibble(label = "n<5", x = Inf, y = Inf))
      }
      fit <- tryCatch(MASS::rlm(model_formula, data = dd), error = function(e) NULL)
      if (is.null(fit)) {
        return(tibble(label = "model failed", x = Inf, y = Inf))
      }
      ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC0"))
      if (!(xvar %in% rownames(ct))) {
        return(tibble(label = "model failed", x = Inf, y = Inf))
      }
      tibble(
        label = paste0(beta_prefix, " = ", signif(ct[xvar, 1], 3), ", p = ", signif(ct[xvar, 4], 3)),
        x = Inf,
        y = Inf
      )
    }) %>%
    ungroup()
  

  # When adjustment covariates are specified, the fitted line corresponds
  # to the adjusted robust model instead of the unadjusted smoother.
  adjusted_line_df <- NULL
  if (!is.null(adjust_covariates) && length(adjust_covariates) > 0) {
    adjusted_line_df <- dat %>%
      dplyr::group_by(Domain) %>%
      dplyr::group_modify(~ {
        dd <- .x
        cc <- stats::complete.cases(dd[, unique(c(xvar, adjust_covariates)), drop = FALSE])
        dd <- dd[cc, , drop = FALSE]
        if (nrow(dd) < 10L || dplyr::n_distinct(dd[[xvar]]) < 3L)
          return(tibble::tibble(x_value = numeric(), model_y = numeric()))
        fit <- tryCatch(MASS::rlm(model_formula, data = dd), error = function(e) NULL)
        if (is.null(fit))
          return(tibble::tibble(x_value = numeric(), model_y = numeric()))
        nd <- dd[rep(1L, 100L), , drop = FALSE]
        for (v in adjust_covariates) {
          if (is.factor(dd[[v]])) {
            lv <- if ("0" %in% levels(dd[[v]])) "0" else levels(dd[[v]])[1]
            nd[[v]] <- factor(lv, levels = levels(dd[[v]]))
          } else {
            nd[[v]] <- stats::median(dd[[v]], na.rm = TRUE)
          }
        }
        nd[[xvar]] <- seq(min(dd[[xvar]]), max(dd[[xvar]]), length.out = nrow(nd))
        tibble::tibble(x_value = nd[[xvar]],
                       model_y = as.numeric(stats::predict(fit, newdata = nd)))
      }) %>% dplyr::ungroup()
  }

  p <- ggplot(dat, aes(x = .data[[xvar]], y = Domain_PC1, color = Domain, fill = Domain)) +
    geom_point(alpha = 0.8) +
    {
      if (is.null(adjusted_line_df)) {
        geom_smooth(method = MASS::rlm, formula = y ~ x, se = TRUE, alpha = 0.20)
      } else {
        geom_line(data = adjusted_line_df,
                  aes(x = x_value, y = model_y, color = Domain, group = Domain),
                  inherit.aes = FALSE, linewidth = 1.2)
      }
    } +
    geom_text(
      data = lab_df,
      aes(x = x, y = y, label = label, color = Domain),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.5,
      size = 4.5,
      fontface = "bold"
    ) +
    facet_wrap(~ Domain, scales = "free_y", nrow = 1) +
    scale_color_manual(values = domain_colors) +
    scale_fill_manual(values = domain_colors) +
    labs(
      x = xlab,
      y = ylab,
      title = title) +
    ggprism::theme_prism(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.caption = element_text(face = "italic", size = 10, hjust = 0),
      legend.position = "none"
    )
  
  if (!is.null(x_limits)) p <- p + coord_cartesian(xlim = x_limits)
  p
}

make_binary_plot <- function(data, xvar, xlab, ylab = "Domain PC1 score",
                             title = NULL, var_levels = NULL,
                             var_labels = NULL, palette = NULL) {
  dat <- data %>% filter(!is.na(.data[[xvar]]), !is.na(Domain_PC1))
  
  if (!is.null(var_levels)) {
    dat[[xvar]] <- factor(dat[[xvar]], levels = var_levels, labels = var_labels)
  } else {
    dat[[xvar]] <- factor(dat[[xvar]])
  }
  
  pval_df <- dat %>%
    group_by(Domain) %>%
    group_modify(~ {
      dd <- .x
      grp <- droplevels(factor(dd[[xvar]]))
      if (nlevels(grp) < 2) {
        return(tibble(label = "n<2", x = 1.5, y = max(dd$Domain_PC1, na.rm = TRUE) * 1.05))
      }
      p <- tryCatch(wilcox.test(dd$Domain_PC1 ~ grp)$p.value, error = function(e) NA_real_)
      tibble(
        label = paste0("p = ", signif(p, 3)),
        x = 1.5,
        y = max(dd$Domain_PC1, na.rm = TRUE) * 1.05
      )
    }) %>%
    ungroup()
  
  p <- ggplot(dat, aes(x = .data[[xvar]], y = Domain_PC1, fill = Domain, color = Domain)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.8, color = "black") +
    ggbeeswarm::geom_beeswarm(alpha = 0.5, size = 2.8, shape = 21, stroke = 0.4, dodge.width = 0.6, cex = 1.4) +
    coord_flip() +
    geom_text(
      data = pval_df,
      aes(x = x, y = y, label = label, color = Domain),
      inherit.aes = FALSE,
      size = 4.5,
      fontface = "bold"
    ) +
    facet_wrap(~ Domain, scales = "free_y", nrow = 1) +
    labs(
      x = xlab,
      y = ylab,
      title = title) +
    ggprism::theme_prism(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.caption = element_text(face = "italic", size = 10, hjust = 0),
      legend.position = "none"
    )
  
  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette) + scale_color_manual(values = palette)
  } else {
    p <- p + scale_fill_manual(values = domain_colors) + scale_color_manual(values = domain_colors)
  }
  p
}

# ------------------------------------------------------------
# Adjusted pathogen-status associations with Domain PC1
# ------------------------------------------------------------

run_adjusted_pathogen_model <- function(data, exposure, adjust_covariates, exposure_label) {
  needed <- unique(c("Domain", "Domain_PC1", exposure, adjust_covariates))

  dat <- data %>%
    dplyr::mutate(
      Domain_PC1 = suppressWarnings(as.numeric(Domain_PC1)),
      age = suppressWarnings(as.numeric(age)),
      illnessduration_enroll = suppressWarnings(as.numeric(illnessduration_enroll)),
      sex = factor(sex),
      hivrdtresult = factor(as.character(hivrdtresult), levels = c("0", "1")),
      malariardtresult = factor(as.character(malariardtresult), levels = c("0", "1")),
      microtbdx = factor(as.character(microtbdx), levels = c("0", "1"))
    ) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(needed), ~ !is.na(.x)))

  model_formula <- reformulate(c(exposure, adjust_covariates), response = "Domain_PC1")
  exposure_term <- paste0(exposure, "1")

  dat %>%
    dplyr::group_by(Domain) %>%
    dplyr::group_modify(~ {
      dd <- .x

      if (nrow(dd) < 10 || dplyr::n_distinct(dd[[exposure]]) < 2) {
        return(tibble::tibble(
          pathogen = exposure_label,
          exposure = exposure,
          beta = NA_real_,
          conf.low = NA_real_,
          conf.high = NA_real_,
          p.value = NA_real_,
          N = nrow(dd),
          covariates = paste(adjust_covariates, collapse = ", ")
        ))
      }

      fit <- tryCatch(
        stats::lm(model_formula, data = dd),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        return(tibble::tibble(
          pathogen = exposure_label,
          exposure = exposure,
          beta = NA_real_,
          conf.low = NA_real_,
          conf.high = NA_real_,
          p.value = NA_real_,
          N = nrow(dd),
          covariates = paste(adjust_covariates, collapse = ", ")
        ))
      }

      td <- broom::tidy(fit, conf.int = TRUE)
      row <- td %>% dplyr::filter(term == exposure_term)

      if (nrow(row) == 0) {
        return(tibble::tibble(
          pathogen = exposure_label,
          exposure = exposure,
          beta = NA_real_,
          conf.low = NA_real_,
          conf.high = NA_real_,
          p.value = NA_real_,
          N = stats::nobs(fit),
          covariates = paste(adjust_covariates, collapse = ", ")
        ))
      }

      tibble::tibble(
        pathogen = exposure_label,
        exposure = exposure,
        beta = row$estimate[1],
        conf.low = row$conf.low[1],
        conf.high = row$conf.high[1],
        p.value = row$p.value[1],
        N = stats::nobs(fit),
        covariates = paste(adjust_covariates, collapse = ", ")
      )
    }) %>%
    dplyr::ungroup()
}


make_pathogen_boxplot <- function(data, xvar, title,
                                  negative_label, positive_label,
                                  palette) {
  dat <- data %>%
    dplyr::mutate(
      "{xvar}" := factor(
        as.character(.data[[xvar]]),
        levels = c("0", "1"),
        labels = c(negative_label, positive_label)
      )
    ) %>%
    dplyr::filter(!is.na(.data[[xvar]]), !is.na(Domain_PC1))

  ggplot(
    dat,
    aes(
      x = .data[[xvar]],
      y = Domain_PC1,
      fill = .data[[xvar]],
      color = .data[[xvar]]
    )
  ) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      linewidth = 0.8,
      color = "black"
    ) +
    ggbeeswarm::geom_beeswarm(
      alpha = 0.5,
      size = 2.8,
      shape = 21,
      stroke = 0.4,
      dodge.width = 0.6,
      cex = 1.4
    ) +
    coord_flip() +
    facet_wrap(~ Domain, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    labs(
      x = "",
      y = "Domain PC1 score",
      title = title
    ) +
    ggprism::theme_prism(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
}


make_domain_pathogen_forest <- function(results, domain_name) {

  # Pathogen colors use slightly darker versions of the pathogen-density palette:
  # HIV = green, Malaria = blue, MTB = pink.
  plot_df <- results %>%
    dplyr::filter(
      as.character(Domain) == domain_name,
      pathogen %in% c("HIV", "Malaria", "MTB")
    ) %>%
    dplyr::mutate(
      pathogen = factor(
        pathogen,
        levels = c("HIV", "Malaria", "MTB")
      )
    )

  finite_ci <- c(plot_df$conf.low, plot_df$conf.high)
  finite_ci <- finite_ci[is.finite(finite_ci)]

  x_span <- if (length(finite_ci) >= 2) {
    diff(range(finite_ci))
  } else {
    1
  }

  if (!is.finite(x_span) || x_span <= 0) x_span <- 1

  plot_df <- plot_df %>%
    dplyr::mutate(
      label = dplyr::case_when(
        is.na(beta) ~ NA_character_,
        TRUE ~ paste0(
          "\u03b2 = ", signif(beta, 3),
          " (", signif(conf.low, 3), " to ", signif(conf.high, 3), ")",
          "\np = ", signif(p.value, 3)
        )
      ),
      label_x = dplyr::if_else(
        is.na(conf.high),
        NA_real_,
        conf.high + 0.06 * x_span
      )
    )

  pathogen_forest_colors <- c(
    "HIV" = "#3A9D82",
    "Malaria" = "#627FB8",
    "MTB" = "#C85D9F"
  )

  ggplot(plot_df, aes(x = beta, y = pathogen, color = pathogen)) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "gray40"
    ) +
    geom_point(
      size = 4.5,
      na.rm = TRUE
    ) +
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      orientation = "y",
      width = 0.15,
      linewidth = 1,
      na.rm = TRUE
    ) +
    geom_text(
      aes(x = label_x, label = label, color = pathogen),
      hjust = 0,
      vjust = 0.5,
      size = 3.7,
      fontface = "bold",
      na.rm = TRUE
    ) +
    scale_color_manual(values = pathogen_forest_colors) +
    scale_x_continuous(
      expand = expansion(mult = c(0.10, 0.42))
    ) +
    labs(
      x = "Adjusted \u03b2 (positive vs negative)",
      y = "",
      title = paste0(
        "Adjusted pathogen associations with ",
        domain_name,
        " PC1"
      )) +
    ggprism::theme_prism(base_size = 14) +
    theme(
      axis.text.y = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.caption = element_text(face = "italic", size = 10, hjust = 0),
      legend.position = "none"
    )
}


make_pathogen_density <- function(data) {
  # Influenza is retained for descriptive density visualization only;
  # it is excluded from all inferential regression models.
  if (!"influenzapcrresult" %in% names(data)) {
    stop("Influenza density plot requires clinical column: influenzapcrresult")
  }
  df_long <- data %>%
    dplyr::mutate(
      dplyr::across(
        c(influenzapcrresult, microtbdx, hivrdtresult, malariardtresult),
        ~ as.integer(as.character(.x))
      )
    ) %>%
    dplyr::select(
      pid, Domain, Domain_PC1,
      influenzapcrresult, microtbdx, hivrdtresult, malariardtresult
    ) %>%
    tidyr::pivot_longer(
      cols = c(influenzapcrresult, microtbdx, hivrdtresult, malariardtresult),
      names_to = "pathogen",
      values_to = "status"
    ) %>%
    dplyr::filter(status == 1) %>%
    dplyr::mutate(
      pathogen = dplyr::recode(
        pathogen,
        influenzapcrresult = "Influenza",
        microtbdx = "MTB",
        hivrdtresult = "HIV",
        malariardtresult = "Malaria"
      ),
      pathogen = factor(
        pathogen,
        levels = c("HIV", "Influenza", "Malaria", "MTB")
      )
    )

  pathogen_cols <- c(
    "HIV" = "#66C2A5",
    "Influenza" = "#FC8D62",
    "Malaria" = "#8DA0CB",
    "MTB" = "#E78AC3"
  )

  ggplot(df_long, aes(x = Domain_PC1, fill = pathogen, color = pathogen)) +
    geom_density(alpha = 0.35, linewidth = 1) +
    facet_wrap(~ Domain, nrow = 1, scales = "free") +
    scale_fill_manual(values = pathogen_cols) +
    scale_color_manual(values = pathogen_cols) +
    labs(
      x = "Domain PC1 score",
      y = "Density",
      fill = "Pathogen",
      color = "Pathogen",
      title = "Domain PC1 distribution by pathogen positivity"
    ) +
    ggprism::theme_prism(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}


# ------------------------------------------------------------
# 5. Generate additional single-domain descriptive plots
# ------------------------------------------------------------

# Continuous
p_age        <- make_continuous_plot(reserve_long, "age", "Age (years)", title = "Domain PC1 vs age")
safe_print_plot(p_age, "p_age")
save_plot(
  p_age,
  "Entebbe_DRT_vs_age.pdf")

p_illnessdur <- make_continuous_plot(reserve_long, "illnessduration_enroll", "Pre-enrollment illness duration (days)", title = "Domain PC1 vs illness duration")
safe_print_plot(p_illnessdur, "p_illnessdur")
save_plot(
  p_illnessdur,
  "Entebbe_DRT_vs_illness_duration.pdf")



# Binary / categorical
# Sex remains descriptive/unadjusted.
p_sex <- make_binary_plot(
  reserve_long %>% filter(!is.na(sex)),
  "sex",
  "Sex",
  title = "Domain PC1 by sex"
)
safe_print_plot(p_sex, "p_sex")
save_plot(
  p_sex,
  "Entebbe_DRT_by_sex.pdf")

# Pathogen-status analyses:
# 1) Raw/descriptive boxplots with no adjusted-model annotation.
# 2) Separate domain-specific adjusted linear-regression forest plots with HIV, malaria, and MTB as rows.
#
# Linear-regression outcome = Domain PC1 score.
# Primary exposure = pathogen status.
# Covariates = age, sex, illness duration, plus the specified other major pathogen groups.

# ---------------------------
# HIV
# ---------------------------
p_hiv <- make_pathogen_boxplot(
  data = reserve_long,
  xvar = "hivrdtresult",
  title = "Domain PC1 by HIV status",
  negative_label = "HIV negative",
  positive_label = "PLWH",
  palette = c("HIV negative" = "#9ECAE1", "PLWH" = "#115f9a")
)
safe_print_plot(p_hiv, "p_hiv")
save_plot(
  p_hiv,
  "Entebbe_DRT_by_HIV_status.pdf")


# ---------------------------
# Malaria
# ---------------------------
p_malaria <- make_pathogen_boxplot(
  data = reserve_long,
  xvar = "malariardtresult",
  title = "Domain PC1 by malaria status",
  negative_label = "Malaria negative",
  positive_label = "Malaria positive",
  palette = c("Malaria negative" = "#BCBDDC", "Malaria positive" = "#54278F")
)
safe_print_plot(p_malaria, "p_malaria")
save_plot(
  p_malaria,
  "Entebbe_DRT_by_malaria_status.pdf")


# ---------------------------
# MTB
# ---------------------------
p_tb <- make_pathogen_boxplot(
  data = reserve_long,
  xvar = "microtbdx",
  title = "Domain PC1 by MTB status",
  negative_label = "MTB negative",
  positive_label = "MTB positive",
  palette = c("MTB negative" = "#C7E9C0", "MTB positive" = "#006D2C")
)
safe_print_plot(p_tb, "p_tb")
save_plot(
  p_tb,
  "Entebbe_DRT_by_MTB_status.pdf")





# Export adjusted pathogen-PC1 regression summaries (CSV only).
adjusted_pathogen_results <- dplyr::bind_rows(
  run_adjusted_pathogen_model(
    reserve_long,
    exposure = "hivrdtresult",
    adjust_covariates = c(
      "age", "sex", "illnessduration_enroll",
      "malariardtresult", "microtbdx"
    ),
    exposure_label = "HIV"
  ),
  run_adjusted_pathogen_model(
    reserve_long,
    exposure = "malariardtresult",
    adjust_covariates = c(
      "age", "sex", "illnessduration_enroll",
      "hivrdtresult", "microtbdx"
    ),
    exposure_label = "Malaria"
  ),
  run_adjusted_pathogen_model(
    reserve_long,
    exposure = "microtbdx",
    adjust_covariates = c(
      "age", "sex", "illnessduration_enroll",
      "hivrdtresult", "malariardtresult"
    ),
    exposure_label = "MTB"
  )
)

write.csv(
  adjusted_pathogen_results,
  file.path(out_dir, "Entebbe_domain_PC1_adjusted_pathogen_linear_regression.csv"),
  row.names = FALSE
)


# Forest plots for adjusted pathogen -> PC1 associations.
# Each DOMAIN gets its own forest plot; HIV, malaria, and MTB
# are shown as separate rows.

p_pathogen_forest_resistance <- make_domain_pathogen_forest(
  adjusted_pathogen_results,
  "Resistance"
)
safe_print_plot(p_pathogen_forest_resistance, "p_pathogen_forest_resistance")
save_plot(
  p_pathogen_forest_resistance,
  "Entebbe_Resistance_adjusted_pathogen_associations_forest.pdf")

p_pathogen_forest_tolerance <- make_domain_pathogen_forest(
  adjusted_pathogen_results,
  "Tolerance"
)
safe_print_plot(p_pathogen_forest_tolerance, "p_pathogen_forest_tolerance")
save_plot(
  p_pathogen_forest_tolerance,
  "Entebbe_Tolerance_adjusted_pathogen_associations_forest.pdf")

p_pathogen_forest_damage <- make_domain_pathogen_forest(
  adjusted_pathogen_results,
  "Damage"
)
safe_print_plot(p_pathogen_forest_damage, "p_pathogen_forest_damage")
save_plot(
  p_pathogen_forest_damage,
  "Entebbe_Damage_adjusted_pathogen_associations_forest.pdf")

# Additional descriptive pathogen-density plot
p_pathogen_density <- make_pathogen_density(reserve_long)
safe_print_plot(p_pathogen_density, "p_pathogen_density")
save_plot(
  p_pathogen_density,
  "Entebbe_DRT_pathogen_density.pdf")


# Joint multi-omic PC1 and transcriptomic projection outcome models are retained above.

cat("\nPlot PNG files written to:\n", out_dir, "\n")
