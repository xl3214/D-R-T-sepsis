# =============================================================================
# ACESO-Ghana: D/R/T domain clinical and longitudinal analyses
# Script: 06_ACESO_Ghana_DRT_Clinical.R
# =============================================================================
#
# INPUT
#   Updated input schema: clinical SAMPLE_ID vs. scores sample_id (lowercase).
#   Both files use SUBJECT_ID, COLL_TIME and COHORT = S/H.
#   Previously exported group = Sepsis/Healthy also remains supported.
#
# ANALYSES
#   (1) Baseline (00h): sepsis vs healthy controls for Resistance, Tolerance,
#       and Damage using descriptive summaries, boxplots with jittered points,
#       and Wilcoxon rank-sum p-values (NO adjusted sepsis-vs-control models).
#   (2) Baseline sepsis: ordinal MEWS, with all 3 domain scores in the SAME
#       proportional-odds model, adjusted for AGE and SEX only.
#   (3) Sepsis: 28-day mortality modeled by joint baseline domain logistic
#       regression, with all three scores and age/sex in the same model.
#   (4) Sepsis: domain trajectories by eventual 28-day survival over 00-72h,
#       using categorical time, time x survival, age and sex, and a participant
#       random intercept. Report the global time x outcome interaction and
#       export within-time adjusted Died - Survived contrasts. No baseline
#       MEWS adjustment; no restricted 00-24h model.
#   (5) 72-hour LANDMARK: participants with paired 00h and 72h molecular scores.
#       Plot absolute changes (72h - 00h) as horizontal boxplots with
#       jittered observations and Wilcoxon tests;
#       fit joint logistic regression with all three changes and all three
#       enrollment scores, age, and sex. Optional Firth sensitivity.
#   No GEE models.
#
# SCALE
#   Original score scales are PRESERVED; the supplied wide score file has
#   normalized each domain relative to sepsis enrollment samples (mean 0, SD 1).
#   Do NOT restandardize separately within the healthy controls.
#
# IMPORTANT
#   No death times are supplied. Longitudinal later measurements
#   represent participants alive and sampled at those times; these are not
#   unbiased trajectories after death. The 72h landmark analysis conditions on
#   being alive AND sampled at 72h, so estimates apply to that selected group.
# =============================================================================

rm(list = ls())

required_pkgs <- c("dplyr", "tidyr", "purrr", "tibble", "ggplot2", "ggprism",
                   "MASS", "lme4", "emmeans")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                     logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install required packages before running: install.packages(c(",
       paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggprism)
})

# Optional Firth logistic sensitivity. Install once: install.packages("logistf")
RUN_FIRTH_SENSITIVITY <- requireNamespace("logistf", quietly = TRUE)

# ------------------------------------------------------------
# INPUT / OUTPUT PATHS (relative to the project root; edit if needed).
# Scores are produced by publication_code/01_DRT_scores.R.
# ------------------------------------------------------------
score_file    <- "DRT_scores_output/Ghana/Ghana_DRT_scores_wide_perSD.csv"
clinical_file <- "Analyses_Sept_2026/ACESO_Ghana/ACESO_GHANA_COGAPSPUB_PMID38816375_N526FUA_METACLINICAL_18DEC25.csv"

# Edit this directory if you prefer a different output location.
out_dir <- "DRT_scores_output/Clinical analyses/Ghana"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Could not create output directory: ", out_dir)

save_plot <- function(p, stem) {
  # Journal figures: vector PDF only, exactly 8.59 × 6.54 in (landscape).
  p <- p + ggplot2::labs(caption = NULL) +
    ggplot2::theme(plot.caption = ggplot2::element_blank())
  tryCatch(print(p), error = function(e)
    warning("Interactive preview failed for ", stem, ": ", conditionMessage(e)))
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(stem, ".pdf")), plot = p,
    device = "pdf", width = 8.59, height = 6.54, units = "in", bg = "white")
}

write_output <- function(data, basename) {
  utils::write.csv(data, file.path(out_dir, paste0(basename, ".csv")),
                   row.names = FALSE, na = "")
}
# No arbitrary reporting floor: e.g., p = 2.13e-09, not p < 0.0001.
# Values are calculated without rounding in all model result CSV files.
format_p <- function(p) {
  vapply(p, function(value) {
    if (is.na(value)) return("p = NA")
    if (value == 0) return("p = 0 (numerical underflow)")
    paste0("p = ", sprintf("%.3g", value))
  }, character(1))
}

# ----------------------------------------------------------------------------
# A. Read, harmonize COHORT/group codes, and merge by unique SAMPLE_ID
# ----------------------------------------------------------------------------
scores <- utils::read.csv(score_file, check.names = FALSE,
                          stringsAsFactors = FALSE)
clinical <- utils::read.csv(clinical_file, check.names = FALSE,
                            stringsAsFactors = FALSE)

# The refreshed molecular export uses lowercase 'sample_id'; the clinical
# metadata uses uppercase 'SAMPLE_ID'. Map each to the same canonical key.
# Reject conflicting duplicate columns instead of joining on a wrong ID.
harmonize_sample_id <- function(dat, source_name) {
  if (!"SAMPLE_ID" %in% names(dat) && !"sample_id" %in% names(dat)) {
    stop(source_name, " lacks SAMPLE_ID/sample_id for merging.")
  }
  if ("SAMPLE_ID" %in% names(dat) && "sample_id" %in% names(dat)) {
    upper <- trimws(as.character(dat$SAMPLE_ID))
    lower <- trimws(as.character(dat$sample_id))
    if (!identical(upper, lower)) {
      stop(source_name, ": SAMPLE_ID and sample_id columns conflict.")
    }
  }
  if (!"SAMPLE_ID" %in% names(dat)) dat$SAMPLE_ID <- dat$sample_id
  dat$SAMPLE_ID <- trimws(as.character(dat$SAMPLE_ID))
  if (anyNA(dat$SAMPLE_ID) || any(!nzchar(dat$SAMPLE_ID))) {
    stop(source_name, ": missing/blank sample identifiers.")
  }
  if (anyDuplicated(dat$SAMPLE_ID)) {
    stop(source_name, ": duplicate sample identifiers.")
  }
  dat
}
scores <- harmonize_sample_id(scores, "Ghana molecular scores")
clinical <- harmonize_sample_id(clinical, "Ghana clinical metadata")

methods <- tibble::tribble(
  ~method,            ~prefix,
  "Transcriptomic projection",        "rna_surrogate_lasso14_1se__"
)
domains <- c("Resistance", "Tolerance", "Damage")
score_lookup <- tidyr::expand_grid(method = methods$method, domain = domains) %>%
  dplyr::left_join(methods, by = "method") %>%
  dplyr::mutate(score_column = paste0(prefix, domain))

require_columns <- function(x, cols, dataset_name) {
  absent <- setdiff(cols, names(x))
  if (length(absent)) stop(dataset_name, " lacks: ", paste(absent, collapse = ", "))
}
# The refreshed Ghana score file contains BOTH a descriptive 'cohort'
# column ("Ghana") and the categorical 'COHORT' column (S/H). Only the
# uppercase COHORT identifies sepsis vs. healthy controls. Previously
# exported 'group' = Sepsis/Healthy is also supported.
# Normalize these labels once so downstream analyses remain unchanged.
require_columns(scores,
                c("SAMPLE_ID", "SUBJECT_ID", "COLL_TIME",
                  score_lookup$score_column), "Ghana scores file")
if (!"COHORT" %in% names(scores) && !"group" %in% names(scores)) {
  stop("Ghana scores file requires COHORT (S/H) or group (Sepsis/Healthy).")
}
require_columns(clinical,
                c("SAMPLE_ID", "SUBJECT_ID", "COLL_TIME", "Age", "Sex",
                  "Survival (28days)", "MEWS"), "Ghana clinical file")

normalize_group <- function(x, source_name) {
  x <- trimws(as.character(x))
  mapped <- dplyr::case_when(
    x %in% c("S", "Sepsis") ~ "Sepsis",
    x %in% c("H", "Healthy") ~ "Healthy",
    TRUE ~ NA_character_
  )
  if (anyNA(mapped)) {
    stop(source_name, " has unexpected/missing group codes: ",
         paste(unique(x[is.na(mapped)]), collapse = ", "))
  }
  mapped
}

scores$group_harmonized <- if ("COHORT" %in% names(scores)) {
  normalize_group(scores$COHORT, "Ghana scores COHORT")
} else {
  normalize_group(scores$group, "Ghana scores group")
}
# Ensure consistency if a future score export contains BOTH forms.
if ("COHORT" %in% names(scores) && "group" %in% names(scores)) {
  if (any(scores$group_harmonized != normalize_group(scores$group, "Ghana scores group"))) {
    stop("Ghana scores file has conflicting COHORT and group assignments")
  }
}
# Cross-check the cohort assignment against the clinical file whenever present.
clinical$group_clinical <- if ("COHORT" %in% names(clinical)) {
  normalize_group(clinical$COHORT, "Ghana clinical COHORT")
} else {
  rep(NA_character_, nrow(clinical))
}

# The key is already harmonized/validated above.
# Verify that all molecular IDs match a metadata sample, and vice versa.
n_scores_unmatched <- sum(!scores$SAMPLE_ID %in% clinical$SAMPLE_ID)
n_clinical_unmatched <- sum(!clinical$SAMPLE_ID %in% scores$SAMPLE_ID)
cat("Ghana input ID matching: ", nrow(scores), " molecular, ",
    nrow(clinical), " clinical, ",
    n_scores_unmatched, " molecular IDs unmatched, ",
    n_clinical_unmatched, " clinical IDs unmatched.\n", sep = "")

score_core <- scores %>%
  dplyr::transmute(SAMPLE_ID, subject = trimws(as.character(SUBJECT_ID)),
    time = trimws(as.character(COLL_TIME)),
    group = group_harmonized,
    dplyr::across(dplyr::all_of(score_lookup$score_column),
                  ~ suppressWarnings(as.numeric(.x))))
clinical_core <- clinical %>%
  dplyr::transmute(SAMPLE_ID,
    subject_clinical = trimws(as.character(SUBJECT_ID)),
    time_clinical = trimws(as.character(COLL_TIME)),
    group_clinical,
    age = suppressWarnings(as.numeric(Age)),
    sex = factor(as.character(Sex), levels = c("Male", "Female")),
    mortality28 = as.character(.data[["Survival (28days)"]]),
    mews = suppressWarnings(as.numeric(as.character(MEWS))))

merged <- dplyr::inner_join(score_core, clinical_core,
                            by = "SAMPLE_ID")
if (any(!is.na(merged$group_clinical) &
        merged$group != merged$group_clinical)) {
  stop("Scores and clinical files disagree on Sepsis/Healthy status for some SAMPLE_IDs")
}
if (!nrow(merged)) stop("No SAMPLE_ID matches between the input files")
if (nrow(merged) < nrow(scores) || nrow(merged) < nrow(clinical)) {
  warning("Not all source samples matched; check sample availability outputs.")
}
if (any(merged$subject != merged$subject_clinical) ||
    any(merged$time != merged$time_clinical)) {
  stop("Subject or timepoint mismatch between clinical and scores files")
}
if (anyDuplicated(merged[, c("subject", "time")])) {
  stop("More than one sample for the same subject and collection time")
}
if (!all(merged$group %in% c("Healthy", "Sepsis"))) {
  stop("Unexpected group values: ", paste(unique(merged$group), collapse = ", "))
}

score_long <- merged %>%
  dplyr::select(SAMPLE_ID, subject, time, group, age, sex,
                mortality28, mews, dplyr::all_of(score_lookup$score_column)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(score_lookup$score_column),
                      names_to = "score_column", values_to = "score") %>%
  dplyr::left_join(dplyr::select(score_lookup, method, domain, score_column),
                   by = "score_column") %>%
  dplyr::mutate(
    method = factor(method, levels = methods$method),
    domain = factor(domain, levels = domains)
  )

baseline <- score_long %>%
  dplyr::filter(time == "00h", group %in% c("Healthy", "Sepsis")) %>%
  dplyr::mutate(group = factor(group, levels = c("Healthy", "Sepsis")))
if (dplyr::n_distinct(baseline$subject) != 130L) {
  warning("Expected 130 enrollment participants (102 sepsis, 28 controls); found ",
          dplyr::n_distinct(baseline$subject), ". Check sample QC.")
}

sample_qc <- merged %>%
  dplyr::group_by(group, time) %>%
  dplyr::summarise(n_samples = dplyr::n(), n_subjects = dplyr::n_distinct(subject),
                   n_mews = sum(!is.na(mews)), .groups = "drop")
write_output(sample_qc, "Ghana_DRT_sample_availability_QC")
write_output(score_lookup, "Ghana_DRT_score_method_key")
cat("\nInput: ", nrow(merged), " matched samples / ",
    dplyr::n_distinct(merged$subject), " participants.\n", sep = "")

# ----------------------------------------------------------------------------
# 1. Baseline domain: sepsis vs healthy control
#    Primary here is descriptive/unadjusted: boxplots + jitter + Wilcoxon test.
# ----------------------------------------------------------------------------
compare_baseline <- function(one_method, one_domain) {
  dd <- baseline %>%
    dplyr::filter(method == one_method, domain == one_domain,
                  is.finite(score), !is.na(group)) %>%
    droplevels()
  if (dplyr::n_distinct(dd$group) != 2L) stop("Missing baseline group for ", one_method, " / ", one_domain)
  g_s <- dd$score[dd$group == "Sepsis"]
  g_h <- dd$score[dd$group == "Healthy"]
  wt <- stats::wilcox.test(g_s, g_h, exact = FALSE)
  tibble::tibble(
    method = one_method, domain = one_domain,
    n_sepsis = length(g_s), n_healthy = length(g_h),
    mean_sepsis = mean(g_s), sd_sepsis = stats::sd(g_s),
    mean_healthy = mean(g_h), sd_healthy = stats::sd(g_h),
    median_sepsis = stats::median(g_s), median_healthy = stats::median(g_h),
    q1_sepsis = stats::quantile(g_s, probs = 0.25, na.rm = TRUE),
    q3_sepsis = stats::quantile(g_s, probs = 0.75, na.rm = TRUE),
    q1_healthy = stats::quantile(g_h, probs = 0.25, na.rm = TRUE),
    q3_healthy = stats::quantile(g_h, probs = 0.75, na.rm = TRUE),
    wilcoxon_W = unname(wt$statistic),
    p_wilcoxon = wt$p.value
  )
}

baseline_table <- purrr::map_dfr(seq_len(nrow(score_lookup)), function(i) {
  compare_baseline(score_lookup$method[i], score_lookup$domain[i])
}) %>%
  dplyr::group_by(method) %>%
  dplyr::mutate(q_BH_3domains = stats::p.adjust(p_wilcoxon, method = "BH")) %>%
  dplyr::ungroup()
write_output(baseline_table, "Ghana_DRT_baseline_sepsis_vs_healthy_wilcoxon")


baseline_plot_data <- baseline %>%
  dplyr::filter(method == "Transcriptomic projection", is.finite(score))
baseline_ann <- baseline_plot_data %>%
  dplyr::group_by(domain) %>%
  dplyr::summarise(y_max = max(score), y_min = min(score), .groups = "drop") %>%
  dplyr::left_join(
    baseline_table %>%
      dplyr::filter(method == "Transcriptomic projection") %>%
      dplyr::select(domain, p_wilcoxon),
    by = "domain"
  ) %>%
  dplyr::mutate(
    y = y_max + 0.10 * pmax(y_max - y_min, 0.5),
    label = format_p(p_wilcoxon), x = 1.5
  )

p_baseline <- ggplot2::ggplot(baseline_plot_data,
                              ggplot2::aes(x = group, y = score,
                                           fill = group, color = group)) +
  ggplot2::geom_boxplot(width = 0.57, outlier.shape = NA, alpha = 0.55,
                        color = "black", linewidth = 0.55) +
  ggplot2::geom_jitter(width = 0.12, height = 0, alpha = 0.55, size = 1.5) +
  ggplot2::geom_text(data = baseline_ann,
                     ggplot2::aes(x = x, y = y, label = label),
                     inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  ggplot2::facet_wrap(~ domain, nrow = 1, scales = "fixed") +
  ggplot2::scale_fill_manual(values = c(Healthy = "#BDBDBD", Sepsis = "#477FB9")) +
  ggplot2::scale_color_manual(values = c(Healthy = "#777777", Sepsis = "#286AB0")) +
  ggplot2::labs(title = "Baseline transcriptomic projection domain scores: sepsis vs healthy",
                x = NULL, y = "Domain score at baseline (SD units)") +
  ggprism::theme_prism(base_size = 13) +
  ggplot2::theme(legend.position = "none",
                 strip.text = ggplot2::element_text(face = "bold"),
                 plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                 axis.text.x = ggplot2::element_text(face = "bold"))
save_plot(p_baseline, "Ghana_RNA_LASSO_baseline_sepsis_healthy_boxplots_Wilcoxon")

# ----------------------------------------------------------------------------
# 2. Baseline MEWS: joint domain proportional-odds model, age and sex adjusted
#    No MEWS covariate: MEWS IS the outcome.
# ----------------------------------------------------------------------------
fit_baseline_mews <- function(one_method) {
  dd <- baseline %>%
    dplyr::filter(group == "Sepsis", method == one_method) %>%
    dplyr::select(subject, age, sex, mews, domain, score) %>%
    tidyr::pivot_wider(names_from = domain, values_from = score) %>%
    dplyr::filter(stats::complete.cases(dplyr::across(
      dplyr::all_of(c(domains, "age", "sex", "mews"))))) %>%
    droplevels()
  if (anyDuplicated(dd$subject)) stop("Repeated subject in baseline MEWS data")
  if (dplyr::n_distinct(dd$mews) < 3L) stop("Insufficient MEWS levels")
  dd$MEWS_ordinal <- ordered(dd$mews, levels = sort(unique(dd$mews)))
  fit <- MASS::polr(MEWS_ordinal ~ Resistance + Tolerance + Damage + age + sex,
                     data = dd, Hess = TRUE, method = "logistic")
  beta <- stats::coef(fit)[domains]
  se <- sqrt(diag(stats::vcov(fit)))[domains]
  if (any(!is.finite(beta)) || any(!is.finite(se))) stop("Invalid MEWS coefficients")
  output <- tibble::tibble(
    method = one_method, domain = domains, n = stats::nobs(fit),
    beta_logOR = unname(beta), SE = unname(se),
    OR = exp(unname(beta)),
    lower95 = exp(unname(beta - 1.96 * se)),
    upper95 = exp(unname(beta + 1.96 * se)),
    p = unname(2 * stats::pnorm(abs(beta / se), lower.tail = FALSE))
  ) %>% dplyr::mutate(q_BH_3domains = stats::p.adjust(p, method = "BH"))
  # Standard VIFs for the three molecular predictors, conditional on the
  # other molecular domains, age, and sex.
  vifs <- purrr::map_dfr(domains, function(term) {
    others <- setdiff(domains, term)
    aux <- stats::lm(stats::reformulate(c(others, "age", "sex"), response = term),
                     data = dd)
    tibble::tibble(method = one_method, domain = term,
                   n = stats::nobs(aux), VIF = 1 / (1 - summary(aux)$r.squared))
  })
  list(estimates = output, vifs = vifs)
}

mews_fits <- lapply(methods$method, fit_baseline_mews)
mews_table <- dplyr::bind_rows(lapply(mews_fits, `[[`, "estimates"))
mews_vifs <- dplyr::bind_rows(lapply(mews_fits, `[[`, "vifs"))
write_output(mews_table, "Ghana_DRT_baseline_joint_MEWS")
write_output(mews_vifs, "Ghana_DRT_baseline_joint_MEWS_domain_VIF")

# Display ORs, 95% CIs, and p-values from the SAME fitted joint MEWS model.
# The shared annotation x-coordinate keeps all three result labels aligned.
mews_plot_data <- mews_table %>%
  dplyr::filter(method == "Transcriptomic projection") %>%
  dplyr::mutate(
    domain = factor(domain, levels = rev(domains)),
    annotation_x = max(upper95, na.rm = TRUE) * 1.12,
    annotation = paste0(
      "OR = ", sprintf("%.2f", OR),
      " (", sprintf("%.2f", lower95), "–", sprintf("%.2f", upper95), ")",
      "\n", format_p(p)
    )
  )

p_mews <- ggplot2::ggplot(
  mews_plot_data, ggplot2::aes(x = OR, y = domain, color = domain)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = lower95, xmax = upper95),
    orientation = "y", linewidth = 0.9, width = 0.16
  ) +
  ggplot2::geom_point(size = 3.4) +
  ggplot2::geom_text(
    ggplot2::aes(x = annotation_x, label = annotation),
    hjust = 0, vjust = 0.5, size = 3.5, fontface = "bold",
    lineheight = 1.1, show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ method, nrow = 1) +
  ggplot2::scale_color_manual(values = c(
    Resistance = "#2A6EBB", Tolerance = "#2CA25F", Damage = "#E6550D"
  )) +
  ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.08, 0.65))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "Baseline MEWS: mutually adjusted domain associations",
    x = "Adjusted odds ratio for higher MEWS per 1-SD increase", y = NULL) +
  ggprism::theme_prism(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.margin = ggplot2::margin(10, 140, 10, 10)
  )
save_plot(p_mews, "Ghana_RNA_LASSO_baseline_joint_MEWS_forest")

# ----------------------------------------------------------------------------
# 3. Baseline transcriptomic projection domain association with 28-day mortality
#    Binary logistic model: 28-day death ~ Resistance + Tolerance + Damage
#                                 + age + sex
#    All three molecular scores enter SIMULTANEOUSLY (mutual adjustment).
#    Effect = adjusted OR for 28-day death per 1 SD higher baseline score.
# ----------------------------------------------------------------------------
fit_joint_baseline_mortality <- function(one_method) {
  dd <- baseline %>%
    dplyr::filter(group == "Sepsis", method == one_method,
                  mortality28 %in% c("Survived", "Died")) %>%
    dplyr::select(subject, age, sex, mortality28, domain, score) %>%
    tidyr::pivot_wider(names_from = domain, values_from = score) %>%
    dplyr::filter(stats::complete.cases(dplyr::across(
      dplyr::all_of(c(domains, "age", "sex", "mortality28"))))) %>%
    droplevels()
  if (anyDuplicated(dd$subject)) stop("Duplicate baseline subject for ", one_method)
  dd$death28 <- as.integer(dd$mortality28 == "Died")
  if (dplyr::n_distinct(dd$death28) != 2L) {
    stop("Both 28-day mortality outcome groups are needed for ", one_method)
  }

  fit <- stats::glm(
    death28 ~ Resistance + Tolerance + Damage + age + sex,
    data = dd, family = stats::binomial()
  )
  beta <- stats::coef(fit)[domains]
  se <- sqrt(diag(stats::vcov(fit)))[domains]
  if (any(!is.finite(beta)) || any(!is.finite(se)) || any(se <= 0)) {
    stop("Unable to estimate all joint mortality coefficients for ", one_method)
  }
  estimates <- tibble::tibble(
    method = one_method, domain = domains,
    n = stats::nobs(fit), deaths = sum(dd$death28 == 1L),
    survivors = sum(dd$death28 == 0L),
    log_OR = unname(beta), SE = unname(se),
    OR = exp(unname(beta)),
    lower95 = exp(unname(beta - 1.96 * se)),
    upper95 = exp(unname(beta + 1.96 * se)),
    p = unname(2 * stats::pnorm(abs(beta / se), lower.tail = FALSE))
  ) %>%
    dplyr::mutate(q_BH_3domains = stats::p.adjust(p, method = "BH"))

  vifs <- purrr::map_dfr(domains, function(term) {
    aux <- stats::lm(
      stats::reformulate(c(setdiff(domains, term), "age", "sex"),
                         response = term),
      data = dd
    )
    tibble::tibble(
      method = one_method, domain = term,
      n = stats::nobs(aux),
      VIF = 1 / (1 - summary(aux)$r.squared)
    )
  })
  list(estimates = estimates, vifs = vifs)
}

baseline_mortality_fits <- lapply(methods$method, fit_joint_baseline_mortality)
baseline_mortality_table <- dplyr::bind_rows(
  lapply(baseline_mortality_fits, `[[`, "estimates")
)
baseline_mortality_vifs <- dplyr::bind_rows(
  lapply(baseline_mortality_fits, `[[`, "vifs")
)
write_output(baseline_mortality_table,
             "Ghana_DRT_baseline_joint_28day_mortality_logistic")
write_output(baseline_mortality_vifs,
             "Ghana_DRT_baseline_joint_28day_mortality_VIF")

mortality_forest <- baseline_mortality_table %>%
  dplyr::filter(method == "Transcriptomic projection") %>%
  dplyr::mutate(
    domain = factor(domain, levels = rev(domains)),
    annotation = paste0(
      "OR = ", sprintf("%.2f", OR),
      " (", sprintf("%.2f", lower95), "–", sprintf("%.2f", upper95), ")",
      "\n", format_p(p)
    )
  )

p_baseline_mortality <- ggplot2::ggplot(
    mortality_forest, ggplot2::aes(x = OR, y = domain, color = domain)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = lower95, xmax = upper95),
    orientation = "y", width = 0.16, linewidth = 0.9
  ) +
  ggplot2::geom_point(size = 3.7) +
  ggplot2::geom_text(
    ggplot2::aes(x = upper95 * 1.06, label = annotation),
    hjust = 0, vjust = 0.5, size = 3.6, fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = c(
    Resistance = "#2A6EBB", Tolerance = "#2CA25F", Damage = "#E6550D"
  )) +
  ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.09, 0.80))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "Baseline transcriptomic projection domain associations with 28-day mortality",
    x = "Adjusted odds ratio for 28-day mortality per 1-SD increase",
    y = NULL) +
  ggprism::theme_prism(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    axis.text.y = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(10, 130, 10, 10)
  )
save_plot(p_baseline_mortality,
          "Ghana_RNA_LASSO_baseline_joint_28day_mortality_logistic_forest")

# ----------------------------------------------------------------------------
# 4. Longitudinal domain trajectories versus eventual 28-day survival
#    - 00h/06h/24h/48h/72h only (full window)
#    Restrict to sepsis subjects enrolled at 00h, with known 28-day outcome.
#    Use BASELINE age, sex and survival for all the subject's observations.
#    Random intercept accounts for within-person repeated measurements.
# ----------------------------------------------------------------------------
sepsis_at_enrollment <- merged %>%
  dplyr::filter(group == "Sepsis", time == "00h") %>%
  dplyr::select(subject, age, sex, mortality28) %>%
  dplyr::filter(mortality28 %in% c("Survived", "Died"),
                !is.na(age), !is.na(sex)) %>%
  dplyr::mutate(mortality28 = factor(mortality28,
                                    levels = c("Survived", "Died")))
if (anyDuplicated(sepsis_at_enrollment$subject)) stop("Duplicate enrollment records")

acute_times <- c("00h", "06h", "24h", "48h", "72h")
acute_long <- score_long %>%
  dplyr::filter(group == "Sepsis", time %in% acute_times) %>%
  dplyr::select(SAMPLE_ID, subject, time, method, domain, score) %>%
  dplyr::inner_join(sepsis_at_enrollment, by = "subject") %>%
  dplyr::filter(is.finite(score)) %>%
  dplyr::mutate(time = factor(time, levels = acute_times),
                mortality28 = factor(mortality28, levels = c("Survived", "Died")))
cat("\nLongitudinal baseline-eligible sepsis subjects: ",
    dplyr::n_distinct(acute_long$subject), "\n", sep = "")

# Sample availability is counted once per patient/time, not once per score.
followup_counts <- acute_long %>%
  dplyr::filter(method == "Transcriptomic projection", domain == "Resistance") %>%
  dplyr::group_by(time, mortality28) %>%
  dplyr::summarise(n_samples = dplyr::n(),
                   n_participants = dplyr::n_distinct(subject), .groups = "drop") %>%
  tidyr::complete(time = factor(acute_times, levels = acute_times),
                  mortality28 = factor(c("Survived", "Died"),
                                       levels = c("Survived", "Died")),
                  fill = list(n_samples = 0L, n_participants = 0L))
write_output(followup_counts, "Ghana_DRT_0to72h_availability_by_28d_survival")

# One scoring method / molecular domain, full 0-72h window only.
fit_trajectory <- function(one_method, one_domain) {
  valid_times <- acute_times
  dd <- acute_long %>%
    dplyr::filter(method == one_method, domain == one_domain,
                  as.character(time) %in% valid_times) %>%
    dplyr::mutate(time = factor(as.character(time), levels = valid_times)) %>%
    dplyr::arrange(subject, time) %>% droplevels()
  if (dplyr::n_distinct(dd$subject) < 20L ||
      dplyr::n_distinct(dd$time) < length(valid_times)) {
    stop("Insufficient longitudinal data: 0-72h / ", one_method, " / ", one_domain)
  }
  full_formula <- stats::as.formula(
    "score ~ time * mortality28 + age + sex + (1 | subject)")
  reduced_formula <- stats::as.formula(
    "score ~ time + mortality28 + age + sex + (1 | subject)")
  # Maximum likelihood for a valid nested-model likelihood-ratio comparison.
  full <- lme4::lmer(full_formula, data = dd, REML = FALSE,
                     control = lme4::lmerControl(optimizer = "bobyqa"))
  reduced <- lme4::lmer(reduced_formula, data = dd, REML = FALSE,
                        control = lme4::lmerControl(optimizer = "bobyqa"))
  lrt <- stats::anova(reduced, full)
  p_interaction <- unname(lrt[2, "Pr(>Chisq)"])

  # Marginal means at each categorical time, adjusted for age and sex.
  emm <- emmeans::emmeans(full, ~ time * mortality28,
                          lmer.df = "asymptotic", weights = "proportional")
  pred <- as.data.frame(stats::confint(emm))
  low_name <- grep("(lower.CL|asymp.LCL|LCL)$", names(pred), value = TRUE)[1]
  high_name <- grep("(upper.CL|asymp.UCL|UCL)$", names(pred), value = TRUE)[1]
  if (is.na(low_name) || is.na(high_name)) stop("Missing emmeans confidence limits")
  means <- tibble::tibble(
    window = "0-72h", method = one_method, domain = one_domain,
    time = as.character(pred$time),
    mortality28 = as.character(pred$mortality28),
    adjusted_mean = pred$emmean,
    lower95 = pred[[low_name]], upper95 = pred[[high_name]],
    SE = pred$SE
  )
  test <- tibble::tibble(
    window = "0-72h", method = one_method, domain = one_domain,
    n_participants = dplyr::n_distinct(dd$subject),
    n_observations = nrow(dd),
    n_eventual_deaths = dplyr::n_distinct(dd$subject[dd$mortality28 == "Died"]),
    global_time_x_mortality_chisq = unname(lrt[2, "Chisq"]),
    global_df = unname(lrt[2, "Df"]),
    p_global_time_x_mortality = p_interaction,
    singular_random_intercept = lme4::isSingular(full, tol = 1e-4)
  )
  # Optional time-specific contrasts describe adjusted Died - Survived
  # differences at each measured timepoint. These are descriptive because
  # later observations depend on survival and availability for sampling.
  get_contrast <- function(cgrid, contrast_name) {
    result <- as.data.frame(summary(
      emmeans::contrast(cgrid, method = list("Died - Survived" = c(-1, 1))),
      infer = c(TRUE, TRUE), adjust = "none"))
    lower <- grep("^(lower.CL|asymp.LCL|LCL)$", names(result), value = TRUE)[1]
    upper <- grep("^(upper.CL|asymp.UCL|UCL)$", names(result), value = TRUE)[1]
    if (is.na(lower) || is.na(upper))
      stop("Missing confidence limits for ", contrast_name)
    tibble::tibble(window = "0-72h", method = one_method, domain = one_domain,
      contrast_type = contrast_name,
      time = as.character(result$time),
      contrast = as.character(result$contrast), estimate = as.numeric(result$estimate),
      SE = as.numeric(result$SE), lower95 = as.numeric(result[[lower]]),
      upper95 = as.numeric(result[[upper]]), p = as.numeric(result$p.value),
      n_participants = dplyr::n_distinct(dd$subject))
  }
  time_emm <- emmeans::emmeans(full, ~ mortality28 | time,
                               lmer.df = "asymptotic")
  time_differences <- get_contrast(time_emm, "Time-specific Died minus Survived")
  list(test = test, means = means, time_differences = time_differences)
}

trajectory_grid <- tidyr::expand_grid(
  method = methods$method,
  domain = domains
)
trajectory_fits <- purrr::pmap(
  trajectory_grid,
  function(method, domain) fit_trajectory(method, domain)
)
trajectory_tests <- dplyr::bind_rows(lapply(trajectory_fits, `[[`, "test")) %>%
  dplyr::group_by(window, method) %>%
  dplyr::mutate(q_BH_3domains = stats::p.adjust(
    p_global_time_x_mortality, method = "BH")) %>%
  dplyr::ungroup()
trajectory_means <- dplyr::bind_rows(lapply(trajectory_fits, `[[`, "means"))
trajectory_time_contrasts <- dplyr::bind_rows(
  lapply(trajectory_fits, `[[`, "time_differences")) %>%
  dplyr::group_by(method, time) %>%
  dplyr::mutate(q_BH_3domains_at_time = stats::p.adjust(p, method = "BH")) %>%
  dplyr::ungroup()
write_output(trajectory_time_contrasts,
             "Ghana_DRT_longitudinal_time_specific_Died_minus_Survived_LME")
write_output(trajectory_tests, "Ghana_DRT_longitudinal_time_by_mortality_LME_tests")
write_output(trajectory_means, "Ghana_DRT_longitudinal_age_sex_adjusted_estimated_means")

# Transcriptomic projection trajectory figure only: full 0-72h, faceted by domain.
# COMMON y-axis across the 3 domains within each figure, plus raw jittered points.
for (one_window in "0-72h") {
  valid_times <- acute_times
  for (one_method in "Transcriptomic projection") {
    panel_key <- trajectory_tests %>%
      dplyr::filter(window == one_window, method == one_method) %>%
      dplyr::mutate(
        panel_label = paste0(domain, "\nTime × mortality ",
                             format_p(p_global_time_x_mortality))
      ) %>%
      dplyr::select(domain, panel_label)

    dd_mean <- trajectory_means %>%
      dplyr::filter(window == one_window, method == one_method) %>%
      dplyr::left_join(panel_key, by = "domain") %>%
      dplyr::mutate(
        time = factor(time, levels = valid_times),
        mortality28 = factor(mortality28, levels = c("Survived", "Died")),
        panel_label = factor(panel_label, levels = panel_key$panel_label)
      )

    dd_raw <- acute_long %>%
      dplyr::filter(method == one_method, domain %in% domains,
                    as.character(time) %in% valid_times) %>%
      dplyr::left_join(panel_key, by = "domain") %>%
      dplyr::mutate(
        time = factor(as.character(time), levels = valid_times),
        mortality28 = factor(mortality28, levels = c("Survived", "Died")),
        panel_label = factor(panel_label, levels = panel_key$panel_label)
      )

    y_lim <- range(c(dd_raw$score, dd_mean$lower95, dd_mean$upper95), na.rm = TRUE)
    y_pad <- 0.06 * diff(y_lim)
    if (!is.finite(y_pad) || y_pad <= 0) y_pad <- 0.2
    y_lim <- c(y_lim[1] - y_pad, y_lim[2] + y_pad)

    pd <- ggplot2::position_dodge(width = 0.28)
    pjd <- ggplot2::position_jitterdodge(jitter.width = 0.08, jitter.height = 0,
                                         dodge.width = 0.28, seed = 123)

    p <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = dd_raw,
        ggplot2::aes(x = time, y = score, color = mortality28),
        position = pjd, alpha = 0.28, size = 1.2, stroke = 0
      ) +
      ggplot2::geom_line(
        data = dd_mean,
        ggplot2::aes(x = time, y = adjusted_mean, color = mortality28,
                     group = mortality28),
        position = pd, linewidth = 1.05
      ) +
      ggplot2::geom_point(
        data = dd_mean,
        ggplot2::aes(x = time, y = adjusted_mean, color = mortality28,
                     group = mortality28),
        position = pd, size = 2.7
      ) +
      ggplot2::geom_errorbar(
        data = dd_mean,
        ggplot2::aes(x = time, y = adjusted_mean, ymin = lower95, ymax = upper95,
                     color = mortality28, group = mortality28),
        position = pd, width = 0.10, linewidth = 0.70
      ) +
      ggplot2::facet_wrap(~panel_label, nrow = 1, scales = "fixed") +
      ggplot2::coord_cartesian(ylim = y_lim) +
      ggplot2::scale_color_manual(values = c(Survived = "#2A6EBB",
                                             Died = "#D65F44")) +
      ggplot2::labs(
        title = paste0("Ghana sepsis: ", one_method, " domain trajectories (", one_window, ")"),
        x = "Hours after enrollment", y = "Adjusted domain score (SD units)",
        color = "28-day outcome") +
      ggprism::theme_prism(base_size = 13) +
      ggplot2::theme(legend.position = "top",
                     strip.text = ggplot2::element_text(face = "bold", size = 11),
                     plot.title = ggplot2::element_text(face = "bold", hjust = .5),
                     plot.caption = ggplot2::element_text(size = 9, hjust = 0))

    file_stub <- paste0("Ghana_DRT_trajectories_commonY_jitter_",
                        gsub("[^A-Za-z0-9]+", "_", one_method), "_",
                        gsub("[^0-9A-Za-z]+", "", one_window))
    save_plot(p, file_stub)
  }
}


# ----------------------------------------------------------------------------
# 5. 72-hour fixed landmark: 00h-to-72h ABSOLUTE domain changes and mortality
# ----------------------------------------------------------------------------
# Include ONLY patients with both 00h and 72h samples and known 28-day status.
# A 72h sample implies survival to that timepoint. The outcome is subsequent
# death by day 28 in the selected 72h-alive-and-sampled group.
# This is NOT a first-to-last slope, and there is NO elapsed-time divisor.
# For each domain: delta = score(72h) - score(00h).
# Joint logistic model: death28 ~ delta_R + delta_T + delta_D +
#                               baseline_R + baseline_T + baseline_D + age + sex.
# The three change variables are standardized using the SD of changes within
# the paired complete-case sample. This is separate from original baseline
# scores, whose scale remains unchanged.
# Firth penalized logistic regression is an OPTIONAL small-sample sensitivity.
# ----------------------------------------------------------------------------

paired_input <- merged %>%
  dplyr::filter(group == "Sepsis", time %in% c("00h", "72h")) %>%
  dplyr::select(subject, time, dplyr::all_of(score_lookup$score_column)) %>%
  tidyr::pivot_wider(
    names_from = time, values_from = dplyr::all_of(score_lookup$score_column),
    names_glue = "{.value}___{time}"
  ) %>%
  dplyr::inner_join(sepsis_at_enrollment, by = "subject")

build_landmark_data <- function(one_method) {
  prefix <- methods$prefix[match(one_method, methods$method)]
  dd <- paired_input %>% dplyr::mutate(
    death28 = as.integer(mortality28 == "Died")
  )
  for (d in domains) {
    col_00 <- paste0(prefix, d, "___00h")
    col_72 <- paste0(prefix, d, "___72h")
    dd[[paste0(d, "_baseline")]] <- as.numeric(dd[[col_00]])
    dd[[paste0(d, "_delta")]] <- as.numeric(dd[[col_72]]) - as.numeric(dd[[col_00]])
  }
  required <- c("death28", "age", "sex",
                paste0(domains, "_baseline"), paste0(domains, "_delta"))
  dd <- dd[stats::complete.cases(dd[, required, drop = FALSE]), , drop = FALSE]
  dd <- droplevels(dd)
  if (nrow(dd) < 20L || length(unique(dd$death28)) < 2L) {
    stop("Insufficient complete paired 00h/72h samples for ", one_method)
  }
  for (d in domains) {
    delta <- dd[[paste0(d, "_delta")]]
    sd_delta <- stats::sd(delta)
    if (!is.finite(sd_delta) || sd_delta <= 0)
      stop("Undefined SD of 0-72h change for ", d, " / ", one_method)
    dd[[paste0(d, "_delta_z")]] <- (delta - mean(delta)) / sd_delta
  }
  dd$method <- one_method
  dd
}

landmark_data_list <- lapply(methods$method, build_landmark_data)
names(landmark_data_list) <- methods$method
landmark_patient_data <- dplyr::bind_rows(landmark_data_list)
write_output(landmark_patient_data, "Ghana_72h_landmark_paired_patient_level_changes")
landmark_sample_counts <- landmark_patient_data %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(n = dplyr::n(), deaths = sum(death28 == 1),
                   survivors = sum(death28 == 0), .groups = "drop")
write_output(landmark_sample_counts, "Ghana_72h_landmark_sample_sizes")

# Observed deltas by 28-day outcome. Wilcoxon is descriptive and unadjusted.
landmark_long <- landmark_patient_data %>%
  dplyr::select(subject, method, mortality28, death28,
                dplyr::all_of(paste0(domains, "_delta"))) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(paste0(domains, "_delta")),
    names_to = "domain", values_to = "delta"
  ) %>%
  dplyr::mutate(
    domain = sub("_delta$", "", domain),
    domain = factor(domain, levels = domains),
    mortality28 = factor(mortality28, levels = c("Survived", "Died"))
  )

landmark_descriptive <- purrr::map_dfr(methods$method, function(m) {
  purrr::map_dfr(domains, function(d) {
    dd <- landmark_long %>%
      dplyr::filter(method == m, as.character(domain) == d)
    a <- dd$delta[dd$mortality28 == "Survived"]
    b <- dd$delta[dd$mortality28 == "Died"]
    wt <- stats::wilcox.test(b, a, exact = FALSE)
    tibble::tibble(
      method = m, domain = d, n = nrow(dd),
      n_survived = length(a), n_died = length(b),
      mean_survived = mean(a), mean_died = mean(b),
      median_survived = stats::median(a), median_died = stats::median(b),
      q1_survived = unname(stats::quantile(a, 0.25)),
      q3_survived = unname(stats::quantile(a, 0.75)),
      q1_died = unname(stats::quantile(b, 0.25)),
      q3_died = unname(stats::quantile(b, 0.75)),
      Wilcoxon_W = unname(wt$statistic), p_wilcoxon = wt$p.value
    )
  })
}) %>%
  dplyr::group_by(method) %>%
  dplyr::mutate(q_BH_3domains = stats::p.adjust(p_wilcoxon, method = "BH")) %>%
  dplyr::ungroup()
write_output(landmark_descriptive, "Ghana_72h_landmark_DRT_changes_by_28d_outcome_Wilcoxon")

# One joint age/sex-adjusted model for all THREE absolute changes AND all THREE
# enrollment scores. All six molecular variables are in the same regression.
landmark_terms <- paste0(domains, "_delta_z")
landmark_baseline_terms <- paste0(domains, "_baseline")
landmark_formula <- stats::reformulate(
  c(landmark_terms, landmark_baseline_terms, "age", "sex"),
  response = "death28"
)

extract_landmark_or <- function(fit, dd, one_method, model_name) {
  beta <- stats::coef(fit)[landmark_terms]
  se <- sqrt(diag(stats::vcov(fit)))[landmark_terms]
  if (any(!is.finite(beta)) || any(!is.finite(se)) || any(se <= 0)) {
    stop("Nonestimable landmark regression coefficients for ", one_method,
         " (", model_name, ")")
  }
  tibble::tibble(
    method = one_method, model = model_name, domain = domains,
    n = nrow(dd), deaths = sum(dd$death28 == 1),
    change_SD = vapply(domains, function(d) stats::sd(dd[[paste0(d, "_delta")]]), numeric(1)),
    log_OR = unname(beta), SE = unname(se), OR = exp(unname(beta)),
    lower95 = exp(unname(beta - 1.96 * se)),
    upper95 = exp(unname(beta + 1.96 * se)),
    p = unname(2 * stats::pnorm(abs(beta / se), lower.tail = FALSE))
  ) %>%
    dplyr::mutate(q_BH_3domains = stats::p.adjust(p, method = "BH"))
}

landmark_logistic <- purrr::map_dfr(methods$method, function(m) {
  dd <- landmark_data_list[[m]]
  fit <- stats::glm(landmark_formula, data = dd, family = stats::binomial())
  if (!fit$converged) warning("Ordinary landmark logistic regression did not converge for ", m)
  extract_landmark_or(fit, dd, m, "Joint logistic")
})
write_output(landmark_logistic, "Ghana_72h_landmark_joint_delta_DRT_mortality_logistic")

# Collinearity diagnostics for the three delta coefficients conditional on all
# other regressors. Particularly relevant in this small landmark sample.
landmark_delta_vifs <- purrr::map_dfr(methods$method, function(m) {
  dd <- landmark_data_list[[m]]
  all_x <- c(landmark_terms, landmark_baseline_terms, "age", "sex")
  purrr::map_dfr(landmark_terms, function(term) {
    auxiliary <- stats::lm(
      stats::reformulate(setdiff(all_x, term), response = term), data = dd
    )
    tibble::tibble(
      method = m, domain = sub("_delta_z$", "", term),
      n = stats::nobs(auxiliary),
      VIF = 1 / (1 - summary(auxiliary)$r.squared)
    )
  })
})
write_output(landmark_delta_vifs, "Ghana_72h_landmark_joint_delta_domain_VIF")

# Firth: report profile penalized-likelihood CI and p-values. If not installed,
# warn and preserve the logistic/other sections, rather than silently switching.
landmark_firth <- tibble::tibble()
if (RUN_FIRTH_SENSITIVITY) {
  landmark_firth <- purrr::map_dfr(methods$method, function(m) {
    dd <- landmark_data_list[[m]]
    fit <- tryCatch(
      logistf::logistf(landmark_formula, data = dd, firth = TRUE, pl = TRUE),
      error = function(e) {
        warning("Firth landmark analysis failed for ", m, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(fit)) {
      return(tibble::tibble(
        method = m, domain = domains, n = nrow(dd),
        deaths = sum(dd$death28 == 1), OR = NA_real_,
        lower95 = NA_real_, upper95 = NA_real_, p = NA_real_,
        note = "Firth model failed"
      ))
    }
    tibble::tibble(
      method = m, domain = domains, n = nrow(dd),
      deaths = sum(dd$death28 == 1),
      OR = exp(unname(fit$coefficients[landmark_terms])),
      lower95 = exp(unname(fit$ci.lower[landmark_terms])),
      upper95 = exp(unname(fit$ci.upper[landmark_terms])),
      p = unname(fit$prob[landmark_terms]),
      note = "Firth; profile penalized-likelihood intervals"
    )
  }) %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(q_BH_3domains = stats::p.adjust(p, method = "BH")) %>%
    dplyr::ungroup()
  write_output(landmark_firth, "Ghana_72h_landmark_joint_delta_DRT_mortality_Firth_sensitivity")
} else {
  warning("Optional 'logistf' package not installed. Skipping Firth landmark sensitivity. ",
          "Install with install.packages('logistf').")
}

# A. Observed 0-72-hour changes by survival: HORIZONTAL boxplots with
#    individual participant points, shared change-score x-axis, zero line,
#    and the Wilcoxon p-value in each domain's facet heading. Transcriptomic projection only.
#    Domain panels are stacked Resistance -> Tolerance -> Damage, avoiding
#    cramped side-by-side horizontal plots. All data and tests are unchanged.
plot_landmark_raw <- landmark_long %>%
  dplyr::filter(method == "Transcriptomic projection") %>%
  dplyr::mutate(domain = factor(domain, levels = domains))

landmark_panel_labels <- landmark_descriptive %>%
  dplyr::filter(method == "Transcriptomic projection") %>%
  dplyr::mutate(
    domain = factor(domain, levels = domains),
    panel_label = paste0(as.character(domain), "   |   Wilcoxon ", format_p(p_wilcoxon))
  ) %>%
  dplyr::arrange(domain) %>%
  dplyr::select(domain, panel_label)

plot_landmark_raw <- plot_landmark_raw %>%
  dplyr::left_join(landmark_panel_labels, by = "domain") %>%
  dplyr::mutate(panel_label = factor(panel_label,
                                    levels = landmark_panel_labels$panel_label))

p_landmark_changes <- ggplot2::ggplot(
  plot_landmark_raw,
  ggplot2::aes(x = delta, y = mortality28,
               fill = mortality28, color = mortality28)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  ggplot2::geom_boxplot(orientation = "y", width = 0.55,
                        outlier.shape = NA, alpha = 0.48,
                        color = "black", linewidth = 0.65) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0, height = 0.12, seed = 2026),
    size = 1.8, alpha = 0.7
  ) +
  ggplot2::facet_wrap(~panel_label, ncol = 1, scales = "fixed") +
  # For a discrete y-axis the last level is displayed at the top.
  ggplot2::scale_y_discrete(limits = c("Died", "Survived")) +
  ggplot2::scale_fill_manual(values = c(Survived = "#2A6EBB", Died = "#D65F44")) +
  ggplot2::scale_color_manual(values = c(Survived = "#2A6EBB", Died = "#D65F44")) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.06))) +
  ggplot2::labs(
    title = "A. Transcriptomic projection: 0–72-hour changes by 28-day outcome",
    x = "Change in domain score (72h − enrollment; SD units)",
    y = NULL) +
  ggprism::theme_prism(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(face = "bold", size = 12),
    strip.placement = "outside",
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    axis.text.y = ggplot2::element_text(face = "bold"),
    axis.title.x = ggplot2::element_text(face = "bold"),
    panel.spacing.y = grid::unit(0.65, "lines"),
    plot.caption = ggplot2::element_text(size = 9, hjust = 0)
  )
save_plot(p_landmark_changes, "Ghana_RNA_LASSO_72h_landmark_changes_boxplots")

# B. Mutually adjusted joint logistic odds ratios per 1 SD higher DELTA.
# All three baseline scores and age/sex are included in the model.
landmark_forest_data <- landmark_logistic %>%
  dplyr::filter(method == "Transcriptomic projection") %>%
  dplyr::mutate(
    domain = factor(domain, levels = rev(domains)),
    annotation = paste0(
      sprintf("%.2f (%.2f–%.2f)", OR, lower95, upper95),
      "\n", format_p(p)
    )
  )

p_landmark_forest <- ggplot2::ggplot(
  landmark_forest_data,
  ggplot2::aes(x = OR, y = domain, color = domain)
) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = lower95, xmax = upper95),
    orientation = "y", width = 0.16, linewidth = 0.85
  ) +
  ggplot2::geom_point(size = 3.8) +
  ggplot2::geom_text(
    ggplot2::aes(x = upper95 * 1.1, label = annotation),
    hjust = 0, size = 3.4, fontface = "bold", show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = c(
    Resistance = "#2A6EBB", Tolerance = "#2CA25F", Damage = "#E6550D"
  )) +
  ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.08, 0.85))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "B. Adjusted 28-day mortality associations",
    x = "Adjusted OR per 1-SD increase in 0–72-hour change",
    y = NULL) +
  ggprism::theme_prism(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    axis.text.y = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(10, 160, 10, 10)
  )
save_plot(p_landmark_forest, "Ghana_RNA_LASSO_72h_landmark_joint_mortality_forest")

# Publication figure combining A and B side by side. Uses grid (base R),
# avoiding any additional graphics package requirements.

save_landmark_combined <- function(left_plot, right_plot, basename) {
  # Save the two-panel figure to one landscape PDF page (8.59 × 6.54 in).
  left_plot <- left_plot + ggplot2::labs(caption = NULL) +
    ggplot2::theme(plot.caption = ggplot2::element_blank())
  right_plot <- right_plot + ggplot2::labs(caption = NULL) +
    ggplot2::theme(plot.caption = ggplot2::element_blank())
  # Compact type for the mandated 8.59 x 6.54 inch two-panel layout.
  left_plot <- left_plot + ggplot2::theme(
    text = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_text(size = 10),
    strip.text = ggplot2::element_text(size = 8))
  right_plot <- right_plot + ggplot2::theme(
    text = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_text(size = 10),
    strip.text = ggplot2::element_text(size = 8))
  draw_combined <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      layout = grid::grid.layout(1, 2, widths = grid::unit(c(1.25, 1), "null"))))
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::grid.draw(ggplot2::ggplotGrob(left_plot))
    grid::upViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grid::grid.draw(ggplot2::ggplotGrob(right_plot))
    grid::upViewport(2)
  }
  grDevices::pdf(file.path(out_dir, paste0(basename, ".pdf")),
                 width = 8.59, height = 6.54, bg = "white", useDingbats = FALSE)
  tryCatch(draw_combined(), finally = grDevices::dev.off())
}

save_landmark_combined(
  p_landmark_changes, p_landmark_forest,
  "Ghana_RNA_LASSO_72h_landmark_two_panel_figure"
)

cat("\n=== 72-hour landmark paired sample sizes ===\n")
print(landmark_sample_counts)
cat("\n=== 72-hour landmark descriptive changes (Wilcoxon) ===\n")
print(landmark_descriptive)
cat("\n=== 72-hour landmark joint mortality logistic regression ===\n")
print(landmark_logistic)
cat("\n=== 72-hour landmark delta VIF diagnostics ===\n")
print(landmark_delta_vifs)
if (nrow(landmark_firth)) {
  cat("\n=== 72-hour landmark Firth sensitivity ===\n")
  print(landmark_firth)
}

cat("\n=== Baseline sepsis versus healthy (Wilcoxon/descriptive) ===\n")
print(baseline_table)
cat("\n=== Baseline MEWS: joint domain (age/sex adjusted) ===\n")
print(mews_table)
cat("\n=== Baseline joint logistic regression for 28-day mortality (age/sex adjusted) ===\n")
print(baseline_mortality_table)
cat("\n=== Longitudinal 28-day survival: global time x mortality tests ===\n")
print(trajectory_tests)
cat("\n=== Time-specific adjusted differences ===\n")
print(trajectory_time_contrasts)
cat("\nAnalysis files saved in: ", out_dir, "\n", sep = "")
