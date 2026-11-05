# =============================================================================
# RESERVE-U: SUPPORTIVE pooled joint domain analysis
# Script: 05_RESERVE_U_Pooled_DRT_Clinical.R
# =============================================================================
#
# ANALYSES
#   Cohort-specific primary analyses live in the two separate scripts.
#   Harmonized pooled 30-day mortality logistic regression AND ordinal MEWS.
#   Original multi-omic PC1 scores; within-cohort standardization.
#   All 3 domain scores in each joint model + 6 clinical covariates + cohort.
#   No separate-domain mortality or MEWS clinical-outcome models.
#   Combined cohort/pooled infection association forests (HIV, malaria, MTB).
#   POOLED 30-day mortality effect modification by HIV, malaria, MTB:
#     one domain x infection interaction per joint model, cohort adjusted only.
#   No pooled descriptive boxplots or pathogen-density plots in the main analyses.
#   Exploratory urine LAM grade (0-4) versus multi-omic D/R/T scores:
#     cohort-specific and pooled HC3-robust regressions (nominal p-values only),
#     one adjusted forest figure.
#
# INPUT
#   Paths for the 4 source CSVs are set below (relative to the project root).
#   Accepts updated sample_id / multiomic_PC1__ score-file columns and older aliases.
#
# OUTPUT
#   DRT_scores_output/Clinical analyses/Pooled
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggprism)
  library(tibble)
  library(purrr)
  library(MASS)
})

# ------------------------------------------------------------
# INPUT / OUTPUT PATHS (relative to the project root; edit if needed).
# Scores are produced by publication_code/01_DRT_scores.R.
# ------------------------------------------------------------
tor_clinical_path <- "data/tororo_all_298_clinical_data.csv"
tor_scores_path   <- "DRT_scores_output/TOR/TOR_DRT_scores_wide_perSD.csv"
ebb_clinical_path <- "data/entebbe_all_301_clinical_data.csv"
ebb_scores_path   <- "DRT_scores_output/EBB/EBB_DRT_scores_wide_perSD.csv"

out_dir <- "DRT_scores_output/Clinical analyses/Pooled"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(out_dir)) stop("Could not create output directory: ", out_dir)

read_csv <- function(path) utils::read.csv(path, check.names = FALSE,
                                           stringsAsFactors = FALSE)

clean_binary <- function(x) {
  out <- trimws(as.character(x))
  out[!out %in% c("0", "1")] <- NA_character_
  factor(out, levels = c("0", "1"))
}

prepare_one_cohort <- function(clinical_file, score_file, cohort_name) {
  clinical <- read_csv(clinical_file)
  scores <- read_csv(score_file)
  # Some RESERVE-U-2-TOR clinical exports use 'gender' rather than 'sex'.
  if (!"sex" %in% names(clinical) && "gender" %in% names(clinical)) {
    names(clinical)[names(clinical) == "gender"] <- "sex"
  }
  # Harmonize new domain export schema with the original pooled-model names.
  # Clinical IDs remain 'pid'; scores now use 'sample_id' (older exports: 'pid').
  if (!"pid" %in% names(clinical))
    stop(cohort_name, ": clinical CSV must include pid")
  if (!"pid" %in% names(scores)) {
    if (!"sample_id" %in% names(scores))
      stop(cohort_name, ": scores must include sample_id or pid")
    scores$pid <- scores$sample_id
  } else if ("sample_id" %in% names(scores)) {
    if (!identical(trimws(as.character(scores$pid)),
                   trimws(as.character(scores$sample_id))))
      stop(cohort_name, ": pid and sample_id disagree in score file")
  }
  clinical$pid <- trimws(as.character(clinical$pid))
  scores$pid <- trimws(as.character(scores$pid))
  if (anyNA(clinical$pid) || anyNA(scores$pid) ||
      any(!nzchar(clinical$pid)) || any(!nzchar(scores$pid)))
    stop(cohort_name, ": missing or empty participant IDs")
  if (anyDuplicated(clinical$pid) || anyDuplicated(scores$pid))
    stop(cohort_name, ": duplicate participant IDs")

  # The new multi-omic columns use 'multiomic_PC1__'. Retain the old internal
  # names so all downstream joint models, standardized scores and plots agree.
  domains <- c("Resistance", "Tolerance", "Damage")
  required_scores <- paste0("multiomic_option2_PC1__", domains)
  updated_scores <- paste0("multiomic_PC1__", domains)
  for (i in seq_along(domains)) {
    old_col <- required_scores[i]
    new_col <- updated_scores[i]
    if (!(old_col %in% names(scores)) && !(new_col %in% names(scores)))
      stop(cohort_name, ": missing PC1 score for ", domains[i],
           " (expected ", new_col, " or ", old_col, ")")
    if (new_col %in% names(scores) && old_col %in% names(scores)) {
      x <- suppressWarnings(as.numeric(scores[[new_col]]))
      y <- suppressWarnings(as.numeric(scores[[old_col]]))
      same_values <- isTRUE(all.equal(x, y, tolerance = 1e-8,
                                      check.attributes = FALSE))
      if (!same_values)
        stop(cohort_name, ": conflicting old/new PC1 scores for ", domains[i])
    }
    scores[[old_col]] <- suppressWarnings(as.numeric(scores[[
      if (new_col %in% names(scores)) new_col else old_col
    ]]))
  }
  if ("cohort" %in% names(scores)) {
    cohort_code <- if (cohort_name == "RESERVE-U-2-TOR") "TOR" else "EBB"
    scores <- scores[toupper(trimws(as.character(scores$cohort))) == cohort_code,
                     , drop = FALSE]
    if (nrow(scores) == 0L)
      stop(cohort_name, ": no molecular rows with cohort = ", cohort_code)
  }
  matched_n <- sum(scores$pid %in% clinical$pid)
  if (matched_n == 0L)
    stop(cohort_name, ": no matching participant IDs between files")
  cat(cohort_name, ": score rows = ", nrow(scores),
      "; matched clinical IDs = ", matched_n, "\n", sep = "")
  if (matched_n < nrow(scores))
    warning(cohort_name, ": ", nrow(scores) - matched_n,
            " molecular IDs have no corresponding clinical record")

  # If a premerged clinical export was selected, avoid .x/.y suffixes:
  # the separately selected score file remains authoritative.
  clinical <- clinical[, setdiff(names(clinical), c(required_scores, updated_scores)),
                       drop = FALSE]
  dat <- dplyr::inner_join(clinical,
                           scores[, c("pid", required_scores), drop = FALSE],
                           by = "pid")
  if (cohort_name == "RESERVE-U-2-TOR") {
    needed <- c("age", "sex", "illnessduration_to_enroll", "death30d",
                "mews_score", "hiv_rdt_result", "malaria_rdt_result", "microtbdx")
    if (!all(needed %in% names(dat))) stop("Missing RESERVE-U-2-TOR clinical columns: ",
                                            paste(setdiff(needed, names(dat)), collapse = ", "))
    sex <- dplyr::case_when(as.character(dat$sex) %in% c("1", "Male", "M") ~ "Male",
                            as.character(dat$sex) %in% c("2", "Female", "F") ~ "Female",
                            TRUE ~ NA_character_)
    hiv <- dat$hiv_rdt_result
    malaria <- dat$malaria_rdt_result
    duration <- dat$illnessduration_to_enroll
    mews <- dat$mews_score
    # Urine LAM: grade 0 is a documented negative test, NEVER an untested person.
    lam_fields <- c("urine_lam_done", "urine_lam_result", "urine_lam_band_grade")
    if (!all(lam_fields %in% names(dat)))
      stop(cohort_name, ": missing urine LAM columns: ",
           paste(setdiff(lam_fields, names(dat)), collapse = ", "))
    lam_done <- dat$urine_lam_done
    lam_result <- dat$urine_lam_result
    lam_band <- dat$urine_lam_band_grade
  } else {
    needed <- c("age", "sex", "illnessduration_enroll", "death30d",
                "MEWS_score", "hivrdtresult", "malariardtresult", "microtbdx")
    if (!all(needed %in% names(dat))) stop("Missing RESERVE-U-1-EBB clinical columns: ",
                                            paste(setdiff(needed, names(dat)), collapse = ", "))
    sex <- dplyr::case_when(as.character(dat$sex) %in% c("0", "Male", "M") ~ "Male",
                            as.character(dat$sex) %in% c("1", "Female", "F") ~ "Female",
                            TRUE ~ NA_character_)
    hiv <- dat$hivrdtresult
    malaria <- dat$malariardtresult
    duration <- dat$illnessduration_enroll
    mews <- dat$MEWS_score
    lam_fields <- c("urinetblam", "urinetblamresult", "urinetblamgrade")
    if (!all(lam_fields %in% names(dat)))
      stop(cohort_name, ": missing urine LAM columns: ",
           paste(setdiff(lam_fields, names(dat)), collapse = ", "))
    lam_done <- dat$urinetblam
    lam_result <- dat$urinetblamresult
    lam_band <- dat$urinetblamgrade
  }
  lam_done <- suppressWarnings(as.numeric(as.character(lam_done)))
  lam_result <- suppressWarnings(as.numeric(as.character(lam_result)))
  lam_band <- suppressWarnings(as.numeric(as.character(lam_band)))
  # Only documented tests with valid positive/negative results are usable.
  if (any(!is.na(lam_result) & !(lam_result %in% c(0, 1))))
    stop(cohort_name, ": unexpected urine LAM result encoding (expected 0/1)")
  if (any(!is.na(lam_result) & (is.na(lam_done) | lam_done != 1)))
    stop(cohort_name, ": LAM result recorded without completed LAM test")
  if (any(lam_result == 1 & !is.na(lam_band) &
          !(lam_band %in% 1:4), na.rm = TRUE))
    stop(cohort_name, ": unexpected positive urine LAM grade (expected 1-4)")
  if (any(lam_result == 0 & !is.na(lam_band), na.rm = TRUE))
    stop(cohort_name, ": LAM-negative patient unexpectedly has a positive band grade")
  lam_positive <- ifelse(lam_done == 1 & lam_result %in% c(0, 1),
                         lam_result, NA_real_)
  lam_grade <- ifelse(lam_positive == 0, 0,
                      ifelse(lam_positive == 1 & lam_band %in% 1:4,
                             lam_band, NA_real_))
  numeric_scores <- lapply(required_scores, function(col) suppressWarnings(as.numeric(dat[[col]])))
  names(numeric_scores) <- c("Resistance", "Tolerance", "Damage")
  out <- tibble::tibble(
    pid = dat$pid, cohort = cohort_name,
    age = suppressWarnings(as.numeric(dat$age)),
    sex = factor(sex, levels = c("Male", "Female")),
    illness_duration = suppressWarnings(as.numeric(duration)),
    hiv = clean_binary(hiv), malaria = clean_binary(malaria),
    tb = clean_binary(dat$microtbdx),
    death30d = suppressWarnings(as.numeric(as.character(dat$death30d))),
    mews = suppressWarnings(as.numeric(as.character(mews))),
    lam_positive = lam_positive,
    lam_grade = lam_grade,
    lam_band_grade = ifelse(lam_positive == 1, lam_band, NA_real_),
    Resistance = numeric_scores$Resistance,
    Tolerance = numeric_scores$Tolerance,
    Damage = numeric_scores$Damage
  )
  out <- out %>%
    dplyr::filter(stats::complete.cases(dplyr::across(
      dplyr::all_of(c("Resistance", "Tolerance", "Damage")))))
  # Scale relative to the cohort-specific distribution of original scores,
  # BEFORE limiting to participants with complete clinical/outcome data.
  for (domain in c("Resistance", "Tolerance", "Damage")) {
    sx <- stats::sd(out[[domain]])
    if (!is.finite(sx) || sx <= 0) stop(cohort_name, ": invalid SD for ", domain)
    out[[paste0(domain, "_z")]] <- (out[[domain]] - mean(out[[domain]])) / sx
  }
  cat("\n", cohort_name, ": matched multi-omic participants = ", nrow(out), "\n", sep = "")
  out
}

tor <- prepare_one_cohort(tor_clinical_path, tor_scores_path, "RESERVE-U-2-TOR")
ebb <- prepare_one_cohort(ebb_clinical_path, ebb_scores_path, "RESERVE-U-1-EBB")
pooled <- dplyr::bind_rows(tor, ebb) %>%
  dplyr::mutate(cohort = factor(cohort, levels = c("RESERVE-U-2-TOR", "RESERVE-U-1-EBB")))

predictors <- c("Resistance_z", "Tolerance_z", "Damage_z", "age", "sex",
                "illness_duration", "hiv", "malaria", "tb")
complete_vars <- c("death30d", "cohort", predictors)
pooled_cc <- pooled[stats::complete.cases(pooled[, complete_vars]), , drop = FALSE]
if (any(!pooled_cc$death30d %in% c(0, 1))) stop("Invalid 30-day mortality encoding")
if (length(unique(pooled_cc$death30d)) != 2L) stop("Pooled mortality has only one outcome group")

utils::write.csv(pooled_cc, file.path(out_dir, "RESERVE_U_POOLED_DRT_complete_case_joint_30d_data.csv"), row.names = FALSE)

# MEWS uses its own complete-case population and the same molecular/clinical
# covariates as mortality, including the cohort fixed effect.
mews_complete_vars <- c("mews", "cohort", predictors)
pooled_mews_cc <- pooled[stats::complete.cases(pooled[, mews_complete_vars]), , drop = FALSE]
if (any(!is.finite(pooled_mews_cc$mews)) || length(unique(pooled_mews_cc$mews)) < 3L)
  stop("MEWS must have at least 3 observed ordered categories and finite values")
pooled_mews_cc <- droplevels(pooled_mews_cc)
pooled_mews_cc$mews_ordinal <- ordered(pooled_mews_cc$mews,
                                       levels = sort(unique(pooled_mews_cc$mews)))
utils::write.csv(pooled_mews_cc,
  file.path(out_dir, "RESERVE_U_POOLED_DRT_complete_case_joint_MEWS_data.csv"),
  row.names = FALSE)

# Supportive pooled model with cohort fixed effect; no cohort-by-domain interactions.
pooled_formula <- stats::reformulate(c(predictors, "cohort"), response = "death30d")
pooled_fit <- stats::glm(pooled_formula, family = stats::binomial(), data = pooled_cc)

# Same specification separately by cohort, for visual comparison only.
cohort_formula <- stats::reformulate(predictors, response = "death30d")
fit_tor <- stats::glm(cohort_formula, family = stats::binomial(),
                      data = droplevels(pooled_cc[pooled_cc$cohort == "RESERVE-U-2-TOR", , drop = FALSE]))
fit_ebb <- stats::glm(cohort_formula, family = stats::binomial(),
                      data = droplevels(pooled_cc[pooled_cc$cohort == "RESERVE-U-1-EBB", , drop = FALSE]))

extract_or <- function(fit, cohort_label, n_deaths) {
  terms <- paste0(c("Resistance", "Tolerance", "Damage"), "_z")
  b <- stats::coef(fit)[terms]
  se <- sqrt(diag(stats::vcov(fit)))[terms]
  tibble::tibble(
    source = cohort_label,
    domain = c("Resistance", "Tolerance", "Damage"),
    n = stats::nobs(fit), deaths = n_deaths,
    OR = exp(unname(b)), lower95 = exp(unname(b - 1.96 * se)),
    upper95 = exp(unname(b + 1.96 * se)),
    p = 2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
  ) %>%
    dplyr::mutate(q_BH = stats::p.adjust(p, method = "BH"),
                  estimate_95CI = sprintf("%.2f (%.2f–%.2f)", OR, lower95, upper95))
}

pooled_estimates <- dplyr::bind_rows(
  extract_or(pooled_fit, "Pooled (cohort-adjusted)", sum(pooled_cc$death30d == 1)),
  extract_or(fit_tor, "RESERVE-U-2-TOR", sum(pooled_cc$death30d[pooled_cc$cohort == "RESERVE-U-2-TOR"] == 1)),
  extract_or(fit_ebb, "RESERVE-U-1-EBB", sum(pooled_cc$death30d[pooled_cc$cohort == "RESERVE-U-1-EBB"] == 1))
)
utils::write.csv(pooled_estimates,
  file.path(out_dir, "RESERVE_U_POOLED_supportive_joint_DRT_30d_estimates.csv"), row.names = FALSE)

# ============================================================
# SUPPORTIVE pooled joint MEWS proportional-odds model.
# Three domains enter SIMULTANEOUSLY; the same six clinical covariates and
# cohort fixed effect are used as in the 30-day mortality analysis.
# ============================================================
mews_formula <- stats::reformulate(c(predictors, "cohort"), response = "mews_ordinal")
mews_pooled_fit <- MASS::polr(mews_formula, data = pooled_mews_cc,
                              Hess = TRUE, method = "logistic")
mews_cohort_formula <- stats::reformulate(predictors, response = "mews_ordinal")
mews_tor_data <- droplevels(pooled_mews_cc[pooled_mews_cc$cohort == "RESERVE-U-2-TOR", , drop=FALSE])
mews_ebb_data <- droplevels(pooled_mews_cc[pooled_mews_cc$cohort == "RESERVE-U-1-EBB", , drop=FALSE])
# Drop unused ordered-outcome levels within each cohort.
mews_tor_data$mews_ordinal <- ordered(mews_tor_data$mews,
                                      levels = sort(unique(mews_tor_data$mews)))
mews_ebb_data$mews_ordinal <- ordered(mews_ebb_data$mews,
                                      levels = sort(unique(mews_ebb_data$mews)))
if (nlevels(mews_tor_data$mews_ordinal) < 3L ||
    nlevels(mews_ebb_data$mews_ordinal) < 3L)
  stop("Not enough distinct MEWS levels in one of the two cohorts")
mews_tor_fit <- MASS::polr(mews_cohort_formula, data = mews_tor_data,
                           Hess = TRUE, method = "logistic")
mews_ebb_fit <- MASS::polr(mews_cohort_formula, data = mews_ebb_data,
                           Hess = TRUE, method = "logistic")
# Extract coefficient-only odds ratios; polr vcov also includes thresholds.
mews_estimates <- dplyr::bind_rows(
  extract_or(mews_pooled_fit, "Pooled (cohort-adjusted)", NA_integer_),
  extract_or(mews_tor_fit, "RESERVE-U-2-TOR", NA_integer_),
  extract_or(mews_ebb_fit, "RESERVE-U-1-EBB", NA_integer_)
)
utils::write.csv(mews_estimates,
  file.path(out_dir, "RESERVE_U_POOLED_supportive_joint_DRT_MEWS_estimates.csv"),
  row.names = FALSE)

mews_counts <- pooled_mews_cc %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarise(n = dplyr::n(), n_mews_levels = dplyr::n_distinct(mews), .groups="drop")
utils::write.csv(mews_counts, file.path(out_dir, "RESERVE_U_POOLED_joint_DRT_MEWS_sample_sizes.csv"),
  row.names = FALSE)

# Sample sizes, deaths, and events-per-parameter transparency.
pooled_counts <- pooled_cc %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarise(n = dplyr::n(), deaths = sum(death30d == 1), .groups = "drop")
pooled_counts <- dplyr::bind_rows(pooled_counts,
 tibble::tibble(cohort = factor("Pooled", levels = c(levels(pooled_cc$cohort), "Pooled")),
                n = nrow(pooled_cc), deaths = sum(pooled_cc$death30d == 1)))
utils::write.csv(pooled_counts,
  file.path(out_dir, "RESERVE_U_POOLED_joint_DRT_sample_sizes.csv"), row.names = FALSE)

# Domain VIFs in the FULL pooled model (including cohort and all other covariates).
pooled_domain_vif <- purrr::map_dfr(c("Resistance_z", "Tolerance_z", "Damage_z"), function(term) {
  fit <- stats::lm(stats::reformulate(setdiff(c(predictors, "cohort"), term), response = term), data = pooled_cc)
  tibble::tibble(domain = sub("_z$", "", term), VIF = 1/(1-summary(fit)$r.squared), n = stats::nobs(fit))
})
utils::write.csv(pooled_domain_vif,
  file.path(out_dir, "RESERVE_U_POOLED_joint_DRT_domain_VIFs.csv"), row.names = FALSE)

mews_domain_vif <- purrr::map_dfr(c("Resistance_z", "Tolerance_z", "Damage_z"), function(term) {
  fit <- stats::lm(stats::reformulate(setdiff(c(predictors, "cohort"), term),
                                      response = term), data = pooled_mews_cc)
  tibble::tibble(domain = sub("_z$", "", term),
                 VIF = 1/(1-summary(fit)$r.squared), n = stats::nobs(fit))
})
utils::write.csv(mews_domain_vif,
  file.path(out_dir, "RESERVE_U_POOLED_joint_DRT_MEWS_domain_VIFs.csv"), row.names = FALSE)

# ============================================================
# JOURNAL-STYLE FOREST FIGURES: shared display function
# ------------------------------------------------------------
# Three domain groups; within each group: Tororo, Entebbe, pooled.
# The RIGHT SIDE is a fixed, aligned numerical table, not labels positioned
# relative to each confidence interval. This avoids overlap/clipping.
# ORs are plotted on a logarithmic axis; adjusted beta values use a linear
# axis. All figures are saved as vector PDF at 8.59 x 6.54 inches.
# ============================================================
domain_colors <- c(Resistance = "#2A6EBB", Tolerance = "#2CA25F", Damage = "#E6550D")
forest_sources <- c("RESERVE-U-2-TOR", "RESERVE-U-1-EBB", "Pooled (cohort-adjusted)")
forest_domains <- c("Resistance", "Tolerance", "Damage")

journal_forest <- function(dat, figure_title, axis_title,
                           measure = c("OR", "beta")) {
  measure <- match.arg(measure)
  needed <- c(if (measure == "OR") "OR" else "beta",
              "lower95", "upper95", "p", "domain", "source")
  if (!all(needed %in% names(dat))) {
    stop("Journal forest is missing: ", paste(setdiff(needed, names(dat)), collapse = ", "))
  }
  if (anyDuplicated(dat[, c("domain", "source")])) {
    stop("Journal forest requires one row per domain and source")
  }

  # One common, fixed layout for all five figures. Everything is drawn INSIDE
  # the plotting area: no clipped labels, no oversized right margins, and no
  # dependence on the width of an OR/CI. 'sans' matches ggprism throughout.
  # The headers have their own full-height rows rather than sharing the
  # Resistance/Tolerance/Damage estimate rows.
  cohort_x <- 0.025
  forest_left <- 0.255
  forest_right <- 0.607
  estimate_x <- 0.663
  pvalue_x <- 0.929
  group_names <- forest_domains
  group_header_y <- c(10.95, 7.53, 4.11)
  cohort_y <- list(c(10.05, 9.35, 8.65),
                   c(6.63, 5.93, 5.23),
                   c(3.21, 2.51, 1.81))

  dd <- dat %>%
    dplyr::mutate(
      domain = as.character(domain),
      source = as.character(source),
      estimate = if (measure == "OR") as.numeric(OR) else as.numeric(beta),
      group_number = match(domain, group_names),
      source_number = match(source, forest_sources)
    )
  if (any(is.na(dd$group_number)) || any(is.na(dd$source_number))) {
    stop("Unexpected domain/source label in journal forest")
  }
  if (nrow(dd) != 9L) warning("Expected 9 forest rows; found ", nrow(dd))
  dd$y <- mapply(function(g, c) cohort_y[[g]][c],
                  dd$group_number, dd$source_number)
  dd$source_label <- ifelse(dd$source == "Pooled (cohort-adjusted)",
                            "Pooled", dd$source)
  dd$source_face <- ifelse(dd$source == "Pooled (cohort-adjusted)",
                           "bold", "plain")
  dd$effect_label <- ifelse(
    is.finite(dd$estimate) & is.finite(dd$lower95) & is.finite(dd$upper95),
    sprintf("%.2f (%.2f, %.2f)", dd$estimate, dd$lower95, dd$upper95),
    "Not estimable"
  )
  dd$p_label <- ifelse(is.finite(dd$p),
                       ifelse(dd$p < 0.0001, "<0.0001", sprintf("%.3g", dd$p)),
                       "—")

  # Preserve actual log-OR vs. linear-beta scale while mapping its visual
  # positions into the forest-only column of the common table layout.
  observed <- c(dd$lower95, dd$upper95, dd$estimate)
  observed <- observed[is.finite(observed)]
  if (!length(observed)) stop("No finite estimates for ", figure_title)
  null <- if (measure == "OR") 1 else 0
  if (measure == "OR") {
    if (any(observed <= 0)) stop("OR and OR confidence limits must be positive")
    transform_x <- function(x) log(x)
    tick_candidates <- c(0.1, 0.2, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)
  } else {
    transform_x <- function(x) x
    tick_candidates <- pretty(range(c(observed, null)), n = 5)
  }
  tr <- range(transform_x(c(observed, null)))
  pad <- max(diff(tr) * 0.065, if (measure == "OR") 0.09 else 0.05)
  tr <- tr + c(-pad, pad)
  map_x <- function(x) {
    forest_left + (transform_x(x) - tr[1]) / diff(tr) *
      (forest_right - forest_left)
  }
  dd$x_mid <- map_x(dd$estimate)
  dd$x_lo <- map_x(dd$lower95)
  dd$x_hi <- map_x(dd$upper95)
  ticks <- tick_candidates[
    is.finite(tick_candidates) &
    transform_x(tick_candidates) >= tr[1] &
    transform_x(tick_candidates) <= tr[2]
  ]
  tick_df <- tibble::tibble(x = map_x(ticks),
                            label = format(ticks, trim = TRUE, scientific = FALSE))
  header_df <- tibble::tibble(domain = group_names,
                              y = group_header_y,
                              label = group_names)
  # The same sans-serif typeface and sizes are used in both numeric columns,
  # cohort names, tick labels, and headings; no monospace text or pseudo-table.
  text_face <- "sans"
  text_size <- 3.65

  p <- ggplot2::ggplot() +
    # Compact, restrained separators between the three domain groups.
    ggplot2::geom_segment(ggplot2::aes(x = cohort_x, xend = 0.985,
                                       y = 8.13, yend = 8.13),
                          color = "gray87", linewidth = 0.42) +
    ggplot2::geom_segment(ggplot2::aes(x = cohort_x, xend = 0.985,
                                       y = 4.71, yend = 4.71),
                          color = "gray87", linewidth = 0.42) +
    ggplot2::geom_segment(
      ggplot2::aes(x = map_x(null), xend = map_x(null),
                   y = 1.47, yend = 10.38),
      color = "gray48", linetype = "dashed", linewidth = 0.55) +
    # Forest confidence intervals, caps, and point estimates.
    ggplot2::geom_segment(
      data = dd,
      ggplot2::aes(x = x_lo, xend = x_hi, y = y, yend = y, color = domain),
      linewidth = 0.82, na.rm = TRUE) +
    ggplot2::geom_segment(
      data = dd,
      ggplot2::aes(x = x_lo, xend = x_lo, y = y - .115,
                   yend = y + .115, color = domain),
      linewidth = 0.82, na.rm = TRUE) +
    ggplot2::geom_segment(
      data = dd,
      ggplot2::aes(x = x_hi, xend = x_hi, y = y - .115,
                   yend = y + .115, color = domain),
      linewidth = 0.82, na.rm = TRUE) +
    ggplot2::geom_point(
      data = dd, ggplot2::aes(x = x_mid, y = y, color = domain),
      size = 3.0, na.rm = TRUE) +
    # Group labels occupy their own rows and align with the left text column.
    ggplot2::geom_text(
      data = header_df,
      ggplot2::aes(x = cohort_x, y = y, label = label, color = domain),
      hjust = 0, size = 4.35, fontface = "bold", family = text_face,
      show.legend = FALSE) +
    ggplot2::geom_text(
      data = dd,
      ggplot2::aes(x = cohort_x + 0.010, y = y,
                   label = source_label, fontface = source_face),
      hjust = 0, size = text_size, family = text_face, color = "black") +
    # Dedicated right-side numeric columns: fixed horizontal positions.
    ggplot2::geom_text(
      data = dd, ggplot2::aes(x = estimate_x, y = y, label = effect_label),
      hjust = 0, size = text_size, family = text_face, color = "black") +
    ggplot2::geom_text(
      data = dd, ggplot2::aes(x = pvalue_x, y = y, label = p_label),
      hjust = 0, size = text_size, family = text_face, color = "black") +
    ggplot2::annotate("text", x = cohort_x, y = 11.96, label = "Cohort",
                      hjust = 0, size = text_size, family = text_face,
                      fontface = "bold") +
    ggplot2::annotate("text", x = estimate_x, y = 11.96,
                      label = if (measure == "OR") "OR (95% CI)" else "Adjusted β (95% CI)",
                      hjust = 0, size = text_size, family = text_face,
                      fontface = "bold") +
    ggplot2::annotate("text", x = pvalue_x, y = 11.96, label = "p-value",
                      hjust = 0, size = text_size, family = text_face,
                      fontface = "bold") +
    ggplot2::geom_segment(
      ggplot2::aes(x = forest_left, xend = forest_right,
                   y = 1.13, yend = 1.13),
      color = "black", linewidth = 0.55) +
    ggplot2::geom_segment(
      data = tick_df,
      ggplot2::aes(x = x, xend = x, y = 1.13, yend = 1.03),
      color = "black", linewidth = 0.55) +
    ggplot2::geom_text(
      data = tick_df,
      ggplot2::aes(x = x, y = .73, label = label),
      size = text_size, family = text_face, color = "black") +
    ggplot2::annotate(
      "text", x = (forest_left + forest_right) / 2,
      y = 0.08, label = axis_title,
      size = 3.50, family = text_face, fontface = "bold") +
    ggplot2::scale_color_manual(values = domain_colors) +
    ggplot2::scale_x_continuous(limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_y_continuous(limits = c(-0.35, 12.45),
                                expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(title = figure_title) +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(
      text = ggplot2::element_text(family = text_face),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0,
                                         margin = ggplot2::margin(b = 12)),
      plot.caption = ggplot2::element_text(size = 10.5, color = "gray30",
                                           hjust = 0, margin = ggplot2::margin(t = 12)),
      plot.margin = ggplot2::margin(t = 10, r = 12, b = 10, l = 12)
    )
  p
}

save_journal_forest <- function(p, stem) {
  # Journal figures: vector PDF only, exactly 8.59 × 6.54 in (landscape).
  p <- p + ggplot2::labs(caption = NULL) +
    ggplot2::theme(plot.caption = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 13))
  tryCatch(print(p), error = function(e) {
    warning("Interactive plot preview failed for ", stem, ": ",
            conditionMessage(e), "; proceeding with PDF export")
  })
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(stem, ".pdf")), plot = p,
    device = "pdf", width = 8.59, height = 6.54, units = "in", bg = "white")
}


# 30-day mortality: three domains estimated in the SAME joint model.
p_pooled_joint <- journal_forest(
  pooled_estimates,
  figure_title = "Joint domain associations with 30-day mortality",
  axis_title = "Adjusted odds ratio for 30-day mortality per 1-SD increase",
  measure = "OR"
)
save_journal_forest(p_pooled_joint, "RESERVE_U_POOLED_supportive_joint_DRT_30d_forest")

cat("\nSupportive pooled joint results:\n")
print(pooled_estimates)
cat("\nVIFs:\n")
print(pooled_domain_vif)
cat("\nAll pooled files saved in:\n", out_dir, "\n")


# ============================================================
# Joint ordinal MEWS: same display and effect annotation as mortality.
# ============================================================
p_pooled_mews <- journal_forest(
  mews_estimates,
  figure_title = "Joint domain associations with MEWS",
  axis_title = "Adjusted odds ratio for MEWS per 1-SD increase",
  measure = "OR"
)
save_journal_forest(p_pooled_mews, "RESERVE_U_POOLED_supportive_joint_DRT_MEWS_forest")

# ============================================================
# ADJUSTED INFECTION-TYPE ASSOCIATIONS
# Each molecular domain is the dependent variable (per within-cohort SD).
# Exposure is HIV, malaria, or MTB (positive vs
# negative), adjusted for age, sex, illness duration, the OTHER major
# infections, and a cohort fixed effect when data are pooled.
# No other infection types are modeled.
# ============================================================
domain_names <- c("Resistance", "Tolerance", "Damage")
infection_vars <- c("hiv", "malaria", "tb")
infection_display <- c(hiv="HIV", malaria="Malaria", tb="MTB")

run_infection_association <- function(data, exposure, domain, source_label, add_cohort=FALSE) {
  outcome <- paste0(domain, "_z")
  adjustment <- c("age", "sex", "illness_duration", setdiff(infection_vars, exposure))
  if (add_cohort) adjustment <- c(adjustment, "cohort")
  required <- unique(c(outcome, exposure, adjustment))
  dd <- data[stats::complete.cases(data[, required, drop=FALSE]), , drop=FALSE]
  dd <- droplevels(dd)
  # Contrasts require both infection-status groups.
  if (nrow(dd) < 15L || length(unique(dd[[exposure]])) < 2L) {
    warning("Skipping insufficient ", source_label, ": ", exposure, " / ", domain)
    return(tibble::tibble(source=source_label, domain=domain, infection=infection_display[[exposure]],
      n=nrow(dd), beta=NA_real_, lower95=NA_real_, upper95=NA_real_, p=NA_real_))
  }
  form <- stats::reformulate(c(exposure, adjustment), response=outcome)
  fitted <- tryCatch(stats::lm(form, data=dd), error=function(e) {
    warning("Infection association fit failed: ", conditionMessage(e)); NULL
  })
  term <- paste0(exposure, "1")
  if (is.null(fitted) || !(term %in% names(stats::coef(fitted)))) {
    warning("Infection association coefficient missing for ", source_label, ": ", domain, " / ", exposure)
    return(tibble::tibble(source=source_label, domain=domain, infection=infection_display[[exposure]],
      n=nrow(dd), beta=NA_real_, lower95=NA_real_, upper95=NA_real_, p=NA_real_))
  }
  sm <- summary(fitted)$coefficients
  b <- unname(sm[term, "Estimate"])
  se <- unname(sm[term, "Std. Error"])
  p <- unname(sm[term, "Pr(>|t|)"])
  q <- stats::qt(.975, df=stats::df.residual(fitted))
  tibble::tibble(source=source_label, domain=domain,
    infection=infection_display[[exposure]], n=stats::nobs(fitted),
    beta=b, lower95=b-q*se, upper95=b+q*se, p=p)
}

infection_associations <- purrr::map_dfr(domain_names, function(domain) {
  purrr::map_dfr(infection_vars, function(exposure) {
    dplyr::bind_rows(
      run_infection_association(tor, exposure, domain, "RESERVE-U-2-TOR"),
      run_infection_association(ebb, exposure, domain, "RESERVE-U-1-EBB"),
      run_infection_association(pooled, exposure, domain, "Pooled (cohort-adjusted)", add_cohort=TRUE)
    )
  })
}) %>%
  dplyr::group_by(source, domain) %>%
  dplyr::mutate(q_BH = stats::p.adjust(p, method="BH")) %>%
  dplyr::ungroup()
utils::write.csv(infection_associations,
  file.path(out_dir, "RESERVE_U_POOLED_DRT_adjusted_infection_type_associations.csv"),
  row.names=FALSE)

# ============================================================
# INFECTION-TYPE COMPARISON FORESTS: HIV, malaria, MTB
# The forest/table layout matches both joint clinical-outcome plots.
# Adjusted β is correct here (the molecular domain score is the outcome),
# whereas the mortality/MEWS regressions report odds ratios.
# Each infection association adjusts for age, sex, illness duration, and
# the OTHER two infections; pooled estimates also adjust for cohort.
# ============================================================
infection_plot_order <- c("HIV", "Malaria", "MTB")
for (infection_name in infection_plot_order) {
  infection_df <- infection_associations %>%
    dplyr::filter(.data$infection == .env$infection_name)
  p_infection <- journal_forest(
    infection_df,
    figure_title = paste0("Adjusted ", infection_name, " associations with domain scores"),
    axis_title = "Adjusted difference in domain score (SD): positive vs negative",
    measure = "beta"
  )
  figure_stem <- paste0("RESERVE_U_POOLED_DRT_", infection_name,
                        "_cohort_and_pooled_forest")
  save_journal_forest(p_infection, figure_stem)
}


# ============================================================
# EXPLORATORY POOLED EFFECT MODIFICATION BY INFECTION TYPE
# 30-day mortality x (one domain x one infection per model)
# ------------------------------------------------------------
# For EACH exposure (HIV, malaria, microbiologically confirmed TB) and
# EACH molecular domain, fit ONE additional interaction coefficient:
#
# death30d ~ Resistance_z + Tolerance_z + Damage_z + cohort +
#            infection + [one selected domain]:infection
#
# There are 9 separate joint models (3 exposures x 3 candidate modifiers).
# No age, sex, illness-duration, or other infection adjustment in THIS block.
# These covariates are retained, unchanged, in the main outcome models above.
# Only the infection of interest enters a given effect-modification model.
# The infection main effect is ALWAYS included alongside its interaction.
#
# The null model for each 1-df likelihood-ratio test has the identical
# predictors and complete-case samples, but without the single interaction.
# Scale: within-cohort 1 SD of each multi-omic PC1 domain (already computed).
# Pooled analyses include a fixed effect for cohort; NO cohort interactions.
#
# IMPORTANT: Infection-specific subsets are based ONLY on mortality, ALL
# three domains, cohort and the infection being tested. Missing values in
# unrelated clinical variables do NOT exclude participants from this block.
# Results are exploratory; sparse exposed deaths and multiple testing limit
# precision. BH q-values are calculated across all NINE interaction tests.
# The interaction OR ratio compares the per-SD domain OR in infection-positive
# versus infection-negative patients; it is NOT an infection-risk OR.
# ============================================================

hte_domains <- c("Resistance", "Tolerance", "Damage")
hte_domain_terms <- paste0(hte_domains, "_z")
hte_exposures <- c(hiv = "HIV", malaria = "Malaria", tb = "MTB")

hte_empty_estimate <- function(exposure, domain, dd, msg) {
  tibble::tibble(
    infection = unname(hte_exposures[[exposure]]), exposure = exposure,
    domain = domain, n = nrow(dd), deaths = sum(dd$death30d == 1, na.rm = TRUE),
    n_negative = sum(as.character(dd[[exposure]]) == "0", na.rm = TRUE),
    n_positive = sum(as.character(dd[[exposure]]) == "1", na.rm = TRUE),
    deaths_negative = sum(dd$death30d == 1 & as.character(dd[[exposure]]) == "0", na.rm = TRUE),
    deaths_positive = sum(dd$death30d == 1 & as.character(dd[[exposure]]) == "1", na.rm = TRUE),
    beta_interaction = NA_real_, SE_interaction = NA_real_,
    OR_ratio = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
    p_wald = NA_real_, p_interaction_LRT = NA_real_,
    LRT_chisq = NA_real_, LRT_df = NA_real_,
    sparse_deaths = NA, note = msg
  )
}

# One entry for every domain and exposure in both cohort-eligible samples.
hte_results <- list()
hte_subgroup_results <- list()
hte_qc <- list()
hte_vifs <- list()

for (exposure in names(hte_exposures)) {
  required <- c("death30d", "cohort", hte_domain_terms, exposure)
  dd <- pooled[stats::complete.cases(pooled[, required, drop = FALSE]),
               , drop = FALSE]
  dd <- droplevels(dd)
  if (any(!dd$death30d %in% c(0, 1)))
    stop("Invalid mortality coding in infection interaction dataset: ", exposure)

  # Full inclusion audit, including number of deaths per infection stratum
  # within each cohort. Cohort/infection values are kept unchanged.
  for (this_cohort in levels(pooled$cohort)) {
    for (infection_value in c("0", "1")) {
      sub <- dd[as.character(dd$cohort) == this_cohort &
                  as.character(dd[[exposure]]) == infection_value, , drop = FALSE]
      hte_qc[[length(hte_qc) + 1L]] <- tibble::tibble(
        infection = unname(hte_exposures[[exposure]]),
        cohort = this_cohort,
        status = if (infection_value == "1") "Positive" else "Negative",
        n = nrow(sub), deaths = sum(sub$death30d == 1, na.rm = TRUE),
        survivors = sum(sub$death30d == 0, na.rm = TRUE)
      )
    }
  }

  if (nrow(dd) == 0L || length(unique(dd$death30d)) != 2L ||
      nlevels(dd$cohort) != 2L || nlevels(dd[[exposure]]) != 2L) {
    warning("Insufficient pooled exposure/cohort/outcome variation for ", exposure)
    for (domain in hte_domains) {
      hte_results[[length(hte_results) + 1L]] <-
        hte_empty_estimate(exposure, domain, dd, "Insufficient variation")
    }
    next
  }

  # Exactly the same baseline model for each of the three 1-df tests
  # within an infection type. Include three domain main effects, infection
  # status as a main effect, and cohort; no clinical covariates.
  null_formula <- stats::reformulate(
    c(hte_domain_terms, "cohort", exposure), response = "death30d")
  null_fit <- stats::glm(null_formula, data = dd, family = stats::binomial())

  for (domain in hte_domains) {
    term <- paste0(domain, "_z")
    full_formula <- stats::reformulate(
      c(hte_domain_terms, "cohort", exposure, paste0(term, ":", exposure)),
      response = "death30d")
    full_fit <- tryCatch(
      stats::glm(full_formula, data = dd, family = stats::binomial()),
      error = function(e) {
        warning("Interaction fit failed, ", exposure, " x ", domain, ": ",
                conditionMessage(e))
        NULL
      })

    if (is.null(full_fit)) {
      hte_results[[length(hte_results) + 1L]] <-
        hte_empty_estimate(exposure, domain, dd, "Full interaction model failed")
      next
    }

    # Factor levels 0/1 give an interaction coefficient for infection=1
    # relative to infection=0; we also check both term orderings.
    interaction_term <- paste0(term, ":", exposure, "1")
    alternate_term <- paste0(exposure, "1:", term)
    if (!interaction_term %in% names(stats::coef(full_fit)) &&
        alternate_term %in% names(stats::coef(full_fit)))
      interaction_term <- alternate_term

    cf <- stats::coef(full_fit)
    vv <- tryCatch(stats::vcov(full_fit), error = function(e) NULL)
    estimable <- !is.null(vv) &&
      all(c(term, interaction_term) %in% names(cf)) &&
      all(c(term, interaction_term) %in% rownames(vv)) &&
      all(is.finite(cf[c(term, interaction_term)])) &&
      all(is.finite(vv[c(term, interaction_term), c(term, interaction_term)]))
    if (!estimable) {
      hte_results[[length(hte_results) + 1L]] <-
        hte_empty_estimate(exposure, domain, dd, "Non-estimable interaction")
      next
    }

    comparison <- tryCatch(
      stats::anova(null_fit, full_fit, test = "Chisq"),
      error = function(e) NULL)
    lrt_p <- if (!is.null(comparison))
      suppressWarnings(as.numeric(comparison[2L, "Pr(>Chi)"])) else NA_real_
    lrt_chisq <- if (!is.null(comparison))
      suppressWarnings(as.numeric(comparison[2L, "Deviance"])) else NA_real_
    lrt_df <- if (!is.null(comparison))
      suppressWarnings(as.numeric(comparison[2L, "Df"])) else NA_real_

    b <- unname(cf[interaction_term]); se <- sqrt(unname(vv[interaction_term, interaction_term]))
    sparse <- sum(dd$death30d == 1 & as.character(dd[[exposure]]) == "1") < 5 ||
              sum(dd$death30d == 1 & as.character(dd[[exposure]]) == "0") < 5
    notes <- character(0)
    if (sparse) notes <- c(notes, "Fewer than 5 deaths in one infection stratum")
    if (!full_fit$converged) notes <- c(notes, "Full logistic model did not converge")
    if (!is.finite(se) || se > 5) notes <- c(notes, "Highly imprecise interaction coefficient")

    hte_results[[length(hte_results) + 1L]] <- tibble::tibble(
      infection = unname(hte_exposures[[exposure]]), exposure = exposure,
      domain = domain, n = stats::nobs(full_fit),
      deaths = sum(dd$death30d == 1),
      n_negative = sum(as.character(dd[[exposure]]) == "0"),
      n_positive = sum(as.character(dd[[exposure]]) == "1"),
      deaths_negative = sum(dd$death30d == 1 & as.character(dd[[exposure]]) == "0"),
      deaths_positive = sum(dd$death30d == 1 & as.character(dd[[exposure]]) == "1"),
      beta_interaction = b, SE_interaction = se,
      OR_ratio = exp(b), lower95 = exp(b - 1.96 * se),
      upper95 = exp(b + 1.96 * se),
      p_wald = 2 * stats::pnorm(abs(b / se), lower.tail = FALSE),
      p_interaction_LRT = lrt_p, LRT_chisq = lrt_chisq, LRT_df = lrt_df,
      sparse_deaths = sparse,
      note = paste(notes, collapse = "; ")
    )

    # Calculate both infection-stratum domain ORs from the SAME joint model.
    # In the reference stratum the domain effect is the main coefficient.
    # In the positive stratum it is main coefficient + interaction.
    for (status in c("Negative", "Positive")) {
      if (status == "Negative") {
        beta_s <- unname(cf[term])
        var_s <- unname(vv[term, term])
      } else {
        beta_s <- unname(cf[term] + cf[interaction_term])
        var_s <- unname(vv[term, term] + vv[interaction_term, interaction_term] +
                          2 * vv[term, interaction_term])
      }
      se_s <- if (is.finite(var_s) && var_s >= 0) sqrt(var_s) else NA_real_
      infection_value <- if (status == "Positive") "1" else "0"
      in_stratum <- as.character(dd[[exposure]]) == infection_value
      hte_subgroup_results[[length(hte_subgroup_results) + 1L]] <- tibble::tibble(
        infection = unname(hte_exposures[[exposure]]), domain = domain,
        status = status, n = sum(in_stratum),
        deaths = sum(dd$death30d == 1 & in_stratum),
        OR_per_SD = exp(beta_s), lower95 = exp(beta_s - 1.96 * se_s),
        upper95 = exp(beta_s + 1.96 * se_s),
        p_wald = 2 * stats::pnorm(abs(beta_s / se_s), lower.tail = FALSE),
        p_interaction_LRT = lrt_p
      )
    }

    # Collinearity diagnostics for the three domain main-effect predictors,
    # conditional on the cohort, exposure main effect and this interaction.
    # The additional interaction term can inflate these VIFs; interpret
    # them separately from VIFs for the primary joint mortality model.
    interaction_column <- paste0("hte_interaction_", domain)
    dd[[interaction_column]] <- dd[[term]] * as.integer(as.character(dd[[exposure]]) == "1")
    for (vif_term in hte_domain_terms) {
      other_terms <- setdiff(hte_domain_terms, vif_term)
      aux <- stats::lm(stats::reformulate(
        c(other_terms, "cohort", exposure, interaction_column),
        response = vif_term), data = dd)
      rsq <- summary(aux)$r.squared
      hte_vifs[[length(hte_vifs) + 1L]] <- tibble::tibble(
        infection = unname(hte_exposures[[exposure]]),
        interaction_domain = domain,
        domain = sub("_z$", "", vif_term),
        n = stats::nobs(aux), VIF = 1 / (1 - rsq)
      )
    }
  }
}

hte_interactions <- dplyr::bind_rows(hte_results) %>%
  dplyr::mutate(
    q_BH_all_9_interactions = stats::p.adjust(p_interaction_LRT, method = "BH"),
    OR_ratio_95CI = ifelse(
      is.finite(OR_ratio) & is.finite(lower95) & is.finite(upper95),
      sprintf("%.2f (%.2f–%.2f)", OR_ratio, lower95, upper95), NA_character_)
  )
hte_subgroups <- dplyr::bind_rows(hte_subgroup_results)
hte_availability <- dplyr::bind_rows(hte_qc)
hte_domain_vifs <- dplyr::bind_rows(hte_vifs)

utils::write.csv(hte_interactions,
  file.path(out_dir, "RESERVE_U_POOLED_DRT_infection_effect_modification_30d_9_single_interactions.csv"),
  row.names = FALSE)
utils::write.csv(hte_subgroups,
  file.path(out_dir, "RESERVE_U_POOLED_DRT_infection_effect_modification_subgroup_ORs.csv"),
  row.names = FALSE)
utils::write.csv(hte_availability,
  file.path(out_dir, "RESERVE_U_POOLED_DRT_infection_effect_modification_subgroup_counts.csv"),
  row.names = FALSE)
utils::write.csv(hte_domain_vifs,
  file.path(out_dir, "RESERVE_U_POOLED_DRT_infection_effect_modification_domain_VIFs.csv"),
  row.names = FALSE)

# One forest chart for the 9 ratios of ORs, with explicit null at 1.
# Missing/non-estimable coefficients remain in the CSV and are not plotted.
hte_plot <- hte_interactions %>%
  dplyr::filter(is.finite(OR_ratio), is.finite(lower95),
                is.finite(upper95), lower95 > 0, upper95 > 0) %>%
  dplyr::mutate(
    infection = factor(infection, levels = unname(hte_exposures)),
    domain = factor(domain, levels = rev(hte_domains)),
    p_label = ifelse(is.na(p_interaction_LRT),
                     "p = NA", paste0("p = ", sprintf("%.3g", p_interaction_LRT)))
  )
if (nrow(hte_plot)) {
  p_hte <- ggplot2::ggplot(hte_plot,
    ggplot2::aes(x = OR_ratio, y = domain, color = domain)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower95, xmax = upper95),
      orientation = "y", width = 0.18, linewidth = 0.9) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::geom_text(ggplot2::aes(x = upper95 * 1.12, label = p_label),
                       hjust = 0, size = 3.4, show.legend = FALSE) +
    ggplot2::facet_wrap(~infection, nrow = 1, scales = "free_x") +
    ggplot2::scale_color_manual(values = domain_colors) +
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.08, 0.72))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Exploratory infection-type modification of domain–mortality associations",
      x = "Interaction OR ratio: infection-positive vs negative (per 1-SD domain)",
      y = NULL) +
    ggprism::theme_prism(base_size = 12) +
    ggplot2::theme(legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.margin = ggplot2::margin(10, 110, 10, 10))
  save_journal_forest(p_hte,
    "RESERVE_U_POOLED_DRT_infection_effect_modification_30d_interaction_forest")
}

# Separate infection-specific joint-model conditional OR figures (three domains
# x two infection strata). Each domain row is also annotated with the formal
# 1-df likelihood-ratio interaction p-value from its corresponding joint model.
# The P-value tests DIFFERENCE BETWEEN STRATA, not an OR within either stratum.
# Annotations use a fixed right-hand column, outside the confidence intervals.
for (infection_name in unname(hte_exposures)) {
  sdf <- hte_subgroups %>%
    dplyr::filter(.data$infection == infection_name,
                  is.finite(OR_per_SD), is.finite(lower95),
                  is.finite(upper95), lower95 > 0, upper95 > 0) %>%
    dplyr::mutate(
      domain = factor(domain, levels = rev(hte_domains)),
      status = factor(status, levels = c("Negative", "Positive"))
    )
  if (!nrow(sdf)) next

  # One annotation per molecular domain. The p-values were computed above
  # by nested-model likelihood-ratio tests on each infection-specific sample.
  ann <- hte_interactions %>%
    dplyr::filter(.data$infection == infection_name,
                  .data$domain %in% as.character(sdf$domain)) %>%
    dplyr::mutate(
      domain = factor(domain, levels = rev(hte_domains)),
      label = ifelse(
        is.finite(p_interaction_LRT),
        paste0("Interaction p = ", sprintf("%.3g", p_interaction_LRT)),
        "Interaction p = NA"
      )
    )
  if (anyDuplicated(ann$domain))
    stop("Multiple interaction p-values for the same domain: ", infection_name)

  # Log-axis positions: create an explicit annotation column to the RIGHT
  # of every 95% CI so long confidence intervals never cover the p-values.
  ci_max <- max(sdf$upper95, na.rm = TRUE)
  ci_min <- min(sdf$lower95, na.rm = TRUE)
  annotation_x <- ci_max * 1.16
  x_left <- min(ci_min, 1) / 1.18
  x_right <- annotation_x * 2.25
  ann$annotation_x <- annotation_x

  p_subgroup <- ggplot2::ggplot(sdf, ggplot2::aes(
    x = OR_per_SD, y = domain, color = status, shape = status)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower95, xmax = upper95),
      orientation = "y", position = ggplot2::position_dodge(width = 0.45),
      width = 0.18, linewidth = 0.85) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.45), size = 3.2) +
    ggplot2::geom_text(
      data = ann,
      ggplot2::aes(x = annotation_x, y = domain, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 0.5,
      size = 3.75, fontface = "bold", color = "gray20"
    ) +
    ggplot2::scale_color_manual(values = c(Negative = "#2A6EBB", Positive = "#D65F44")) +
    ggplot2::scale_x_log10(limits = c(x_left, x_right)) +
    ggplot2::labs(
      title = paste0(infection_name, ": domain associations with 30-day mortality"),
      x = "Adjusted odds ratio for death per 1-SD domain score", y = NULL,
      color = "Infection status", shape = "Infection status") +
    ggprism::theme_prism(base_size = 12) +
    ggplot2::theme(legend.position = "top",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
  save_journal_forest(p_subgroup,
    paste0("RESERVE_U_POOLED_DRT_", infection_name,
           "_30d_joint_interaction_subgroup_OR_forest"))
}

cat("\n=== EXPLORATORY POOLED INFECTION x DOMAIN EFFECT MODIFICATION ===\n")
cat("Each of nine models: all three domain main effects + cohort + infection + ONE interaction.\n")
cat("No age/sex/illness-duration/other-infection adjustment in this block.\n")
print(hte_interactions)
cat("\nEvents and sample sizes by cohort and infection status:\n")
print(hte_availability)
cat("\nSubgroup-specific domain odds ratios:\n")
print(hte_subgroups)
cat("\nInteraction-model domain VIFs:\n")
print(hte_domain_vifs)


cat("\nJoint MEWS pooled and cohort estimates:\n")
print(mews_estimates)
cat("\nPooled adjusted infection-type associations:\n")
print(infection_associations %>% dplyr::filter(.data$source == "Pooled (cohort-adjusted)"))
cat("\nAll pooled mortality, MEWS, infection-type and effect-modification outputs saved to:\n", out_dir, "\n")


# ============================================================
# EXPLORATORY URINE LAM GRADE AND MULTI-OMIC D/R/T SIGNATURES
# ------------------------------------------------------------
# Only the ORIGINAL multi-omic PC1 domains are evaluated (not transcriptomic projections).
# The source score files are per-SD; the matching Resistance_z, Tolerance_z,
# and Damage_z fields were standardized using the complete multi-omic cohort
# (not the urine LAM-tested subset) in prepare_one_cohort() above.
#
# Clinical variable harmonization done above:
# RESERVE-U-2-TOR: urine_lam_done/result/band_grade
# RESERVE-U-1-EBB: urinetblam/result/grade
# Documented urine LAM-negative tests = grade 0; positives = grades 1-4.
# Untested participants are NOT classified as grade 0. Positive tests with
# missing band grade are excluded from ordinal-grade analyses.
# All participants with documented urine LAM results in these data are PLWH;
# therefore HIV status is not included as a covariate.
#
# Primary EXPLORATORY model for each domain, separately and pooled:
#   domain_z ~ lam_grade + age + sex + illness_duration + malaria [+ cohort]
# This follows the covariate strategy of the main infection-association block:
# age, sex, illness duration, and relevant concomitant infection (malaria).
# MEWS is NOT included. HIV does not vary among participants with LAM testing.
# microtbdx is NOT included since urine LAM positivity contributes to that
# composite TB diagnosis variable in both cohorts.
# HC3 robust SE, two-sided asymptotic z p-values, 95% Wald CI.
# Report NOMINAL p-values only; no multiplicity corrections or sensitivities
# are performed in this urine LAM section.
# This is an association, not a causal estimate of mycobacterial burden.
# ============================================================

lam_domain_names <- c("Resistance", "Tolerance", "Damage")
lam_covars <- c("age", "sex", "illness_duration", "malaria")

# Audit the LAM population and any positive tests without an observed grade.
lam_counts <- pooled %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarise(
    n_multiomic = dplyr::n(),
    n_lam_tested = sum(!is.na(lam_positive)),
    n_lam_negative = sum(lam_positive == 0, na.rm = TRUE),
    n_lam_positive = sum(lam_positive == 1, na.rm = TRUE),
    n_lam_grade_observed = sum(!is.na(lam_grade)),
    n_positive_missing_grade = sum(lam_positive == 1 & is.na(lam_grade), na.rm = TRUE),
    n_grade_0 = sum(lam_grade == 0, na.rm = TRUE),
    n_grade_1 = sum(lam_grade == 1, na.rm = TRUE),
    n_grade_2 = sum(lam_grade == 2, na.rm = TRUE),
    n_grade_3 = sum(lam_grade == 3, na.rm = TRUE),
    n_grade_4 = sum(lam_grade == 4, na.rm = TRUE),
    .groups = "drop"
  )
utils::write.csv(lam_counts,
  file.path(out_dir, "RESERVE_U_POOLED_DRT_LAM_sample_availability.csv"),
  row.names = FALSE)

# HC3 (MacKinnon-White) sandwich variance for lm without additional packages.
# Identical coefficient inference to lm + sandwich::vcovHC(type = "HC3").
lam_hc3_coefficient <- function(fit, term) {
  xx <- stats::model.matrix(fit)
  bb <- stats::coef(fit)
  if (fit$rank != ncol(xx) || !term %in% names(bb) ||
      !is.finite(bb[[term]]))
    stop("LAM model has a nonestimable coefficient: ", term)
  leverage <- stats::hatvalues(fit)
  if (any(!is.finite(leverage)) || any(leverage >= 1 - 1e-10))
    stop("LAM model has leverage near 1; HC3 variance is unstable")
  rr <- as.numeric(stats::residuals(fit)) / (1 - leverage)
  bread <- solve(crossprod(xx))
  meat <- crossprod(xx * rr)
  vc_hc3 <- bread %*% meat %*% bread
  se <- sqrt(unname(vc_hc3[term, term]))
  if (!is.finite(se) || se <= 0)
    stop("Invalid HC3 standard error for ", term)
  b <- unname(bb[[term]])
  z975 <- stats::qnorm(0.975)
  tibble::tibble(beta = b, lower95 = b - z975 * se,
                 upper95 = b + z975 * se,
                 p = 2 * stats::pnorm(abs(b / se), lower.tail = FALSE),
                 SE_HC3 = se)
}

lam_empty_row <- function(domain, source, dd, note) {
  tibble::tibble(
    source = source, domain = domain, exposure = "lam_grade",
    n = nrow(dd), n_negative = sum(dd$lam_positive == 0, na.rm = TRUE),
    n_positive = sum(dd$lam_positive == 1, na.rm = TRUE),
    beta = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
    p = NA_real_, SE_HC3 = NA_real_, R2 = NA_real_, note = note)
}

run_lam_association <- function(data, domain, source_label,
                                add_cohort = FALSE) {
  stopifnot(domain %in% lam_domain_names)
  outcome <- paste0(domain, "_z")
  adjustment <- c(lam_covars, if (add_cohort) "cohort")
  needed <- c(outcome, "lam_grade", adjustment)
  dd <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  dd <- droplevels(dd)
  if (nrow(dd) < 15L || dplyr::n_distinct(dd$lam_grade) < 2L ||
      (add_cohort && dplyr::n_distinct(dd$cohort) < 2L))
    return(lam_empty_row(domain, source_label, dd,
                         "Insufficient grade or cohort variation"))

  form <- stats::reformulate(c("lam_grade", adjustment), response = outcome)
  res <- tryCatch({
    fit <- stats::lm(form, data = dd)
    dplyr::bind_cols(lam_hc3_coefficient(fit, "lam_grade"),
                     tibble::tibble(R2 = summary(fit)$r.squared))
  }, error = function(e) {
    warning("LAM grade model failed for ", source_label, " / ", domain,
            ": ", conditionMessage(e))
    NULL
  })
  if (is.null(res))
    return(lam_empty_row(domain, source_label, dd,
                         "Nonestimable HC3 model"))

  dplyr::bind_cols(
    tibble::tibble(
      source = source_label, domain = domain, exposure = "lam_grade",
      n = nrow(dd),
      n_negative = sum(dd$lam_positive == 0, na.rm = TRUE),
      n_positive = sum(dd$lam_positive == 1, na.rm = TRUE)),
    res, tibble::tibble(note = ""))
}

# The three cohorts/source groups must have the same variable specification.
lam_sources <- list(
  "RESERVE-U-2-TOR" = tor,
  "RESERVE-U-1-EBB" = ebb,
  "Pooled (cohort-adjusted)" = pooled
)
lam_grade_results <- purrr::map_dfr(lam_domain_names, function(domain) {
  purrr::imap_dfr(lam_sources, function(dat, source_label) {
    run_lam_association(
      dat, domain, source_label = source_label,
      add_cohort = identical(source_label, "Pooled (cohort-adjusted)"))
  })
})
utils::write.csv(lam_grade_results,
  file.path(out_dir, "RESERVE_U_POOLED_DRT_LAM_grade_cohort_and_pooled_HC3.csv"),
  row.names = FALSE)

# ------------------------------------------------------------
# FIGURE: Adjusted urine LAM grade regression forest.
# Reuse the existing journal_forest() layout: for EACH domain, cohort-specific
# and pooled estimates, 95% HC3 CIs, and nominal p-values. No transcriptomic projections.
# ------------------------------------------------------------
if (any(!is.finite(lam_grade_results$beta))) {
  warning("Some LAM grade coefficients were not estimable; verify CSV first")
} else {
  p_lam_grade_forest <- journal_forest(
    lam_grade_results,
    figure_title = "Urine LAM grade and multi-omic host-response domains",
    axis_title = "Adjusted change in domain score (SD) per one-grade increase",
    measure = "beta"
  )
  save_journal_forest(p_lam_grade_forest,
    "RESERVE_U_POOLED_DRT_LAM_grade_adjusted_cohort_and_pooled_forest")
}

cat("\n=== EXPLORATORY URINE LAM GRADE / MULTI-OMIC D/R/T ANALYSES ===\n")
cat("Urine LAM sample availability by cohort:\n")
print(lam_counts)
cat("\nCohort-specific and pooled HC3-adjusted ordinal grade results:\n")
print(lam_grade_results)
cat("\nLAM association results CSV and adjusted forest PDF saved to:\n", out_dir, "\n")
