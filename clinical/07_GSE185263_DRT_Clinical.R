# =============================================================================
# GSE185263: Resistance / Tolerance / Damage (domain) clinical analyses
# Script: 07_GSE185263_DRT_Clinical.R
# =============================================================================
#
# SCORES
#   Locked 14-mediator transcriptomic projections (sepsis-standardized).
#   No MoSS inputs or models. The three domains enter adjusted clinical-outcome
#   models SIMULTANEOUSLY, controlling for one another.
#
# ANALYSES
#   (1) Healthy vs sepsis: unadjusted Wilcoxon comparisons by domain.
#   (2) In-hospital mortality in sepsis: joint logistic regression with all three
#       domain domains + age, sex, and recruitment location fixed effects.
#       Descriptive predicted-risk curves vary one domain at a time; other domains
#       are held at 0 and prediction is averaged over observed site frequencies.
#       Exploratory domain-by-site interaction LRTs appear in tables only.
#   (3) 24-hour SOFA in sepsis: joint median (tau=0.5) quantile regression with
#       the same adjustment set; site-weighted curves with bootstrap CIs.
#       Domain VIFs exported separately for the mortality and SOFA complete cases.
#   (4) Sex comparisons: unadjusted Wilcoxon and adjusted HC3-robust linear
#       models, overall (age, disease state, site) and sepsis-only (age, site).
#
# NOTES
#   All transcriptomic projection plots saved as PDF only; tables exported as CSV.
#   No score restandardization: score-file units are kept.
#   Mortality is IN-HOSPITAL, not 28-day or 30-day mortality.
# =============================================================================

rm(list = ls())

required_pkgs <- c("dplyr", "tidyr", "tibble", "purrr", "ggplot2", "ggprism",
                   "MASS", "quantreg", "sandwich", "lmtest")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                     logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install required R packages: install.packages(c(",
       paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(ggprism)
})

# ------------------------------------------------------------
# INPUT / OUTPUT PATHS (relative to the project root; edit if needed).
# Scores are produced by publication_code/01_DRT_scores.R.
# ------------------------------------------------------------
score_file    <- "DRT_scores_output/GSE185263/GSE185263_DRT_scores_wide_perSD.csv"
metadata_file <- "Multinational_sepsis_GSE185263/GSE185263_metadata.csv"

out_dir <- "DRT_scores_output/Clinical analyses/GSE185263"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(out_dir)) stop("Output directory unavailable: ", out_dir)

# Settings for primary plot uncertainty: reproducible simulation / bootstrap.
N_LOGISTIC_DRAWS <- 500L
N_SOFA_BOOTSTRAPS <- 500L
RUN_SITE_INTERACTION_SENSITIVITY <- TRUE
set.seed(185263)

domains <- c("Resistance", "Tolerance", "Damage")
# A single, prespecified scoring method is used throughout this analysis.
method_key <- tibble::tibble(
  method = "Transcriptomic projection",
  prefix = "rna_surrogate_lasso14_1se__"
)
domain_colors <- c(Resistance = "#2A6EBB", Tolerance = "#2CA25F",
                   Damage = "#E6550D")

save_csv <- function(x, stem) {
  utils::write.csv(x, file.path(out_dir, paste0(stem, ".csv")),
                   row.names = FALSE, na = "")
}

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

format_p <- function(p) {
  vapply(p, function(v) {
    if (is.na(v) || !is.finite(v)) return("p = NA")
    if (v == 0) return("p = 0 (numerical underflow)")
    paste0("p = ", sprintf("%.3g", v))
  }, character(1))
}
check_cols <- function(x, required, source) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(source, " is missing: ", paste(missing, collapse = ", "))
}

# ----------------------------------------------------------------------------
# 1. Input metadata: use definitive GEO columns, not parsed free-text fields.
# ----------------------------------------------------------------------------
scores <- utils::read.csv(score_file, check.names = FALSE,
                          stringsAsFactors = FALSE)
meta <- utils::read.csv(metadata_file, check.names = FALSE,
                        stringsAsFactors = FALSE)
score_cols <- unlist(lapply(method_key$prefix, function(prefix)
  paste0(prefix, domains)), use.names = FALSE)
check_cols(scores, c("sample_id", score_cols), "GSE185263 score CSV")
check_cols(meta, c("geo_accession", "age:ch1", "Sex:ch1",
                   "collection location:ch1", "disease state:ch1",
                   "in hospital mortality:ch1", "sofa 24h post admisssion:ch1"),
           "GSE185263 metadata CSV")

scores$sample_id <- trimws(as.character(scores$sample_id))
meta$geo_accession <- trimws(as.character(meta$geo_accession))
if (anyNA(scores$sample_id) || anyNA(meta$geo_accession) ||
    anyDuplicated(scores$sample_id) || anyDuplicated(meta$geo_accession)) {
  stop("Missing or duplicate GEO identifiers in one of the source files")
}

# Check potentially duplicated clinical fields in the score file, without
# allowing their names to create .x/.y suffixes during the metadata join.
shared_checks <- c("disease state:ch1", "in hospital mortality:ch1")
for (field in shared_checks) {
  if (field %in% names(scores)) {
    ix <- match(scores$sample_id, meta$geo_accession)
    a <- trimws(as.character(scores[[field]]))
    b <- trimws(as.character(meta[[field]][ix]))
    inconsistent <- !is.na(ix) & !is.na(a) & !is.na(b) & a != b
    if (any(inconsistent)) stop("Score/metadata discrepancy for ", field,
                               ": ", sum(inconsistent), " records")
  }
}

clinical <- meta %>%
  dplyr::transmute(
    sample_id = geo_accession,
    age = suppressWarnings(as.numeric(.data[["age:ch1"]])),
    sex = toupper(trimws(as.character(.data[["Sex:ch1"]]))),
    disease_state = tolower(trimws(as.character(.data[["disease state:ch1"]]))),
    RecruitmentSite_raw = tolower(trimws(as.character(.data[["collection location:ch1"]]))),
    mortality_raw = trimws(as.character(.data[["in hospital mortality:ch1"]])),
    sofa24 = suppressWarnings(as.numeric(.data[["sofa 24h post admisssion:ch1"]]))
  ) %>%
  dplyr::mutate(
    sex = dplyr::case_when(sex %in% c("M", "MALE") ~ "M",
                            sex %in% c("F", "FEMALE") ~ "F", TRUE ~ NA_character_),
    sex = factor(sex, levels = c("F", "M")),
    disease_state = dplyr::case_when(
      disease_state %in% c("sepsis", "infected", "infection") ~ "Sepsis",
      disease_state %in% c("healthy", "healthy control", "control") ~ "Healthy",
      TRUE ~ NA_character_),
    disease_state = factor(disease_state, levels = c("Healthy", "Sepsis")),
    # Replicates the original MoSS workflow's Toronto/Vancouver consolidation.
    RecruitmentSite = dplyr::case_when(
      RecruitmentSite_raw %in% c("toronto", "vancouver") ~ "canada",
      RecruitmentSite_raw %in% c("", "unknown") ~ NA_character_,
      TRUE ~ RecruitmentSite_raw),
    mortality01 = dplyr::case_when(
      tolower(mortality_raw) %in% c("died", "dead", "1") ~ 1L,
      tolower(mortality_raw) %in% c("survived", "alive", "0") ~ 0L,
      TRUE ~ NA_integer_)
  )

score_values <- scores %>%
  dplyr::select(sample_id, dplyr::all_of(score_cols)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(score_cols),
                              ~ suppressWarnings(as.numeric(.x))))
dat <- dplyr::inner_join(clinical, score_values, by = "sample_id")
if (!nrow(dat)) stop("No sample IDs matched between files")
if (nrow(dat) != nrow(scores) || nrow(dat) != nrow(meta)) {
  warning("Not all input samples matched; inspect the QC table")
}
if (any(is.na(dat$disease_state))) stop("Unexpected missing/unknown disease state")

sample_qc <- tibble::tibble(
  metric = c("Score-file rows", "Metadata rows", "Matched samples",
             "Healthy controls", "Sepsis participants", "Sepsis with known mortality",
             "Sepsis in-hospital deaths", "Sepsis with 24h SOFA"),
  n = c(nrow(scores), nrow(meta), nrow(dat),
        sum(dat$disease_state == "Healthy"),
        sum(dat$disease_state == "Sepsis"),
        sum(dat$disease_state == "Sepsis" & !is.na(dat$mortality01)),
        sum(dat$disease_state == "Sepsis" & dat$mortality01 == 1, na.rm = TRUE),
        sum(dat$disease_state == "Sepsis" & !is.na(dat$sofa24)))
)
save_csv(sample_qc, "GSE185263_DRT_sample_availability_QC")
save_csv(method_key, "GSE185263_DRT_score_methods")
save_csv(dat %>% dplyr::count(disease_state, RecruitmentSite, name = "n"),
         "GSE185263_DRT_recruitment_site_counts")
print(sample_qc)

# Retain original perSD values; never restandardize controls or sepsis scores.
wide_method <- function(one_method) {
  prefix <- method_key$prefix[match(one_method, method_key$method)]
  dd <- dat %>% dplyr::select(sample_id, age, sex, disease_state,
                             RecruitmentSite, mortality01, sofa24)
  for (domain in domains) dd[[domain]] <- dat[[paste0(prefix, domain)]]
  dd$method <- one_method
  dd
}
method_data <- lapply(method_key$method, wide_method)
names(method_data) <- method_key$method

# ----------------------------------------------------------------------------
# 2. Healthy controls vs sepsis: unadjusted Wilcoxon per domain.
# ----------------------------------------------------------------------------
baseline_comparison <- purrr::map_dfr(method_key$method, function(m) {
  dd <- method_data[[m]]
  purrr::map_dfr(domains, function(d) {
    x <- dd[[d]][dd$disease_state == "Sepsis" & is.finite(dd[[d]])]
    y <- dd[[d]][dd$disease_state == "Healthy" & is.finite(dd[[d]])]
    if (!length(x) || !length(y)) stop("Missing groups for ", m, ": ", d)
    test <- stats::wilcox.test(x, y, exact = FALSE)
    tibble::tibble(method = m, domain = d, n_sepsis = length(x), n_healthy = length(y),
                   median_sepsis = stats::median(x), median_healthy = stats::median(y),
                   q1_sepsis = unname(stats::quantile(x, .25)),
                   q3_sepsis = unname(stats::quantile(x, .75)),
                   q1_healthy = unname(stats::quantile(y, .25)),
                   q3_healthy = unname(stats::quantile(y, .75)),
                   W = unname(test$statistic), p = test$p.value)
  })
}) %>%
  dplyr::group_by(method) %>%
  dplyr::mutate(q_BH_3domains = stats::p.adjust(p, method = "BH")) %>%
  dplyr::ungroup()
save_csv(baseline_comparison, "GSE185263_DRT_healthy_vs_sepsis_Wilcoxon")

rna <- method_data[["Transcriptomic projection"]]
plot_disease <- rna %>%
  tidyr::pivot_longer(cols = dplyr::all_of(domains),
                      names_to = "domain", values_to = "score") %>%
  dplyr::filter(is.finite(score)) %>%
  dplyr::mutate(domain = factor(domain, levels = domains))
disease_labels <- baseline_comparison %>%
  dplyr::filter(method == "Transcriptomic projection") %>%
  dplyr::mutate(domain = factor(domain, levels = domains),
                label = paste0("Wilcoxon ", format_p(p)))
p_disease <- ggplot2::ggplot(plot_disease,
                             ggplot2::aes(x = disease_state, y = score,
                                          fill = disease_state, color = disease_state)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.58,
                        linewidth = 0.65, color = "black", alpha = 0.46) +
  ggplot2::geom_jitter(width = 0.15, height = 0,
                       size = 1.5, alpha = 0.48) +
  ggplot2::geom_text(data = disease_labels,
    ggplot2::aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE,
    vjust = 1.25, fontface = "bold", size = 3.6) +
  ggplot2::facet_wrap(~domain, nrow = 1, scales = "fixed") +
  ggplot2::scale_fill_manual(values = c(Healthy = "#A2A2A2", Sepsis = "#2B8C89")) +
  ggplot2::scale_color_manual(values = c(Healthy = "#777777", Sepsis = "#1A5354")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.06, .20))) +
  ggplot2::labs(title = "GSE185263: domain scores in sepsis and healthy controls",
                x = NULL, y = "Transcriptomic projection domain score (sepsis SD units)") +
  ggprism::theme_prism(base_size = 13) +
  ggplot2::theme(legend.position = "none",
                 strip.text = ggplot2::element_text(face = "bold"),
                 plot.title = ggplot2::element_text(face = "bold", hjust = .5))
save_plot(p_disease, "GSE185263_RNA_LASSO_healthy_vs_sepsis_boxplots")

# ----------------------------------------------------------------------------
# 3. In-hospital mortality: JOINT 3-domain logistic regression.
#    Site fixed effects, age, sex; confidence intervals are Wald 95%.
# ----------------------------------------------------------------------------
clinical_adjustment <- c("age", "sex", "RecruitmentSite")
model_predictors <- c(domains, clinical_adjustment)
logit_formula <- stats::reformulate(model_predictors, response = "mortality01")
sofa_formula <- stats::reformulate(model_predictors, response = "sofa24")

extract_domain_model <- function(fit, one_method, outcome) {
  b <- stats::coef(fit)[domains]
  v <- stats::vcov(fit)
  se <- sqrt(diag(v))[domains]
  if (any(!is.finite(b)) || any(!is.finite(se)) || any(se <= 0))
    stop("Nonestimable domain coefficient in ", one_method, " / ", outcome)
  p <- 2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
  tibble::tibble(method = one_method, outcome = outcome,
                 domain = domains, n = stats::nobs(fit),
                 log_OR = unname(b), SE = unname(se),
                 OR = exp(unname(b)),
                 lower95 = exp(unname(b - 1.96 * se)),
                 upper95 = exp(unname(b + 1.96 * se)),
                 p = unname(p), q_BH_3domains = stats::p.adjust(p, "BH"))
}

mortality_fits <- lapply(method_key$method, function(m) {
  dd <- method_data[[m]] %>%
    dplyr::filter(disease_state == "Sepsis") %>%
    dplyr::filter(stats::complete.cases(dplyr::across(
      dplyr::all_of(c("mortality01", model_predictors))))) %>%
    dplyr::mutate(RecruitmentSite = factor(RecruitmentSite)) %>% droplevels()
  if (nrow(dd) < 40L || dplyr::n_distinct(dd$mortality01) < 2L)
    stop("Insufficient mortality observations: ", m)
  fit <- stats::glm(logit_formula, data = dd, family = stats::binomial())
  if (!fit$converged) warning("Mortality GLM did not converge for ", m)
  result <- extract_domain_model(fit, m, "In-hospital mortality") %>%
    dplyr::mutate(deaths = sum(dd$mortality01 == 1L),
                  survivors = sum(dd$mortality01 == 0L))
  list(data = dd, fit = fit, estimates = result)
})
names(mortality_fits) <- method_key$method
mortality_table <- dplyr::bind_rows(lapply(mortality_fits, `[[`, "estimates"))
save_csv(mortality_table, "GSE185263_DRT_joint_inhospital_mortality_logistic")

# Explore one-domain-at-a-time interactions with recruitment site, adjusting
# for the other two domains (sensitivity; never used to select primary model).
site_interactions <- tibble::tibble()
if (RUN_SITE_INTERACTION_SENSITIVITY) {
  site_interactions <- purrr::map_dfr(method_key$method, function(m) {
    dd <- mortality_fits[[m]]$data
    base_fit <- mortality_fits[[m]]$fit
    purrr::map_dfr(domains, function(d) {
      int_form <- stats::as.formula(paste0(
        "mortality01 ~ ", paste(domains, collapse = " + "),
        " + age + sex + RecruitmentSite + ", d, ":RecruitmentSite"))
      fi <- tryCatch(suppressWarnings(stats::glm(
        int_form, data = dd, family = stats::binomial())), error = function(e) NULL)
      if (is.null(fi) || !fi$converged) {
        return(tibble::tibble(method = m, domain = d, n = nrow(dd),
                              degrees_freedom = NA_real_, p_LRT = NA_real_,
                              note = "Interaction fit unavailable/nonconvergent"))
      }
      lrt <- stats::anova(base_fit, fi, test = "Chisq")
      tibble::tibble(method = m, domain = d, n = nrow(dd),
                     degrees_freedom = as.numeric(lrt[2L, "Df"]),
                     p_LRT = as.numeric(lrt[2L, "Pr(>Chi)"]),
                     note = "Exploratory domain-by-site likelihood-ratio test")
    })
  }) %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(q_BH_3domains = stats::p.adjust(p_LRT, method = "BH")) %>%
    dplyr::ungroup()
  save_csv(site_interactions, "GSE185263_DRT_mortality_domain_by_site_interactions_exploratory")

}

# Tidy regression forest plot for the locked transcriptomic projection scores.
forest_mort <- mortality_table %>%
  dplyr::filter(method == "Transcriptomic projection") %>%
  dplyr::mutate(domain = factor(domain, levels = rev(domains)),
                annotation_x = max(upper95) * 1.10,
                annotation = sprintf("OR %.2f (%.2f–%.2f)\n%s",
                                     OR, lower95, upper95, format_p(p)))
p_mort_forest <- ggplot2::ggplot(forest_mort,
  ggplot2::aes(x = OR, y = domain, color = domain)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = lower95, xmax = upper95),
                         orientation = "y", width = .15, linewidth = .85) +
  ggplot2::geom_point(size = 3.5) +
  ggplot2::geom_text(ggplot2::aes(x = annotation_x, label = annotation),
                     hjust = 0, fontface = "bold", size = 3.5,
                     show.legend = FALSE) +
  ggplot2::scale_color_manual(values = domain_colors) +
  ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(.08, .78))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(title = "GSE185263: joint domain and in-hospital mortality",
                x = "Adjusted odds ratio per 1-SD unit", y = NULL) +
  ggprism::theme_prism(base_size = 13) +
  ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(face = "bold"),
                 axis.text.y = ggplot2::element_text(face = "bold"),
                 plot.margin = ggplot2::margin(10, 155, 10, 10))
save_plot(p_mort_forest, "GSE185263_RNA_LASSO_joint_inhospital_mortality_forest")

# Helper: prepare predictions for each domain at 5th–95th percentiles of its
# observed sepsis distribution, at median age / modal sex, weighted across
# actual site proportions. The other two domain scores are fixed at 0 SD.
make_site_grid <- function(dd, n_per_domain = 101L) {
  sites <- levels(dd$RecruitmentSite)
  site_w <- dd %>% dplyr::count(RecruitmentSite, name = "site_n") %>%
    dplyr::mutate(w = site_n / sum(site_n)) %>%
    dplyr::select(RecruitmentSite, w)
  modal_sex <- names(sort(table(dd$sex), decreasing = TRUE))[1L]
  out <- purrr::map_dfr(domains, function(d) {
    x <- seq(stats::quantile(dd[[d]], .05),
             stats::quantile(dd[[d]], .95), length.out = n_per_domain)
    nd <- tidyr::expand_grid(grid_ix = seq_along(x),
                             RecruitmentSite = factor(sites, levels = sites)) %>%
      dplyr::mutate(domain = d, score_x = x[grid_ix],
                    age = stats::median(dd$age),
                    sex = factor(modal_sex, levels = levels(dd$sex)))
    for (domain in domains) nd[[domain]] <- if (domain == d) nd$score_x else 0
    dplyr::left_join(nd, site_w, by = "RecruitmentSite")
  })
  out$domain <- factor(out$domain, levels = domains)
  if (anyNA(out$w)) stop("Missing site weights in prediction grid")
  out
}

# One index per (domain, grid point), each containing the rows for all sites.
site_group_index <- function(grid) {
  split(seq_len(nrow(grid)), interaction(grid$domain, grid$grid_ix,
                                         drop = TRUE, lex.order = TRUE))
}
site_weight_predictions <- function(pred, grid, groups) {
  vapply(groups, function(ix) sum(pred[ix] * grid$w[ix]), numeric(1))
}
site_curve_key <- function(grid, groups) {
  grid[vapply(groups, `[`, integer(1), 1L), c("domain", "grid_ix", "score_x")]
}

# For confidence bounds, draw model coefficients jointly from their estimated
# covariance, marginalize each draw over sites, THEN take pointwise quantiles.
# This avoids averaging site-specific CI endpoints (not a valid pooled CI).
fit_m <- mortality_fits[["Transcriptomic projection"]]
mort_grid <- make_site_grid(fit_m$data)
mort_groups <- site_group_index(mort_grid)
Xm <- stats::model.matrix(stats::delete.response(stats::terms(fit_m$fit)),
                          data = mort_grid)
beta <- stats::coef(fit_m$fit)
vc <- stats::vcov(fit_m$fit)
if (any(!is.finite(beta)) || any(!is.finite(vc))) stop("Mortality covariance invalid")
stopifnot(identical(colnames(Xm), names(beta)))
pred_mort <- stats::plogis(drop(Xm %*% beta))
mort_curve <- site_curve_key(mort_grid, mort_groups)
mort_curve$probability <- site_weight_predictions(pred_mort, mort_grid, mort_groups)
set.seed(185264)
bdraws <- MASS::mvrnorm(n = N_LOGISTIC_DRAWS, mu = beta, Sigma = vc)
curve_draws <- matrix(NA_real_, nrow = nrow(mort_curve), ncol = N_LOGISTIC_DRAWS)
for (b in seq_len(N_LOGISTIC_DRAWS)) {
  site_pred <- stats::plogis(drop(Xm %*% bdraws[b, ]))
  curve_draws[, b] <- site_weight_predictions(site_pred, mort_grid, mort_groups)
}
mort_curve$lower95 <- apply(curve_draws, 1L, stats::quantile, probs = .025)
mort_curve$upper95 <- apply(curve_draws, 1L, stats::quantile, probs = .975)
save_csv(mort_curve, "GSE185263_RNA_LASSO_joint_inhospital_mortality_siteweighted_predictions")

mort_ann <- mortality_table %>%
  dplyr::filter(method == "Transcriptomic projection") %>%
  dplyr::mutate(domain = factor(domain, levels = domains),
                label = sprintf("OR %.2f (%.2f–%.2f)\n%s",
                                OR, lower95, upper95, format_p(p)))
p_mort_curve <- ggplot2::ggplot(mort_curve,
                                 ggplot2::aes(x = score_x, y = probability,
                                              color = domain, fill = domain)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower95, ymax = upper95),
                       alpha = .18, color = NA) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_text(data = mort_ann,
    ggplot2::aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -.04, vjust = 1.15,
    fontface = "bold", size = 3.5) +
  ggplot2::facet_wrap(~domain, nrow = 1, scales = "free_x") +
  ggplot2::scale_color_manual(values = domain_colors) +
  ggplot2::scale_fill_manual(values = domain_colors) +
  ggplot2::scale_y_continuous(limits = c(0, 1),
                              expand = ggplot2::expansion(mult = c(0, .02))) +
  ggplot2::labs(title = "GSE185263: adjusted in-hospital mortality by transcriptomic projection domain",
                x = "Domain score (SD units)",
                y = "Adjusted probability of in-hospital death") +
  ggprism::theme_prism(base_size = 13) +
  ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(face = "bold"),
                 plot.title = ggplot2::element_text(face = "bold", hjust = .5),
                 plot.caption = ggplot2::element_text(size = 9, hjust = 0))
save_plot(p_mort_curve, "GSE185263_RNA_LASSO_joint_inhospital_mortality_probability_curves")

# ----------------------------------------------------------------------------
# 4. SOFA at 24 hours: JOINT median quantile regression, age/sex/site adjusted.
# ----------------------------------------------------------------------------
fit_joint_sofa <- function(m) {
  dd <- method_data[[m]] %>%
    dplyr::filter(disease_state == "Sepsis") %>%
    dplyr::filter(stats::complete.cases(dplyr::across(
      dplyr::all_of(c("sofa24", model_predictors))))) %>%
    dplyr::mutate(RecruitmentSite = factor(RecruitmentSite)) %>% droplevels()
  if (nrow(dd) < 40L) stop("Insufficient 24-hour SOFA observations: ", m)
  fit <- quantreg::rq(sofa_formula, data = dd, tau = .5)
  sm <- tryCatch(summary(fit, se = "nid"), error = function(e) {
    warning("Quantile NID standard errors failed for ", m, ": ", conditionMessage(e),
            ". Falling back to bootstrap coefficient SE.")
    summary(fit, se = "boot", R = 499L)
  })
  coefs <- as.matrix(sm$coefficients)
  if (!all(domains %in% rownames(coefs))) stop("SOFA quantile coefficients missing")
  b <- coefs[domains, 1L]
  se <- coefs[domains, 2L]
  if (any(!is.finite(b)) || any(!is.finite(se)) || any(se <= 0))
    stop("Invalid SOFA quantile regression SE for ", m)
  # The conventional summary(rq) fourth column is its two-sided p-value.
  p <- if (ncol(coefs) >= 4L) coefs[domains, 4L] else
    2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
  out <- tibble::tibble(method = m, outcome = "24-hour SOFA (median)",
    domain = domains, n = nrow(dd), beta = unname(b), SE = unname(se),
    lower95 = unname(b - 1.96 * se), upper95 = unname(b + 1.96 * se),
    p = unname(p), q_BH_3domains = stats::p.adjust(p, "BH"),
    tau = .5, adjusted_for = "Other DRT domains, age, sex, recruitment location")
  list(data = dd, fit = fit, estimates = out)
}
sofa_fits <- lapply(method_key$method, fit_joint_sofa)
names(sofa_fits) <- method_key$method
sofa_table <- dplyr::bind_rows(lapply(sofa_fits, `[[`, "estimates"))
save_csv(sofa_table, "GSE185263_DRT_joint_24h_SOFA_quantile_regression")

# ------------------------------------------------------------
# Joint-model domain VIFs, computed separately on each model's ACTUAL
# complete-case population, using the same three domains + age/sex/site.
# VIF uses auxiliary linear regressions of each domain on all other predictors;
# this applies equally to GLM and median quantile-regression design matrices.
# ------------------------------------------------------------
calculate_joint_vif <- function(dd, outcome, method) {
  purrr::map_dfr(domains, function(domain) {
    others <- setdiff(model_predictors, domain)
    fit_aux <- stats::lm(stats::reformulate(others, response = domain), data = dd)
    r2 <- summary(fit_aux)$r.squared
    tibble::tibble(method = method, outcome = outcome, domain = domain,
                   n = stats::nobs(fit_aux), VIF = 1 / (1 - r2))
  })
}
joint_vif <- dplyr::bind_rows(
  purrr::map_dfr(method_key$method, function(m)
    calculate_joint_vif(mortality_fits[[m]]$data, "In-hospital mortality", m)),
  purrr::map_dfr(method_key$method, function(m)
    calculate_joint_vif(sofa_fits[[m]]$data, "24-hour SOFA (median)", m))
)
save_csv(joint_vif, "GSE185263_DRT_joint_mortality_SOFA_domain_VIF")

# Bootstrap the entire joint QR, resampling participants WITHIN each site to
# preserve fixed-effect categories in every replicate. Compute the site-weighted
# curve within each replicate and percentile CIs across those curves.
fit_s <- sofa_fits[["Transcriptomic projection"]]
sofa_grid <- make_site_grid(fit_s$data)
sofa_groups <- site_group_index(sofa_grid)
sofa_curve <- site_curve_key(sofa_grid, sofa_groups)
sofa_curve$median_sofa <- site_weight_predictions(
  as.numeric(stats::predict(fit_s$fit, newdata = sofa_grid)),
  sofa_grid, sofa_groups)

set.seed(185265)
site_indices <- split(seq_len(nrow(fit_s$data)), fit_s$data$RecruitmentSite,
                      drop = TRUE)
sofa_boot <- matrix(NA_real_, nrow = nrow(sofa_curve),
                    ncol = N_SOFA_BOOTSTRAPS)
successes <- 0L
for (b in seq_len(N_SOFA_BOOTSTRAPS)) {
  sampled <- unlist(lapply(site_indices, function(idx)
    sample(idx, size = length(idx), replace = TRUE)), use.names = FALSE)
  boot_dat <- fit_s$data[sampled, , drop = FALSE]
  fb <- tryCatch(suppressWarnings(quantreg::rq(
    sofa_formula, data = boot_dat, tau = .5)), error = function(e) NULL)
  if (is.null(fb)) next
  pr <- tryCatch(as.numeric(stats::predict(fb, newdata = sofa_grid)),
                 error = function(e) NULL)
  if (is.null(pr) || any(!is.finite(pr))) next
  successes <- successes + 1L
  sofa_boot[, b] <- site_weight_predictions(pr, sofa_grid, sofa_groups)
}
cat("\nSuccessful SOFA bootstrap replicates:", successes,
    "of", N_SOFA_BOOTSTRAPS, "\n")
if (successes < 50L) stop("Too few successful SOFA bootstrap replicates (<50)")
if (successes < .8 * N_SOFA_BOOTSTRAPS)
  warning("Fewer than 80% successful SOFA bootstrap replicates")
sofa_curve$lower95 <- apply(sofa_boot, 1L, stats::quantile,
                            probs = .025, na.rm = TRUE)
sofa_curve$upper95 <- apply(sofa_boot, 1L, stats::quantile,
                            probs = .975, na.rm = TRUE)
sofa_curve$bootstrap_successes <- successes
save_csv(sofa_curve, "GSE185263_RNA_LASSO_joint_24h_SOFA_siteweighted_median_curves")

raw_sofa_long <- fit_s$data %>%
  tidyr::pivot_longer(cols = dplyr::all_of(domains),
                      names_to = "domain", values_to = "score_x") %>%
  dplyr::mutate(domain = factor(domain, levels = domains))
sofa_ann <- sofa_table %>%
  dplyr::filter(method == "Transcriptomic projection") %>%
  dplyr::mutate(domain = factor(domain, levels = domains),
    label = sprintf("Adjusted β %.2f (%.2f–%.2f)\n%s",
                    beta, lower95, upper95, format_p(p)))
p_sofa <- ggplot2::ggplot() +
  ggplot2::geom_point(data = raw_sofa_long,
    ggplot2::aes(x = score_x, y = sofa24, color = domain),
    alpha = .24, size = 1.25) +
  ggplot2::geom_ribbon(data = sofa_curve,
    ggplot2::aes(x = score_x, ymin = lower95, ymax = upper95, fill = domain),
    alpha = .18) +
  ggplot2::geom_line(data = sofa_curve,
    ggplot2::aes(x = score_x, y = median_sofa, color = domain),
    linewidth = 1.15) +
  ggplot2::geom_text(data = sofa_ann,
    ggplot2::aes(x = -Inf, y = Inf, label = label),
    hjust = -.04, vjust = 1.17, fontface = "bold", size = 3.35) +
  ggplot2::facet_wrap(~domain, nrow = 1, scales = "free_x") +
  ggplot2::scale_color_manual(values = domain_colors) +
  ggplot2::scale_fill_manual(values = domain_colors) +
  ggplot2::labs(title = "GSE185263: domain associations with 24-hour SOFA",
                x = "Domain score (SD units)",
                y = "SOFA score at 24 hours") +
  ggprism::theme_prism(base_size = 13) +
  ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(face = "bold"),
                 plot.title = ggplot2::element_text(face = "bold", hjust = .5),
                 plot.caption = ggplot2::element_text(size = 9, hjust = 0))
save_plot(p_sofa, "GSE185263_RNA_LASSO_joint_24h_SOFA_quantile_curves")

# ----------------------------------------------------------------------------
# 5. Sex comparisons: unadjusted Wilcoxon + adjusted HC3-robust linear models.
#    Overall: age, disease state, recruitment site.
#    Sepsis only: age, recruitment site. Outcome = one domain score.
# ----------------------------------------------------------------------------
sex_comparisons <- purrr::map_dfr(method_key$method, function(m) {
  original <- method_data[[m]]
  purrr::map_dfr(c("Overall", "Sepsis only"), function(pop) {
    dd <- if (pop == "Overall") original else
      dplyr::filter(original, disease_state == "Sepsis")
    purrr::map_dfr(domains, function(d) {
      use <- dd[!is.na(dd$sex) & is.finite(dd[[d]]), , drop = FALSE]
      f <- use[[d]][use$sex == "F"]
      male <- use[[d]][use$sex == "M"]
      if (!length(f) || !length(male)) stop("Sex groups not present: ", m, "/", d)
      w <- stats::wilcox.test(f, male, exact = FALSE)
      adjust <- if (pop == "Overall")
        c("sex", "age", "disease_state", "RecruitmentSite") else
        c("sex", "age", "RecruitmentSite")
      dd_adj <- use[stats::complete.cases(use[, c(d, adjust), drop = FALSE]), , drop = FALSE]
      dd_adj$RecruitmentSite <- factor(dd_adj$RecruitmentSite)
      lmfit <- stats::lm(stats::reformulate(adjust, response = d), data = dd_adj)
      coef_hc3 <- lmtest::coeftest(lmfit,
        vcov. = sandwich::vcovHC(lmfit, type = "HC3"))
      if (!"sexM" %in% rownames(coef_hc3)) stop("Could not estimate male-vs-female effect")
      bet <- unname(coef_hc3["sexM", "Estimate"])
      se <- unname(coef_hc3["sexM", "Std. Error"])
      p_adj <- unname(coef_hc3["sexM", 4L])
      crit <- stats::qt(.975, df = stats::df.residual(lmfit))
      tibble::tibble(method = m, population = pop, domain = d,
        n_F = length(f), n_M = length(male),
        median_F = stats::median(f), median_M = stats::median(male),
        wilcoxon_W = unname(w$statistic), p_wilcoxon = w$p.value,
        n_adjusted = stats::nobs(lmfit),
        adjusted_beta_M_vs_F = bet, adjusted_lower95 = bet - crit * se,
        adjusted_upper95 = bet + crit * se, p_adjusted_HC3 = p_adj,
        adjusted_covariates = if (pop == "Overall")
          "Age, disease state, recruitment location" else "Age, recruitment location")
    })
  })
}) %>%
  dplyr::group_by(method, population) %>%
  dplyr::mutate(q_BH_Wilcoxon = stats::p.adjust(p_wilcoxon, "BH"),
                q_BH_adjusted = stats::p.adjust(p_adjusted_HC3, "BH")) %>%
  dplyr::ungroup()
save_csv(sex_comparisons, "GSE185263_DRT_sex_overall_and_sepsis_Wilcoxon_adjusted")

plot_sex <- function(pop) {
  dd <- if (pop == "Overall") rna else dplyr::filter(rna, disease_state == "Sepsis")
  pp <- dd %>%
    tidyr::pivot_longer(cols = dplyr::all_of(domains),
                        names_to = "domain", values_to = "score") %>%
    dplyr::filter(is.finite(score), !is.na(sex)) %>%
    dplyr::mutate(domain = factor(domain, levels = domains))
  labs <- sex_comparisons %>%
    dplyr::filter(method == "Transcriptomic projection", population == pop) %>%
    dplyr::mutate(domain = factor(domain, levels = domains),
                  label = paste0("Wilcoxon ", format_p(p_wilcoxon)))
  ggplot2::ggplot(pp,
    ggplot2::aes(x = sex, y = score, color = sex, fill = sex)) +
    ggplot2::geom_boxplot(width = .56, linewidth = .65,
                          outlier.shape = NA, alpha = .40, color = "black") +
    ggplot2::geom_jitter(width = .14, height = 0, alpha = .48, size = 1.5) +
    ggplot2::geom_text(data = labs, ggplot2::aes(x = 1.5, y = Inf, label = label),
                       inherit.aes = FALSE, vjust = 1.22,
                       fontface = "bold", size = 3.6) +
    ggplot2::facet_wrap(~domain, nrow = 1, scales = "fixed") +
    ggplot2::scale_color_manual(values = c(F = "#B56576", M = "#4A6FA5")) +
    ggplot2::scale_fill_manual(values = c(F = "#B56576", M = "#4A6FA5")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.06, .19))) +
    ggplot2::labs(
      title = if (pop == "Overall") "GSE185263: domains by sex (all participants)" else
        "GSE185263: domains by sex (sepsis only)",
      x = "Sex", y = "Transcriptomic projection domain score (sepsis SD units)") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(legend.position = "none",
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold", hjust = .5))
}
save_plot(plot_sex("Overall"), "GSE185263_RNA_LASSO_domain_scores_by_sex_overall")
save_plot(plot_sex("Sepsis only"), "GSE185263_RNA_LASSO_domain_scores_by_sex_sepsis")

# ----------------------------------------------------------------------------
# 6. Missingness / model sample audits (sepsis only, transcriptomic projection).
# ----------------------------------------------------------------------------
audit_table <- purrr::map_dfr(method_key$method, function(m) {
  dd <- method_data[[m]] %>% dplyr::filter(disease_state == "Sepsis")
  cc_common <- stats::complete.cases(dd[, model_predictors, drop = FALSE])
  tibble::tibble(method = m, n_total_sepsis = nrow(dd),
    n_with_all_three_domains = sum(stats::complete.cases(dd[, domains, drop = FALSE])),
    n_with_age = sum(!is.na(dd$age)), n_with_sex = sum(!is.na(dd$sex)),
    n_with_recruitment_location = sum(!is.na(dd$RecruitmentSite)),
    n_known_mortality = sum(!is.na(dd$mortality01)),
    n_inhospital_deaths = sum(dd$mortality01 == 1, na.rm = TRUE),
    n_with_SOFA24 = sum(!is.na(dd$sofa24)),
    n_complete_mortality_joint = sum(cc_common & !is.na(dd$mortality01)),
    n_complete_SOFA_joint = sum(cc_common & !is.na(dd$sofa24)))
})
save_csv(audit_table, "GSE185263_DRT_clinical_model_sample_size_audit")

cat("\n=== Baseline healthy/sepsis comparisons ===\n")
print(baseline_comparison)
cat("\n=== Joint domain in-hospital mortality estimates ===\n")
print(mortality_table)
cat("\n=== Joint domain median SOFA24 estimates ===\n")
print(sofa_table)
cat("\n=== Joint domain domain VIFs for mortality and SOFA ===\n")
print(joint_vif)
cat("\n=== Sex comparisons, overall and sepsis only ===\n")
print(sex_comparisons)
cat("\n=== Model sample audits ===\n")
print(audit_table)
cat("\nAll outputs saved to: ", out_dir, "\n", sep = "")



# =============================================================================
# 7. SUPPORTIVE ENROLLMENT-SETTING ANALYSES: EMERGENCY ROOM VS ICU
# =============================================================================
# This section is APPENDED to the original analysis; all preceding code, models,
# plots, score definitions, and output files are unchanged.
#
# Inclusion: sepsis participants with GEO 'collection site:ch1' exactly
# 'Emergency Room' or 'ICU'. Do NOT confuse this field with the geographic
# 'collection location:ch1' used in the main analyses.
#
# For each enrollment setting, fit:
#   (A) In-hospital death ~ Resistance + Tolerance + Damage + age + sex
#       + recruitment location where identifiable (ED only)
#   (B) Median 24-hour SOFA ~ same predictors, tau = 0.5
# Domain scores retain their ORIGINAL supplied per-SD scale, no re-scaling.
#
# CRITICAL INTERPRETATION: all ICU samples are from Toronto; ED samples were
# recruited elsewhere. Thus enrollment setting and recruiting geography are
# confounded. No formal ED-vs-ICU domain-by-setting interaction is estimated.
# ICU analyses are exploratory with approximately 20 mortality events; their
# confidence intervals and estimates may be unstable.
# =============================================================================

if (!"collection site:ch1" %in% names(meta)) {
  stop("ED/ICU analysis requires 'collection site:ch1' in GEO metadata")
}

edicu_key <- meta %>%
  dplyr::transmute(
    sample_id = trimws(as.character(geo_accession)),
    enrollment_setting_raw = trimws(as.character(.data[["collection site:ch1"]]))
  )
if (anyDuplicated(edicu_key$sample_id)) {
  stop("Repeated GEO accession in ED/ICU metadata")
}

# The score export may also carry the enrollment setting. Require agreement
# when it does, but continue to use the GEO metadata as authoritative.
if ("collection site:ch1" %in% names(scores)) {
  ix <- match(scores$sample_id, edicu_key$sample_id)
  score_setting <- trimws(as.character(scores[["collection site:ch1"]]))
  meta_setting <- edicu_key$enrollment_setting_raw[ix]
  conflict <- !is.na(ix) & !is.na(score_setting) & !is.na(meta_setting) &
    score_setting != meta_setting
  if (any(conflict)) stop("Metadata and score file disagree on collection site")
}

edicu_data <- method_data[["Transcriptomic projection"]] %>%
  dplyr::left_join(edicu_key, by = "sample_id") %>%
  dplyr::filter(disease_state == "Sepsis") %>%
  dplyr::mutate(
    enrollment_setting = factor(
      enrollment_setting_raw, levels = c("Emergency Room", "ICU"))
  )

if (anyNA(edicu_data$enrollment_setting)) {
  warning("Some sepsis samples have an unexpected enrollment setting; ",
          "they are excluded from the ED/ICU comparison")
}
edicu_data <- edicu_data %>%
  dplyr::filter(!is.na(enrollment_setting)) %>%
  droplevels()

# Sample-flow and geographic recruitment-site diagnostics.
cc_mortality_vars <- c(domains, "age", "sex", "RecruitmentSite", "mortality01")
cc_sofa_vars <- c(domains, "age", "sex", "RecruitmentSite", "sofa24")
edicu_data$cc_mortality <- stats::complete.cases(
  edicu_data[, cc_mortality_vars, drop = FALSE])
edicu_data$cc_sofa <- stats::complete.cases(
  edicu_data[, cc_sofa_vars, drop = FALSE])

edicu_sample_qc <- edicu_data %>%
  dplyr::group_by(enrollment_setting) %>%
  dplyr::summarise(
    n_sepsis = dplyr::n(),
    n_known_mortality = sum(!is.na(mortality01)),
    inhospital_deaths = sum(mortality01 == 1L, na.rm = TRUE),
    n_with_sofa24 = sum(is.finite(sofa24)),
    median_sofa24 = stats::median(sofa24, na.rm = TRUE),
    n_complete_mortality = sum(cc_mortality),
    n_complete_sofa = sum(cc_sofa),
    n_geographic_locations = dplyr::n_distinct(RecruitmentSite, na.rm = TRUE),
    .groups = "drop"
  )

edicu_site_qc <- edicu_data %>%
  dplyr::group_by(enrollment_setting, RecruitmentSite) %>%
  dplyr::summarise(
    n_sepsis = dplyr::n(),
    n_known_mortality = sum(!is.na(mortality01)),
    deaths = sum(mortality01 == 1L, na.rm = TRUE),
    n_with_sofa24 = sum(is.finite(sofa24)),
    .groups = "drop"
  )

save_csv(edicu_sample_qc, "GSE185263_DRT_ED_ICU_sample_size_QC")
save_csv(edicu_site_qc, "GSE185263_DRT_ED_ICU_recruitment_location_QC")

# One JOINT three-domain model in each subgroup for each outcome. The complete-
# case samples are calculated separately for mortality and SOFA.
# ED: age + sex + recruitment-location fixed effects.
# ICU: age + sex (recruitment location is constant and is not estimable).
fit_edicu_joint <- function(setting, response = c("mortality01", "sofa24")) {
  response <- match.arg(response)
  subgroup <- edicu_data[
    as.character(edicu_data$enrollment_setting) == setting, , drop = FALSE]

  use_location <- dplyr::n_distinct(subgroup$RecruitmentSite, na.rm = TRUE) > 1L
  regressors <- c(domains, "age", "sex",
                  if (use_location) "RecruitmentSite" else character(0))
  required <- unique(c(response, regressors))
  dd <- subgroup[stats::complete.cases(subgroup[, required, drop = FALSE]),
                 , drop = FALSE]
  dd <- droplevels(dd)

  # Detect geographic locations disappearing after outcome complete-case filtering.
  if (use_location && nlevels(factor(dd$RecruitmentSite)) < 2L) {
    warning(setting, ": only one recruitment location in complete cases; ",
            "fitting age/sex-only model")
    use_location <- FALSE
    regressors <- setdiff(regressors, "RecruitmentSite")
  }
  if (nrow(dd) < 25L || length(unique(dd$sex)) < 2L) {
    stop("Insufficient complete-case sample for ", setting, " / ", response)
  }
  if (response == "mortality01" && length(unique(dd$mortality01)) < 2L) {
    stop("Mortality outcome has only one class for ", setting)
  }

  if (use_location) {
    dd$RecruitmentSite <- droplevels(factor(dd$RecruitmentSite))
    # Some locations have zero deaths (e.g., Vancouver in the ED subset).
    # Warn about sparse strata; do not silently change the prespecified model.
    if (response == "mortality01") {
      events_by_site <- tapply(dd$mortality01, dd$RecruitmentSite,
                               function(x) sum(x == 1L))
      if (any(events_by_site == 0L)) {
        warning(setting, " includes a recruitment site with zero deaths. ",
                "Inspect logistic-model stability and the site QC table.")
      }
    }
  }

  form <- stats::reformulate(regressors, response = response)
  if (response == "mortality01") {
    fit <- stats::glm(form, data = dd, family = stats::binomial())
    if (!isTRUE(fit$converged)) warning(setting, ": mortality GLM did not converge")
    bb <- stats::coef(fit)[domains]
    se <- sqrt(diag(stats::vcov(fit)))[domains]
    pvalue <- 2 * stats::pnorm(abs(bb / se), lower.tail = FALSE)
    est <- exp(bb)
    lo <- exp(bb - 1.96 * se)
    hi <- exp(bb + 1.96 * se)
    outcome_label <- "In-hospital mortality"
    measure <- "OR"
  } else {
    fit <- quantreg::rq(form, data = dd, tau = 0.5)
    summ <- tryCatch(summary(fit, se = "nid"), error = function(e) {
      warning(setting, ": SOFA NID SE failed; using bootstrap: ",
              conditionMessage(e))
      summary(fit, se = "boot", R = 499L)
    })
    ctab <- as.matrix(summ$coefficients)
    bb <- ctab[domains, 1L]
    se <- ctab[domains, 2L]
    pvalue <- if (ncol(ctab) >= 4L) ctab[domains, 4L] else
      2 * stats::pnorm(abs(bb / se), lower.tail = FALSE)
    est <- bb
    lo <- bb - 1.96 * se
    hi <- bb + 1.96 * se
    outcome_label <- "24-hour SOFA (median)"
    measure <- "beta"
  }
  if (any(!is.finite(est)) || any(!is.finite(se)) || any(se <= 0)) {
    warning("Non-finite or non-estimable domain estimates: ", setting, " / ", response)
  }

  estimates <- tibble::tibble(
    enrollment_setting = setting,
    outcome = outcome_label,
    measure = measure,
    domain = domains,
    n = nrow(dd),
    deaths = if (response == "mortality01") sum(dd$mortality01 == 1L) else NA_integer_,
    estimate = as.numeric(est), SE = as.numeric(se),
    lower95 = as.numeric(lo), upper95 = as.numeric(hi),
    p = as.numeric(pvalue),
    adjustment = if (use_location) "Age, sex, recruitment location" else "Age, sex"
  ) %>%
    dplyr::mutate(q_BH_3domains = stats::p.adjust(p, method = "BH"))

  # Domain-specific VIFs conditional on all remaining regressors, computed on
  # exactly the same complete-case sample as the corresponding outcome model.
  vifs <- purrr::map_dfr(domains, function(domain) {
    aux <- stats::lm(stats::reformulate(setdiff(regressors, domain),
                                        response = domain), data = dd)
    r2 <- summary(aux)$r.squared
    tibble::tibble(enrollment_setting = setting, outcome = outcome_label,
                   domain = domain, n = stats::nobs(aux), VIF = 1 / (1 - r2))
  })
  list(estimates = estimates, vifs = vifs)
}

edicu_grid <- tidyr::expand_grid(
  enrollment_setting = c("Emergency Room", "ICU"),
  response = c("mortality01", "sofa24")
)
edicu_models <- purrr::pmap(
  edicu_grid,
  function(enrollment_setting, response) fit_edicu_joint(enrollment_setting, response)
)
edicu_subgroup_estimates <- dplyr::bind_rows(lapply(edicu_models, `[[`, "estimates"))
edicu_subgroup_vifs <- dplyr::bind_rows(lapply(edicu_models, `[[`, "vifs"))

# Include the unchanged published full-cohort model output for descriptive
# comparison. No extra full-cohort fit is performed here.
edicu_full_mortality <- mortality_table %>%
  dplyr::filter(method == "Transcriptomic projection") %>%
  dplyr::transmute(
    enrollment_setting = "All locations", outcome = "In-hospital mortality",
    measure = "OR", domain, n, deaths,
    estimate = OR, SE, lower95, upper95, p,
    adjustment = "Age, sex, recruitment location", q_BH_3domains
  )
edicu_full_sofa <- sofa_table %>%
  dplyr::filter(method == "Transcriptomic projection") %>%
  dplyr::transmute(
    enrollment_setting = "All locations", outcome = "24-hour SOFA (median)",
    measure = "beta", domain, n, deaths = NA_integer_,
    estimate = beta, SE, lower95, upper95, p,
    adjustment = "Age, sex, recruitment location", q_BH_3domains
  )

edicu_comparison <- dplyr::bind_rows(
  edicu_full_mortality, edicu_full_sofa, edicu_subgroup_estimates) %>%
  dplyr::mutate(
    estimate_95CI = sprintf("%.2f (%.2f to %.2f)",
                            estimate, lower95, upper95)
  )

save_csv(edicu_comparison, "GSE185263_DRT_Full_ED_ICU_joint_mortality_SOFA_estimates")
save_csv(edicu_subgroup_vifs, "GSE185263_DRT_ED_ICU_joint_domain_VIFs")

# Side-by-side forests for All locations, ED only, and ICU only. Comparisons are
# DESCRIPTIVE; there is no formal setting interaction because enrollment
# setting is confounded with geographic recruitment location in these data.
plot_edicu_comparison <- function(outcome_name, measure = c("OR", "beta")) {
  measure <- match.arg(measure)
  dd <- edicu_comparison %>%
    dplyr::filter(outcome == outcome_name,
                  is.finite(estimate), is.finite(lower95), is.finite(upper95)) %>%
    dplyr::mutate(
      domain = factor(domain, levels = rev(domains)),
      enrollment_setting = factor(enrollment_setting,
                                  levels = c("All locations", "Emergency Room", "ICU"))
    )
  dodge <- ggplot2::position_dodge(width = 0.58)
  p <- ggplot2::ggplot(dd,
                       ggplot2::aes(x = estimate, y = domain,
                                    color = enrollment_setting,
                                    shape = enrollment_setting)) +
    ggplot2::geom_vline(xintercept = if (measure == "OR") 1 else 0,
                        linetype = "dashed", color = "grey55") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower95, xmax = upper95),
      orientation = "y", position = dodge, width = 0.17, linewidth = 0.85) +
    ggplot2::geom_point(position = dodge, size = 3.4) +
    ggplot2::scale_color_manual(values = c(
      "All locations" = "#646464", "Emergency Room" = "#2A6EBB", "ICU" = "#D65F44")) +
    ggplot2::scale_shape_manual(values = c(
      "All locations" = 16, "Emergency Room" = 17, "ICU" = 15)) +
    ggplot2::labs(
      title = paste0("GSE185263: domain and ", outcome_name, " by enrollment setting"),
      x = if (measure == "OR")
        "Adjusted mortality OR per 1-SD domain score" else
        "Adjusted median-SOFA difference per 1-SD domain score",
      y = NULL, color = "Enrollment setting", shape = "Enrollment setting") +
    ggprism::theme_prism(base_size = 12) +
    ggplot2::theme(
      legend.position = "top", axis.text.y = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.caption = ggplot2::element_text(size = 9, hjust = 0)
    )
  if (measure == "OR") p <- p + ggplot2::scale_x_log10()
  p
}

save_plot(
  plot_edicu_comparison("In-hospital mortality", "OR"),
  "GSE185263_DRT_Full_ED_ICU_joint_mortality_forest")
save_plot(
  plot_edicu_comparison("24-hour SOFA (median)", "beta"),
  "GSE185263_DRT_Full_ED_ICU_joint_SOFA_forest")

cat("\n=== ED / ICU: samples and events ===\n")
print(edicu_sample_qc)
cat("\n=== ED / ICU: joint R/T/D mortality and 24-hour SOFA ===\n")
print(edicu_comparison)
cat("\n=== ED / ICU: domain VIFs ===\n")
print(edicu_subgroup_vifs)
cat("\nED/ICU comparison outputs saved to: ", out_dir, "\n", sep = "")
# =============================================================================
# END of appended ED / ICU comparison; original full-cohort code is unchanged.
# =============================================================================
