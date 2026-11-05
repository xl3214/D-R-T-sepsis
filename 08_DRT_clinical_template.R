# =============================================================================
# Damage / Resistance / Tolerance (D/R/T) domain scores: clinical analysis template
# Script: 08_DRT_clinical_template.R
# =============================================================================
#
# PURPOSE
#   A configurable, standalone version of the clinical analyses in the
#   manuscript, for use in any cohort scored with 01_DRT_scores.R. The
#   cohort-specific scripts in clinical/ reproduce the published results; this
#   template applies the same analytic framework to new data.
#
# MODULES (switch on/off in Section 1)
#   A. Case vs control        Wilcoxon rank-sum test of each domain score
#                             between cases and controls (BH across domains).
#   B. Joint outcome models   D/R/T scores entered together (mutually adjusted)
#                             with user-chosen covariates:
#                               binary  -> logistic regression (OR per SD), with
#                                          Firth penalized logistic sensitivity
#                               ordinal -> proportional-odds model (OR per SD)
#                               median  -> median (tau = 0.5) quantile regression
#                               linear  -> linear regression
#                             VIFs for the three domain scores in each model.
#   C. Exposure -> score      Linear regression of each domain score on each
#                             binary exposure (e.g., an infection), adjusted for
#                             covariates and, optionally, the other exposures;
#                             BH across exposures within each domain. Optional
#                             exploratory exposure x domain interactions in the
#                             joint binary-outcome model (likelihood-ratio test,
#                             BH across all exposure x domain models), with
#                             domain ORs within each exposure stratum.
#   D. Longitudinal           Linear mixed-effects model per domain with a
#                             time x outcome-group interaction (random subject
#                             intercept; LRT for the interaction), time-specific
#                             group contrasts, and a landmark analysis of score
#                             changes from baseline (Wilcoxon; joint logistic
#                             model of changes + baseline scores + covariates,
#                             Firth sensitivity).
#
# INPUT
#   SCORE_FILE     <cohort>_DRT_scores_wide_perSD.csv from 01_DRT_scores.R
#                  (one row per sample; columns <method>__Damage, etc.).
#   CLINICAL_FILE  one row per sample (longitudinal: one row per sample, with a
#                  subject ID and time column; subject-level variables repeated).
#   Binary variables (outcomes, exposures) may be coded 0/1, TRUE/FALSE, or as
#   text; set the level that denotes the event / exposure in Section 1.
#
# SCALE
#   Scores are analysed in the SD units written by 01_DRT_scores.R (mean 0, SD 1
#   in the cohort's reference samples), so estimates are per 1 SD.
#
# REQUIRED PACKAGES
#   dplyr, tidyr, ggplot2; MASS (ordinal); quantreg (median); lme4 and emmeans
#   (longitudinal); logistf (optional, Firth sensitivity).
# =============================================================================


# =============================================================================
# SECTION 1. CONFIGURATION -- edit here
# =============================================================================
SCORE_FILE    <- "DRT_scores_output/MyCohort/MyCohort_DRT_scores_wide_perSD.csv"
CLINICAL_FILE <- "data/MyCohort_clinical.csv"
OUT_DIR       <- "DRT_scores_output/Clinical analyses/MyCohort"
COHORT_LABEL  <- "MyCohort"

SCORE_ID_COL    <- "sample_id"     # ID column in SCORE_FILE
CLINICAL_ID_COL <- "sample_id"     # matching ID column in CLINICAL_FILE
# Score method: "rna_surrogate_lasso14_1se" (transcriptomic projection) or
# "multiomic_PC1" (cohorts with Olink + RNA-seq).
SCORE_METHOD <- "rna_surrogate_lasso14_1se"

COVARIATES <- c("age", "sex")      # adjustment covariates for modules B-D

# Case/control definition. Set CASE_COL <- NULL if all samples are cases.
CASE_COL      <- "disease_state"
CASE_LEVEL    <- "sepsis"
CONTROL_LEVEL <- "healthy"

# Longitudinal design. Set TIME_COL <- NULL for a single time point.
# Modules A-C use baseline samples only.
SUBJECT_COL   <- NULL              # e.g. "subject_id"
TIME_COL      <- NULL              # e.g. "time"
TIME_LEVELS   <- NULL              # ordered time points, e.g. c("0h", "24h", "72h"); first = baseline

RUN_A_CASE_CONTROL  <- TRUE
RUN_B_OUTCOMES      <- TRUE
RUN_C_EXPOSURES     <- FALSE
RUN_D_LONGITUDINAL  <- FALSE

# B. Outcomes: one entry per outcome. type = "binary" | "ordinal" | "median" | "linear".
#    event = value of `col` denoting the event; nonevent (optional) = value(s)
#    denoting no event, other values are then treated as missing (binary only).
OUTCOMES <- list(
  list(name = "Mortality", col = "death", type = "binary", event = "1", nonevent = "0"),
  list(name = "Severity score", col = "severity_score", type = "ordinal")
)
RUN_FIRTH <- TRUE                  # Firth sensitivity for binary outcomes (needs logistf)

# C. Exposures (binary). exposed = value denoting exposure.
EXPOSURES <- list(
  list(name = "HIV", col = "hiv", exposed = "1"),
  list(name = "Malaria", col = "malaria", exposed = "1")
)
EXPOSURE_ADJUST_OTHER   <- TRUE    # adjust each exposure model for the other exposures
RUN_INTERACTIONS        <- TRUE    # exploratory exposure x domain interactions
INTERACTION_OUTCOME     <- "Mortality"   # name of a binary outcome in OUTCOMES

# D. Longitudinal: binary subject-level group (usually the outcome) and landmark time.
LONG_GROUP_COL      <- "death"
LONG_GROUP_EVENT    <- "1"
LONG_GROUP_NONEVENT <- "0"         # NULL: any other non-missing value
LANDMARK_TIME    <- NULL           # e.g. "72h"; NULL to skip the landmark analysis

SEED <- 20261002


# =============================================================================
# SECTION 2. SETUP
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

DOMAINS <- c("Damage", "Resistance", "Tolerance")
DCOL    <- c(Damage = "#E6550D", Resistance = "#2A6EBB", Tolerance = "#2CA25F")
types   <- vapply(OUTCOMES, `[[`, character(1), "type")
need <- c(if (RUN_B_OUTCOMES && "ordinal" %in% types) "MASS",
          if (RUN_B_OUTCOMES && "median" %in% types) "quantreg",
          if (RUN_D_LONGITUDINAL) c("lme4", "emmeans"))
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Please install: ", paste(miss, collapse = ", "))
HAS_FIRTH <- RUN_FIRTH && requireNamespace("logistf", quietly = TRUE)
if (RUN_FIRTH && !HAS_FIRTH) warning("Package 'logistf' not installed; Firth sensitivity skipped.")
if (RUN_D_LONGITUDINAL && (is.null(SUBJECT_COL) || is.null(TIME_COL) || is.null(TIME_LEVELS)))
  stop("Module D requires SUBJECT_COL, TIME_COL and TIME_LEVELS.")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out_csv <- function(x, nm) utils::write.csv(x, file.path(OUT_DIR, paste0(COHORT_LABEL, "_", nm, ".csv")), row.names = FALSE)
out_pdf <- function(p, nm, w = 7, h = 4.5) ggsave(file.path(OUT_DIR, paste0(COHORT_LABEL, "_", nm, ".pdf")), p, width = w, height = h)

# Binary coding: 1 if value equals `level` (compared as text); 0 if it equals
# `nonlevel` (or, when nonlevel = NULL, any other non-missing value); NA otherwise
to01 <- function(x, level, nonlevel = NULL) {
  x <- trimws(as.character(x))
  out <- ifelse(is.na(x) | x == "", NA_real_, as.numeric(x == as.character(level)))
  if (!is.null(nonlevel)) out[!is.na(out) & out == 0 & !(x %in% as.character(nonlevel))] <- NA_real_
  out
}


# =============================================================================
# SECTION 3. DATA
# =============================================================================
scores <- utils::read.csv(SCORE_FILE, check.names = FALSE, stringsAsFactors = FALSE)
clin   <- utils::read.csv(CLINICAL_FILE, check.names = FALSE, stringsAsFactors = FALSE)
score_cols <- paste0(SCORE_METHOD, "__", DOMAINS)
absent <- setdiff(c(SCORE_ID_COL, score_cols), names(scores))
if (length(absent)) stop("Columns not found in SCORE_FILE: ", paste(absent, collapse = ", "))

dat <- scores[, c(SCORE_ID_COL, score_cols)]
names(dat) <- c(".id", DOMAINS)
clin[[".id"]] <- trimws(as.character(clin[[CLINICAL_ID_COL]]))
dat$.id <- trimws(as.character(dat$.id))
clin <- clin[, setdiff(names(clin), DOMAINS), drop = FALSE]
dat <- inner_join(dat, clin[!duplicated(clin$.id), , drop = FALSE], by = ".id")
message(sprintf("%d samples with scores and clinical data", nrow(dat)))

used <- unique(c(COVARIATES, CASE_COL, SUBJECT_COL, TIME_COL,
                 if (RUN_B_OUTCOMES) vapply(OUTCOMES, `[[`, character(1), "col"),
                 if (RUN_C_EXPOSURES) vapply(EXPOSURES, `[[`, character(1), "col"),
                 if (RUN_D_LONGITUDINAL) LONG_GROUP_COL))
absent <- setdiff(used, names(dat))
if (length(absent)) stop("Columns not found in CLINICAL_FILE: ", paste(absent, collapse = ", "))
for (v in COVARIATES) if (is.character(dat[[v]])) dat[[v]] <- factor(dat[[v]])

if (!is.null(TIME_COL)) dat[[".time"]] <- factor(trimws(as.character(dat[[TIME_COL]])), levels = TIME_LEVELS)
is_case <- if (is.null(CASE_COL)) rep(TRUE, nrow(dat)) else trimws(as.character(dat[[CASE_COL]])) == CASE_LEVEL
is_base <- if (is.null(TIME_COL)) rep(TRUE, nrow(dat)) else !is.na(dat$.time) & dat$.time == TIME_LEVELS[1]
base_cases <- dat[is_case & is_base, , drop = FALSE]      # analysis sample for modules B and C
message(sprintf("%d baseline case samples", nrow(base_cases)))


# =============================================================================
# SECTION 4. HELPERS
# =============================================================================
# Wald estimates for the named terms of a fitted model
tidy_terms <- function(fit, terms, model_type) {
  if (model_type == "median") {
    # NID standard errors (bootstrap if NID fails), as in the GSE185263 analysis
    sm <- tryCatch(summary(fit, se = "nid"), error = function(e) { set.seed(SEED); summary(fit, se = "boot", R = 499L) })
    cf <- as.matrix(sm$coefficients)
    est <- cf[terms, 1]; se <- cf[terms, 2]; p <- cf[terms, 4]
  } else {
    est <- stats::coef(fit)[terms]; se <- sqrt(diag(stats::vcov(fit)))[terms]
    z <- est / se
    p <- if (model_type == "linear") 2 * stats::pt(-abs(z), df = stats::df.residual(fit)) else 2 * stats::pnorm(-abs(z))
  }
  ratio <- model_type %in% c("binary", "ordinal")
  tf <- if (ratio) exp else identity
  tibble(term = terms, estimate = tf(est), conf_low = tf(est - 1.96 * se), conf_high = tf(est + 1.96 * se),
         p = unname(p), scale = ifelse(ratio, "OR per SD", "difference per SD"))
}

# VIF of each domain score from auxiliary regressions on the other predictors
domain_vif <- function(d, rhs) {
  vapply(DOMAINS, function(x) {
    r2 <- summary(stats::lm(stats::reformulate(setdiff(rhs, x), response = x), data = d))$r.squared
    1 / (1 - r2)
  }, numeric(1))
}

fit_model <- function(form, d, type) {
  switch(type,
    binary  = stats::glm(form, data = d, family = stats::binomial()),
    ordinal = MASS::polr(form, data = d, Hess = TRUE, method = "logistic"),
    median  = quantreg::rq(form, data = d, tau = 0.5),
    linear  = stats::lm(form, data = d))
}

forest <- function(tab, xlab, log_x, title) {
  tab <- tab %>% mutate(term = factor(term, levels = rev(unique(term))))
  p <- ggplot(tab, aes(estimate, term, colour = domain)) +
    geom_vline(xintercept = if (log_x) 1 else 0, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.2, orientation = "y") +
    geom_point(size = 2.2) +
    scale_colour_manual(values = DCOL, guide = "none") +
    labs(x = xlab, y = NULL, title = title) + theme_bw(base_size = 10)
  if (log_x) p + scale_x_log10() else p
}


# =============================================================================
# SECTION 5. MODULE A -- CASE VS CONTROL
# =============================================================================
if (RUN_A_CASE_CONTROL && !is.null(CASE_COL)) {
  message("\n=== A. Case vs control ===")
  dA <- dat[is_base & trimws(as.character(dat[[CASE_COL]])) %in% c(CASE_LEVEL, CONTROL_LEVEL), , drop = FALSE]
  dA$.group <- factor(trimws(as.character(dA[[CASE_COL]])), levels = c(CONTROL_LEVEL, CASE_LEVEL))
  resA <- bind_rows(lapply(DOMAINS, function(dm) {
    x <- dA[[dm]][dA$.group == CASE_LEVEL]; y <- dA[[dm]][dA$.group == CONTROL_LEVEL]
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    tibble(domain = dm, n_case = length(x), n_control = length(y),
           median_case = stats::median(x), median_control = stats::median(y),
           difference_in_medians = stats::median(x) - stats::median(y),
           p = stats::wilcox.test(x, y, exact = FALSE)$p.value)
  })) %>% mutate(q_BH = stats::p.adjust(p, method = "BH"))
  out_csv(resA, "A_case_vs_control_Wilcoxon")
  print(as.data.frame(resA))
  pA <- dA %>% pivot_longer(all_of(DOMAINS), names_to = "domain", values_to = "score") %>%
    filter(is.finite(score)) %>%
    ggplot(aes(.group, score, fill = domain)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_point(position = position_jitter(width = 0.15, height = 0, seed = SEED), size = 0.6, alpha = 0.5) +
    facet_wrap(~ domain) + scale_fill_manual(values = DCOL, guide = "none") +
    labs(x = NULL, y = "Domain score (SD units)") + theme_bw(base_size = 10)
  out_pdf(pA, "A_case_vs_control_boxplots", 8, 3.5)
}


# =============================================================================
# SECTION 6. MODULE B -- JOINT OUTCOME MODELS
# =============================================================================
if (RUN_B_OUTCOMES) {
  message("\n=== B. Joint outcome models ===")
  rhs <- c(DOMAINS, COVARIATES)
  resB <- list(); vifB <- list(); firthB <- list()
  for (o in OUTCOMES) {
    d <- base_cases
    yv <- d[[o$col]]
    num <- suppressWarnings(as.numeric(as.character(yv)))
    if (o$type == "binary") {
      d$.y <- to01(yv, o$event, o$nonevent)
    } else if (o$type == "ordinal") {
      # numeric scores (e.g. MEWS) are ordered numerically; text categories as given
      d$.y <- if (all(is.na(yv) | !is.na(num))) factor(num, levels = sort(unique(num[!is.na(num)])), ordered = TRUE)
              else factor(yv, ordered = TRUE)
    } else {
      d$.y <- num
    }
    d <- d[stats::complete.cases(d[, c(".y", rhs)]), , drop = FALSE]
    if (o$type == "ordinal") d$.y <- droplevels(d$.y)
    form <- stats::reformulate(rhs, response = ".y")
    fit <- fit_model(form, d, o$type)
    n_ev <- if (o$type == "binary") sum(d$.y == 1) else NA_integer_
    resB[[o$name]] <- tidy_terms(fit, DOMAINS, o$type) %>%
      mutate(outcome = o$name, model = o$type, n = nrow(d), n_events = n_ev, domain = term, .before = 1)
    vifB[[o$name]] <- tibble(outcome = o$name, domain = DOMAINS, VIF = domain_vif(d, rhs))
    if (o$type == "binary" && HAS_FIRTH) {
      ff <- tryCatch(logistf::logistf(form, data = d, firth = TRUE, pl = TRUE), error = function(e) NULL)
      if (!is.null(ff)) {
        ci <- stats::confint(ff)
        firthB[[o$name]] <- tibble(outcome = o$name, domain = DOMAINS, OR = exp(stats::coef(ff)[DOMAINS]),
                                   conf_low = exp(ci[DOMAINS, 1]), conf_high = exp(ci[DOMAINS, 2]),
                                   p = unname(ff$prob[DOMAINS]), n = nrow(d), n_events = n_ev)
      }
    }
    message(sprintf("  %s (%s): n = %d%s", o$name, o$type, nrow(d), if (is.na(n_ev)) "" else paste0(", events = ", n_ev)))
  }
  resB <- bind_rows(resB); out_csv(resB, "B_joint_outcome_models"); print(as.data.frame(resB))
  out_csv(bind_rows(vifB), "B_joint_outcome_models_VIF")
  if (length(firthB)) out_csv(bind_rows(firthB), "B_joint_binary_outcomes_Firth_sensitivity")
  for (o in OUTCOMES) {
    tb <- resB %>% filter(outcome == o$name) %>% mutate(term = domain)
    ratio <- o$type %in% c("binary", "ordinal")
    out_pdf(forest(tb, if (ratio) "Adjusted OR per SD (95% CI)" else "Adjusted difference per SD (95% CI)", ratio,
                   paste0(o$name, ": mutually adjusted D/R/T")), paste0("B_forest_", gsub("[^A-Za-z0-9]+", "_", o$name)), 6, 2.8)
  }
}


# =============================================================================
# SECTION 7. MODULE C -- EXPOSURE -> SCORE, AND EXPLORATORY INTERACTIONS
# =============================================================================
if (RUN_C_EXPOSURES) {
  message("\n=== C. Exposure -> score ===")
  d0 <- base_cases
  ex_names <- vapply(EXPOSURES, `[[`, character(1), "name")
  ex_vars  <- paste0(".ex_", make.names(ex_names))
  for (i in seq_along(EXPOSURES)) d0[[ex_vars[i]]] <- to01(d0[[EXPOSURES[[i]]$col]], EXPOSURES[[i]]$exposed)

  resC <- bind_rows(lapply(DOMAINS, function(dm) bind_rows(lapply(seq_along(EXPOSURES), function(i) {
    rhs <- c(ex_vars[i], if (EXPOSURE_ADJUST_OTHER) ex_vars[-i], COVARIATES)
    d <- d0[stats::complete.cases(d0[, c(dm, rhs)]), , drop = FALSE]
    fit <- stats::lm(stats::reformulate(rhs, response = dm), data = d)
    tidy_terms(fit, ex_vars[i], "linear") %>%
      mutate(domain = dm, exposure = ex_names[i], n = nrow(d), n_exposed = sum(d[[ex_vars[i]]] == 1), .before = 1)
  })) %>% mutate(q_BH_within_domain = stats::p.adjust(p, method = "BH"))))
  out_csv(resC %>% dplyr::select(-term), "C_exposure_to_score_linear")
  print(as.data.frame(resC))
  pC <- forest(resC %>% mutate(term = paste(exposure, domain, sep = ": ")),
               "Adjusted difference in score, SD (95% CI)", FALSE, "Exposure -> domain score")
  out_pdf(pC, "C_exposure_to_score_forest", 6, 1 + 0.35 * nrow(resC))

  if (RUN_INTERACTIONS) {
    o <- Filter(function(x) x$name == INTERACTION_OUTCOME, OUTCOMES)
    if (!length(o) || o[[1]]$type != "binary") stop("INTERACTION_OUTCOME must name a binary outcome in OUTCOMES.")
    o <- o[[1]]
    resI <- list(); subI <- list()
    for (i in seq_along(EXPOSURES)) for (dm in DOMAINS) {
      d <- d0; d$.y <- to01(d[[o$col]], o$event, o$nonevent)
      rhs <- c(DOMAINS, ex_vars[i], COVARIATES)
      d <- d[stats::complete.cases(d[, c(".y", rhs)]), , drop = FALSE]
      int <- paste0(dm, ":", ex_vars[i])
      red  <- stats::glm(stats::reformulate(rhs, response = ".y"), data = d, family = stats::binomial())
      full <- stats::glm(stats::reformulate(c(rhs, int), response = ".y"), data = d, family = stats::binomial())
      lrt  <- stats::anova(red, full, test = "LRT")
      b <- stats::coef(full); V <- stats::vcov(full)
      # stratum-specific domain ORs from the interaction model
      for (lv in 0:1) {
        g <- setNames(rep(0, length(b)), names(b)); g[dm] <- 1; if (lv == 1) g[int] <- 1
        est <- sum(g * b); se <- sqrt(as.numeric(t(g) %*% V %*% g))
        subI[[paste(i, dm, lv)]] <- tibble(exposure = ex_names[i], domain = dm,
          stratum = ifelse(lv == 1, "exposed", "unexposed"), n = sum(d[[ex_vars[i]]] == lv),
          n_events = sum(d$.y[d[[ex_vars[i]]] == lv] == 1),
          OR = exp(est), conf_low = exp(est - 1.96 * se), conf_high = exp(est + 1.96 * se))
      }
      resI[[paste(i, dm)]] <- tibble(outcome = o$name, exposure = ex_names[i], domain = dm, n = nrow(d),
        n_events = sum(d$.y == 1), interaction_ratio_of_ORs = exp(b[int]), p_LRT = lrt[2, "Pr(>Chi)"])
    }
    resI <- bind_rows(resI) %>% mutate(q_BH_all_models = stats::p.adjust(p_LRT, method = "BH"))
    out_csv(resI, "C_exploratory_exposure_by_domain_interactions")
    out_csv(bind_rows(subI), "C_exploratory_interaction_stratum_ORs")
    print(as.data.frame(resI))
  }
}


# =============================================================================
# SECTION 8. MODULE D -- LONGITUDINAL
# =============================================================================
if (RUN_D_LONGITUDINAL) {
  message("\n=== D. Longitudinal ===")
  dL <- dat[is_case & !is.na(dat$.time), , drop = FALSE]
  dL$.subject <- factor(dL[[SUBJECT_COL]])
  dL$.group <- factor(ifelse(to01(dL[[LONG_GROUP_COL]], LONG_GROUP_EVENT, LONG_GROUP_NONEVENT) == 1, "Event", "No event"),
                      levels = c("No event", "Event"))
  # subjects with a baseline sample and a known group
  dL <- dL[!is.na(dL$.group) & dL$.subject %in% dL$.subject[dL$.time == TIME_LEVELS[1]], , drop = FALSE]
  dL$.subject <- droplevels(dL$.subject)

  testsD <- list(); contrD <- list()
  for (dm in DOMAINS) {
    d <- dL[stats::complete.cases(dL[, c(dm, ".time", ".group", ".subject", COVARIATES)]), , drop = FALSE]
    d$.score <- d[[dm]]
    cov_txt <- if (length(COVARIATES)) paste("+", paste(COVARIATES, collapse = " + ")) else ""
    f_full <- stats::as.formula(paste(".score ~ .time * .group", cov_txt, "+ (1 | .subject)"))
    f_red  <- stats::as.formula(paste(".score ~ .time + .group", cov_txt, "+ (1 | .subject)"))
    ctl <- lme4::lmerControl(optimizer = "bobyqa")
    full <- lme4::lmer(f_full, data = d, REML = FALSE, control = ctl)   # ML for the LRT
    red  <- lme4::lmer(f_red,  data = d, REML = FALSE, control = ctl)
    lrt  <- stats::anova(red, full)
    testsD[[dm]] <- tibble(domain = dm, n_samples = nrow(d), n_subjects = n_distinct(d$.subject),
                           chisq = lrt[2, "Chisq"], df = lrt[2, "Df"], p_time_by_group = lrt[2, "Pr(>Chisq)"])
    emm <- emmeans::emmeans(full, ~ .group | .time, lmer.df = "asymptotic")
    cc  <- as.data.frame(summary(emmeans::contrast(emm, method = list("Event - No event" = c(-1, 1))),
                                 infer = c(TRUE, TRUE), adjust = "none"))
    lo <- grep("(lower.CL|asymp.LCL)$", names(cc), value = TRUE)[1]; hi <- grep("(upper.CL|asymp.UCL)$", names(cc), value = TRUE)[1]
    contrD[[dm]] <- tibble(domain = dm, time = as.character(cc$.time), difference = cc$estimate,
                           conf_low = cc[[lo]], conf_high = cc[[hi]], p = cc$p.value)
  }
  out_csv(bind_rows(testsD), "D_longitudinal_LME_time_by_group_LRT")
  out_csv(bind_rows(contrD), "D_longitudinal_time_specific_Event_minus_NoEvent")
  print(as.data.frame(bind_rows(testsD)))
  pD <- dL %>% pivot_longer(all_of(DOMAINS), names_to = "domain", values_to = "score") %>%
    filter(is.finite(score)) %>%
    ggplot(aes(.time, score, colour = .group, group = .group)) +
    stat_summary(fun = mean, geom = "line", position = position_dodge(0.2)) +
    stat_summary(fun.data = function(x) { m <- mean(x); h <- stats::qt(0.975, length(x) - 1) * stats::sd(x) / sqrt(length(x))
                                          data.frame(y = m, ymin = m - h, ymax = m + h) },
                 geom = "pointrange", position = position_dodge(0.2), size = 0.3) +
    facet_wrap(~ domain) + scale_colour_manual(values = c("No event" = "grey40", "Event" = "#B2182B"), name = NULL) +
    labs(x = NULL, y = "Domain score (SD units), mean (95% CI)") + theme_bw(base_size = 10)
  out_pdf(pD, "D_longitudinal_trajectories", 8, 3.5)

  if (!is.null(LANDMARK_TIME)) {
    w <- dL %>% filter(.time %in% c(TIME_LEVELS[1], LANDMARK_TIME)) %>%
      mutate(.tp = ifelse(.time == TIME_LEVELS[1], "base", "land")) %>%
      dplyr::select(.subject, .group, all_of(COVARIATES), .tp, all_of(DOMAINS)) %>%
      pivot_longer(all_of(DOMAINS), names_to = "domain", values_to = "score") %>%
      pivot_wider(names_from = c(.tp, domain), values_from = score, values_fn = mean)
    for (dm in DOMAINS) w[[paste0("delta_", dm)]] <- w[[paste0("land_", dm)]] - w[[paste0("base_", dm)]]
    w <- w[stats::complete.cases(w[, c(paste0("delta_", DOMAINS), paste0("base_", DOMAINS), COVARIATES)]), , drop = FALSE]
    wil <- bind_rows(lapply(DOMAINS, function(dm) {
      x <- w[[paste0("delta_", dm)]]
      tibble(domain = dm, n_event = sum(w$.group == "Event"), n_no_event = sum(w$.group == "No event"),
             median_change_event = stats::median(x[w$.group == "Event"]),
             median_change_no_event = stats::median(x[w$.group == "No event"]),
             p = stats::wilcox.test(x[w$.group == "Event"], x[w$.group == "No event"], exact = FALSE)$p.value)
    }))
    out_csv(wil, paste0("D_landmark_", LANDMARK_TIME, "_changes_Wilcoxon"))
    w$.y <- as.numeric(w$.group == "Event")
    rhs <- c(paste0("delta_", DOMAINS), paste0("base_", DOMAINS), COVARIATES)
    form <- stats::reformulate(rhs, response = ".y")
    fit <- stats::glm(form, data = w, family = stats::binomial())
    lm_tab <- tidy_terms(fit, paste0("delta_", DOMAINS), "binary") %>%
      mutate(domain = sub("^delta_", "", term), n = nrow(w), n_events = sum(w$.y), .before = 1)
    out_csv(lm_tab, paste0("D_landmark_", LANDMARK_TIME, "_joint_change_logistic"))
    if (HAS_FIRTH) {
      ff <- tryCatch(logistf::logistf(form, data = w, firth = TRUE, pl = TRUE), error = function(e) NULL)
      if (!is.null(ff)) {
        tt <- paste0("delta_", DOMAINS); ci <- stats::confint(ff)
        out_csv(tibble(domain = DOMAINS, OR = exp(stats::coef(ff)[tt]), conf_low = exp(ci[tt, 1]),
                       conf_high = exp(ci[tt, 2]), p = unname(ff$prob[tt])),
                paste0("D_landmark_", LANDMARK_TIME, "_joint_change_Firth_sensitivity"))
      }
    }
    print(as.data.frame(lm_tab))
  }
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, paste0(COHORT_LABEL, "_sessionInfo.txt")))
message("\nDone. Outputs in ", OUT_DIR)
