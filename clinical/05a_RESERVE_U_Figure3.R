# =============================================================================
# RESERVE-U: combined journal Figure 3 (three panels, 180 mm wide)
# Script: 05a_RESERVE_U_Figure3.R
# =============================================================================
#
# PANELS
#   A  Joint domain associations with MEWS and 30-day mortality
#      (cohort-specific + pooled aOR, log axis, aligned estimate/p table)
#   B  Adjusted HIV / malaria / MTB differences in domain scores
#      (shared beta axis across pathogens, aligned beta (95% CI) / p tables)
#   C  Domain associations with 30-day mortality by pathogen status
#      (shared log OR axis; filled = positive, open = negative;
#       1-df LRT interaction p-values)
#
# USE
#   source() this file at the END of 05_RESERVE_U_Pooled_DRT_Clinical.R
#   (or run it on its own after that script, from the project root).
#   If the model objects are not in memory, the estimates are read from the
#   CSVs that script writes to out_dir, so the models need not be refit.
#
# REQUIRED PACKAGES
#   ggplot2 (>= 3.4), dplyr, patchwork
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

if (!exists("out_dir")) {
  out_dir <- "DRT_scores_output/Clinical analyses/Pooled"
}
read_out <- function(file) {
  utils::read.csv(file.path(out_dir, file), check.names = FALSE,
                  stringsAsFactors = FALSE)
}
get_or_read <- function(obj, file) {
  if (exists(obj, inherits = TRUE)) get(obj, inherits = TRUE) else read_out(file)
}

fig_mews <- get_or_read("mews_estimates",
                        "RESERVE_U_POOLED_supportive_joint_DRT_MEWS_estimates.csv")
fig_mort <- get_or_read("pooled_estimates",
                        "RESERVE_U_POOLED_supportive_joint_DRT_30d_estimates.csv")
fig_inf  <- get_or_read("infection_associations",
                        "RESERVE_U_POOLED_DRT_adjusted_infection_type_associations.csv")
fig_sub  <- get_or_read("hte_subgroups",
                        "RESERVE_U_POOLED_DRT_infection_effect_modification_subgroup_ORs.csv")
fig_int  <- get_or_read("hte_interactions",
                        "RESERVE_U_POOLED_DRT_infection_effect_modification_30d_9_single_interactions.csv")

# ---------------- shared style ----------------
f3_domains <- c("Resistance", "Tolerance", "Damage")
f3_sources <- c("RESERVE-U-2-TOR", "RESERVE-U-1-EBB", "Pooled (cohort-adjusted)")
f3_source_labels <- c("RESERVE-U-2-TOR", "RESERVE-U-1-EBB", "Pooled")
f3_pathogens <- c("HIV", "Malaria", "MTB")
# Okabe-Ito: colorblind safe, and distinct from any status coding
f3_col <- c(Resistance = "#0072B2", Tolerance = "#009E73", Damage = "#D55E00")

f3_pt <- 7                       # base text size, pt
f3_txt <- f3_pt / ggplot2::.pt   # geom_text size, mm
f3_top <- 1.9                    # space above first row (column headers)

use_cairo <- isTRUE(capabilities("cairo"))
f3_minus <- function(x) {
  lab <- format(x, trim = TRUE, drop0trailing = TRUE)
  if (use_cairo) sub("^-", "\u2212", lab) else lab
}
f3_p <- function(p) {
  ifelse(!is.finite(p), "NA",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# Row layout: each domain gets a header row, then one row per key.
# y is plotted as -y so rows run top to bottom without a reversed scale.
f3_rows <- function(keys, labels, gap = 0.5) {
  y <- 0; rows <- list(); hdr <- list()
  for (d in f3_domains) {
    hdr[[d]] <- data.frame(domain = d, y = y)
    y <- y + 1
    rows[[d]] <- data.frame(domain = d, key = keys, label = labels,
                            y = y + seq_along(keys) - 1)
    y <- y + length(keys) + gap
  }
  list(rows = dplyr::bind_rows(rows), hdr = dplyr::bind_rows(hdr),
       ymax = y - gap)
}

f3_ylim <- function(lay) {
  ggplot2::scale_y_continuous(limits = c(-(lay$ymax + 0.4), f3_top),
                              expand = c(0, 0))
}

f3_theme_forest <- function() {
  ggplot2::theme_classic(base_size = f3_pt, base_family = "sans") +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(linewidth = 0.3),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.3),
      axis.ticks.length.x = grid::unit(1, "mm"),
      axis.text.x = ggplot2::element_text(size = 6.5, colour = "black"),
      axis.title.x = ggplot2::element_text(size = 6.5, margin = ggplot2::margin(t = 2)),
      plot.title = ggplot2::element_text(size = 7.5, face = "bold", hjust = 0.5,
                                         margin = ggplot2::margin(b = 3)),
      plot.margin = ggplot2::margin(2, 3, 2, 3),
      legend.position = "none"
    )
}
f3_theme_void <- function() {
  ggplot2::theme_void(base_size = f3_pt, base_family = "sans") +
    ggplot2::theme(plot.margin = ggplot2::margin(2, 0, 2, 0),
                   plot.title = ggplot2::element_text(size = 7.5, face = "bold",
                                                      margin = ggplot2::margin(b = 3)))
}

# Left label column (domain headers + row labels). An empty title keeps the
# panel area aligned with neighbouring forests that carry a title.
f3_label_col <- function(lay, markers = FALSE) {
  rows <- lay$rows
  p <- ggplot2::ggplot() +
    ggplot2::geom_text(data = lay$hdr,
                       ggplot2::aes(x = 0, y = -y, label = domain, colour = domain),
                       hjust = 0, fontface = "bold", size = f3_txt) +
    ggplot2::scale_colour_manual(values = f3_col, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    f3_ylim(lay) + ggplot2::labs(title = " ") + f3_theme_void()
  if (markers) {
    rows$fill <- ifelse(rows$key == "Positive", "grey30", "white")
    p <- p +
      ggplot2::geom_point(data = rows,
                          ggplot2::aes(x = 0.06, y = -y, fill = fill),
                          shape = 21, colour = "grey30", size = 1.4, stroke = 0.5) +
      ggplot2::scale_fill_identity() +
      ggplot2::geom_text(data = rows,
                         ggplot2::aes(x = 0.11, y = -y, label = label), hjust = 0, size = f3_txt)
  } else {
    rows$face <- ifelse(rows$label == "Pooled", "bold", "plain")
    p <- p + ggplot2::geom_text(data = rows,
                                ggplot2::aes(x = 0.05, y = -y, label = label, fontface = face),
                                hjust = 0, size = f3_txt)
  }
  p
}

# Forest column. dat needs: domain, y, est, lo, hi, pooled (logical),
# filled (logical). Shape 23 = pooled diamond, 21 = circle. CIs have end caps.
f3_forest <- function(dat, lay, title, xlab, log = TRUE, xlim, breaks,
                      ref = if (log) 1 else 0) {
  dat <- dat %>%
    dplyr::filter(is.finite(est), is.finite(lo), is.finite(hi)) %>%
    dplyr::mutate(col = unname(f3_col[domain]),
                  fill = ifelse(filled, col, "white"),
                  shape = ifelse(pooled, 23, 21),
                  size = ifelse(pooled, 2.0, 1.5))
  p <- ggplot2::ggplot(dat) +
    ggplot2::geom_vline(xintercept = ref, colour = "grey55",
                        linewidth = 0.3, linetype = "22") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lo, xmax = hi, y = -y, colour = col),
      orientation = "y", width = 0.35, linewidth = 0.45) +
    ggplot2::geom_point(
      ggplot2::aes(x = est, y = -y, colour = col, fill = fill,
                   shape = shape, size = size), stroke = 0.5) +
    ggplot2::scale_colour_identity() + ggplot2::scale_fill_identity() +
    ggplot2::scale_shape_identity() + ggplot2::scale_size_identity() +
    f3_ylim(lay) +
    ggplot2::labs(title = title, x = xlab) +
    f3_theme_forest()
  if (log) {
    p + ggplot2::scale_x_log10(limits = xlim, breaks = breaks,
                               labels = f3_minus, expand = c(0, 0))
  } else {
    p + ggplot2::scale_x_continuous(limits = xlim, breaks = breaks,
                                    labels = f3_minus, expand = c(0, 0))
  }
}

# Right-hand estimate / p table (panels A and B).
# Use literal UTF-8 text with fontface="bold" for the header rather than
# plotmath bold(beta), which may render Greek beta at regular weight.
# Save this R script as UTF-8 and use a Unicode-capable font/PDF device.
f3_table <- function(dat, lay, header = "aOR (95% CI)", pt = f3_pt) {
  sz <- pt / ggplot2::.pt
  dat <- dat %>% dplyr::mutate(
    face = ifelse(pooled, "bold", "plain"),
    est_lab = ifelse(is.finite(est),
                     sprintf("%.2f (%.2f, %.2f)", est, lo, hi), "Not estimable"),
    est_lab = if (use_cairo) gsub("-", "\u2212", est_lab, fixed = TRUE) else est_lab,
    p_lab = f3_p(p))
  ggplot2::ggplot(dat) +
    ggplot2::geom_text(ggplot2::aes(x = 0.04, y = -y, label = est_lab, fontface = face),
                       hjust = 0, size = sz) +
    ggplot2::geom_text(ggplot2::aes(x = 0.98, y = -y, label = p_lab, fontface = face),
                       hjust = 1, size = sz) +
    ggplot2::annotate("text", x = 0.04, y = 0.9, label = header,
                      hjust = 0, size = sz, fontface = "bold", family = "sans") +
    ggplot2::annotate("text", x = 0.98, y = 0.9, label = "italic(p)", parse = TRUE,
                      hjust = 1, size = sz) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    f3_ylim(lay) + ggplot2::labs(title = " ") + f3_theme_void()
}

# Shared, data-driven log-OR limits and breaks.
f3_or_axis <- function(lo, hi) {
  r <- range(c(lo, hi, 1), na.rm = TRUE)
  lim <- c(r[1] / 1.12, r[2] * 1.12)
  cand <- c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32)
  list(lim = lim, breaks = cand[cand >= lim[1] & cand <= lim[2]])
}

# ============================================================
# PANEL A
# ============================================================
layA <- f3_rows(f3_sources, f3_source_labels)
prepA <- function(d) {
  d %>%
    dplyr::transmute(domain = as.character(domain), key = as.character(source),
                     est = as.numeric(OR), lo = as.numeric(lower95),
                     hi = as.numeric(upper95), p = as.numeric(p)) %>%
    dplyr::inner_join(layA$rows, by = c("domain", "key")) %>%
    dplyr::mutate(pooled = key == "Pooled (cohort-adjusted)", filled = TRUE)
}
A_mews <- prepA(fig_mews)
A_mort <- prepA(fig_mort)
axA <- f3_or_axis(c(A_mews$lo, A_mort$lo), c(A_mews$hi, A_mort$hi))

panelA <- patchwork::wrap_plots(
  f3_label_col(layA),
  f3_forest(A_mews, layA, "MEWS", "aOR per 1 SD", xlim = axA$lim, breaks = axA$breaks),
  f3_table(A_mews, layA),
  f3_forest(A_mort, layA, "30-day mortality", "aOR per 1 SD", xlim = axA$lim, breaks = axA$breaks),
  f3_table(A_mort, layA),
  # Allocate additional space for full cohort labels in the left column.
  nrow = 1, widths = c(1.65, 1.35, 1.25, 1.35, 1.25))

# ============================================================
# PANEL B
# ============================================================
B_all <- fig_inf %>%
  dplyr::transmute(domain = as.character(domain), key = as.character(source),
                   infection = as.character(infection),
                   est = as.numeric(beta), lo = as.numeric(lower95),
                   hi = as.numeric(upper95), p = as.numeric(p)) %>%
  dplyr::inner_join(layA$rows, by = c("domain", "key")) %>%
  dplyr::mutate(pooled = key == "Pooled (cohort-adjusted)", filled = TRUE)

rB <- range(c(B_all$lo, B_all$hi, 0), na.rm = TRUE)
limB <- c(rB[1] - 0.1, rB[2] + 0.1)
brkB <- seq(floor(limB[1] * 2) / 2, ceiling(limB[2] * 2) / 2, by = 0.5)
brkB <- brkB[brkB >= limB[1] & brkB <= limB[2]]

forestB <- function(inf) {
  d <- B_all %>% dplyr::filter(infection == inf)
  f3_forest(d, layA, inf, "Difference (SD)", log = FALSE,
            xlim = limB, breaks = brkB) +
    # direction cue for the beta scale (full wording goes in the legend)
    ggplot2::annotate("segment", x = -0.04, xend = limB[1] + 0.03, y = 0.85, yend = 0.85,
                      linewidth = 0.25, colour = "grey45",
                      arrow = grid::arrow(length = grid::unit(0.8, "mm"), type = "closed")) +
    ggplot2::annotate("segment", x = 0.04, xend = limB[2] - 0.03, y = 0.85, yend = 0.85,
                      linewidth = 0.25, colour = "grey45",
                      arrow = grid::arrow(length = grid::unit(0.8, "mm"), type = "closed")) +
    ggplot2::annotate("text", x = limB[1] / 2, y = 1.4, label = "Lower",
                      size = 5.5 / ggplot2::.pt, colour = "grey35") +
    ggplot2::annotate("text", x = limB[2] / 2, y = 1.4, label = "Higher",
                      size = 5.5 / ggplot2::.pt, colour = "grey35")
}
tableB <- function(inf) {
  f3_table(B_all %>% dplyr::filter(infection == inf), layA,
           header = "\u03b2 (95% CI)", pt = 6)
}
panelB <- patchwork::wrap_plots(
  c(list(f3_label_col(layA)),
    unlist(lapply(f3_pathogens, function(inf) list(forestB(inf), tableB(inf))),
           recursive = FALSE)),
  # Panel B previously truncated labels; widen the cohort-label column
  # without changing the overall 180 mm figure width.
  nrow = 1, widths = c(1.85, 1.1, 1.4, 1.1, 1.4, 1.1, 1.4))

# ============================================================
# PANEL C
# ============================================================
layC <- f3_rows(c("Positive", "Negative"),
                c("Pathogen positive", "Pathogen negative"), gap = 0.4)
C_all <- fig_sub %>%
  dplyr::transmute(domain = as.character(domain), key = as.character(status),
                   infection = as.character(infection),
                   est = as.numeric(OR_per_SD), lo = as.numeric(lower95),
                   hi = as.numeric(upper95)) %>%
  dplyr::inner_join(layC$rows, by = c("domain", "key")) %>%
  dplyr::mutate(pooled = FALSE, filled = key == "Positive")
C_int <- fig_int %>%
  dplyr::transmute(domain = as.character(domain), infection = as.character(infection),
                   p_int = as.numeric(p_interaction_LRT)) %>%
  dplyr::inner_join(layC$hdr, by = "domain") %>%
  dplyr::mutate(lab = ifelse(p_int < 0.001,
                             "italic(p)[int] < 0.001",
                             sprintf("italic(p)[int] == '%s'", f3_p(p_int))))
axC <- f3_or_axis(C_all$lo, C_all$hi)

forestC <- function(inf) {
  f3_forest(C_all %>% dplyr::filter(infection == inf), layC, inf,
            "aOR for death per 1 SD", xlim = axC$lim, breaks = axC$breaks) +
    ggplot2::geom_text(data = C_int %>% dplyr::filter(infection == inf),
                       ggplot2::aes(x = axC$lim[2], y = -y, label = lab), parse = TRUE,
                       hjust = 1, size = 6.2 / ggplot2::.pt)
}
panelC <- patchwork::wrap_plots(
  c(list(f3_label_col(layC, markers = TRUE)), lapply(f3_pathogens, forestC)),
  nrow = 1, widths = c(1.25, 1.55, 1.55, 1.55))

# ============================================================
# ASSEMBLE + SAVE (180 mm full page width)
# ============================================================
tag_theme <- ggplot2::theme(plot.tag = ggplot2::element_text(size = 10, face = "bold"),
                            plot.tag.position = c(0, 1))
fig3 <- (patchwork::wrap_elements(full = panelA) + ggplot2::labs(tag = "A") + tag_theme) /
  (patchwork::wrap_elements(full = panelB) + ggplot2::labs(tag = "B") + tag_theme) /
  (patchwork::wrap_elements(full = panelC) + ggplot2::labs(tag = "C") + tag_theme) +
  patchwork::plot_layout(heights = c(layA$ymax + 3.2, layA$ymax + 3.2, layC$ymax + 3.2))

fig3_stem <- file.path(out_dir, "RESERVE_U_Figure3_combined")
ggplot2::ggsave(paste0(fig3_stem, ".pdf"), fig3, width = 180, height = 222, units = "mm",
                device = if (use_cairo) grDevices::cairo_pdf else "pdf", bg = "white")
ggplot2::ggsave(paste0(fig3_stem, ".png"), fig3, width = 180, height = 222, units = "mm",
                dpi = 600, bg = "white")
cat("\nCombined Figure 3 saved to:\n", fig3_stem, ".pdf / .png\n", sep = "")