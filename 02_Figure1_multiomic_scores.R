# =============================================================================
# Figure 1: multi-omic D/R/T domain scores in RESERVE-U-2-TOR, their application
# to RESERVE-U-1-EBB, and the transcriptomic projections
# Script: 02_Figure1_multiomic_scores.R 
# =============================================================================
#
# LAYOUT
#   Top row    : A  variance explained (TOR PC1; EBB de novo PC1; TOR-weighted
#                   score in EBB)                     | C  domain intercorrelations
#   Bottom row : B  feature-score correlations | D  transcriptomic projection coefficients |
#                E  transcriptomic projection vs multi-omic score and its transcriptome/protein components
#
# INPUT  (outputs of 01_DRT_scores.R with DERIVE = TRUE, read from SCORES_DIR)
#   derivation/multiomic_PC1_variance_explained.csv
#   derivation/multiomic_PC1_feature_score_correlations.csv
#   derivation/surrogate_coefficients.csv
#   derivation/surrogate_crossvalidated_performance.csv      (TOR, out-of-fold)
#   transfer_checks/A_variance_captured_and_tucker_phi.csv
#   transfer_checks/C_feature_score_correlations.csv
#   transfer_checks/D_transcriptomic_projection_external_performance.csv  (EBB)
#   {TOR,EBB}/{TOR,EBB}_DRT_scores_wide_perSD.csv
#
# OUTPUT
#   OUT/Figure1.{pdf,png} and OUT/panels/Figure1{A..E}_*.{pdf,png}
#
# NOTES
#   Run from the project root after 01_DRT_scores.R.
#   Axis labels are coloured by domain with vectorized element_text();
#   ggplot2 may print a harmless warning about this.
#
# REQUIRED PACKAGES
#   dplyr, tidyr, readr, ggplot2, ggprism, patchwork
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(ggprism); library(patchwork)
})

SCORES_DIR <- "DRT_scores_output"
DISCOVERY  <- "TOR"
EXTERNAL   <- "EBB"
OUT        <- file.path(SCORES_DIR, "figures")
PANELS     <- file.path(OUT, "panels")
dir.create(PANELS, recursive = TRUE, showWarnings = FALSE)

DOMAINS <- c("Damage", "Resistance", "Tolerance")
DCOL    <- c(Damage = "#E6550D", Resistance = "#2A6EBB", Tolerance = "#2CA25F")
PROT    <- c(il6 = "IL-6", il8 = "IL-8", icam1 = "ICAM-1", il10 = "IL-10",
             ho1 = "HO-1", angpt1 = "ANGPT1", angpt2 = "ANGPT2")
REVERSED_GENES    <- c("ANGPT1", "SERPINC1")                       # reverse = TRUE in 01_DRT_scores.R
REVERSED_FEATURES <- c("angpt1_prot", "ANGPT1_gene", "SERPINC1_gene")
MED_GENES <- list(Damage     = c("ANGPT2", "SERPINE1", "HMGB1", "SDC1"),
                  Resistance = c("IL6", "CXCL8", "SELE", "ICAM1"),
                  Tolerance  = c("IL10", "HMOX1", "ANGPT1", "TIMP2", "IGFBP7", "SERPINC1"))
BASE <- 9   # base font size

# Shared diverging scale for Pearson r (panels B, C, E)
scale_r <- function(show = FALSE) {
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1), na.value = "white", name = "Pearson r",
                       guide = if (show) guide_colourbar(barwidth = 8, barheight = 0.5,
                                                         title.position = "top") else "none")
}
heat_theme <- function() {
  theme_prism(base_size = BASE, border = FALSE) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "none",
          plot.tag = element_text(face = "bold", size = BASE + 4))
}
txt_col <- function(v, thr = 0.6) ifelse(!is.na(v) & abs(v) > thr, "white", "black")
rd <- function(...) read_csv(file.path(SCORES_DIR, ...), show_col_types = FALSE)

# ----------------------------------------------------------------------------- data
# Panel A: variance explained
tvar <- rd("transfer_checks", "A_variance_captured_and_tucker_phi.csv") %>% filter(cohort == EXTERNAL)
summ <- rd("derivation", "multiomic_PC1_variance_explained.csv") %>%
  transmute(domain, pct_variance_PC1_TOR = pct_variance_PC1) %>%
  inner_join(tvar %>% transmute(domain,
                                A_pct_variance_EBB_own_PC1 = pct_variance_own_PC1,
                                A_pct_variance_EBB_TOR_axis = pct_variance_discovery_axis,
                                A_TOR_axis_as_fraction_of_EBB_best = ratio_discovery_axis_vs_own),
             by = "domain")

# Panel B: feature-score correlations (TOR weights applied in both cohorts)
feat <- rd("derivation", "multiomic_PC1_feature_score_correlations.csv") %>%
  transmute(domain, feature, weight_TOR_used_for_both_cohorts = weight, corr_feature_vs_score_TOR = r_feature_score) %>%
  inner_join(rd("transfer_checks", "C_feature_score_correlations.csv") %>% filter(cohort == EXTERNAL) %>%
               transmute(domain, feature, corr_feature_vs_score_EBB = r_feature_score), by = c("domain", "feature")) %>%
  mutate(
    direction_reversed = feature %in% REVERSED_FEATURES,
    is_prot = grepl("_prot$", feature),
    label = ifelse(is_prot,
                   paste0(PROT[sub("_prot$", "", feature)], " (protein)"),
                   paste0(sub("_gene$", "", feature), " (gene)")),
    label = paste0(label, ifelse(direction_reversed, "†", "")),
    # axis label as plotmath: protein name upright, gene symbol italic
    label_expr = paste0(ifelse(is_prot,
                               paste0('"', PROT[sub("_prot$", "", feature)], ' (protein)"'),
                               paste0('italic("', sub("_gene$", "", feature), '")~"(gene)"')),
                        ifelse(direction_reversed, '*"†"', "")),
    domain = factor(domain, levels = DOMAINS)
  ) %>%
  arrange(domain, desc(is_prot), desc(abs(corr_feature_vs_score_TOR)))

# Panel C: domain intercorrelations (multi-omic scores)
inter <- bind_rows(lapply(c(DISCOVERY, EXTERNAL), function(c) {
  w <- rd(c, paste0(c, "_DRT_scores_wide_perSD.csv"))
  m <- na.omit(as.data.frame(w[, paste0("multiomic_PC1__", DOMAINS)]))
  r <- cor(m); dimnames(r) <- list(DOMAINS, DOMAINS)
  expand_grid(row = DOMAINS, col = DOMAINS) %>%
    mutate(r = r[cbind(row, col)], cohort = sprintf("%s (n = %d)", c, nrow(m)))
})) %>% mutate(cohort = factor(cohort, levels = unique(cohort)))   # TOR first, then EBB

# Panel D: locked LASSO coefficients (unreversed scale in the file)
sel <- rd("derivation", "surrogate_coefficients.csv") %>% transmute(domain, symbol = gene, coefficient)

# Panel E: TOR out-of-fold (first CV repeat) and EBB external performance
perf <- bind_rows(
  rd("derivation", "surrogate_crossvalidated_performance.csv") %>%
    transmute(domain, cohort_short = DISCOVERY, pearson_r = oof_pearson_r,
              r_with_rna_part = oof_r_with_rna_component, r_with_protein_part = oof_r_with_protein_component,
              r_target_vs_protein_part = score_r_with_protein_component),
  rd("transfer_checks", "D_transcriptomic_projection_external_performance.csv") %>% filter(cohort == EXTERNAL) %>%
    transmute(domain, cohort_short = EXTERNAL, pearson_r,
              r_with_rna_part = r_with_rna_component, r_with_protein_part = r_with_protein_component,
              r_target_vs_protein_part = score_r_with_protein_component))

# ----------------------------------------------------------------------------- A
bars <- summ %>%
  transmute(domain = factor(domain, levels = DOMAINS),
            `TOR: TOR PC1` = pct_variance_PC1_TOR,
            `EBB: EBB de novo PC1` = A_pct_variance_EBB_own_PC1,
            `EBB: TOR-weighted score` = A_pct_variance_EBB_TOR_axis,
            ratio = A_TOR_axis_as_fraction_of_EBB_best) %>%
  pivot_longer(c(`TOR: TOR PC1`, `EBB: EBB de novo PC1`, `EBB: TOR-weighted score`),
               names_to = "score", values_to = "pct") %>%
  mutate(score = factor(score, levels = c("TOR: TOR PC1", "EBB: EBB de novo PC1",
                                          "EBB: TOR-weighted score")))
# Group-level annotation above each domain (avoids overlap with bars/value labels)
ratio_lab <- bars %>% filter(score == "EBB: TOR-weighted score") %>%
  mutate(lab = sprintf("TOR-weighted score:\n%.0f%% of EBB de novo PC1", 100 * ratio))
dodge <- position_dodge(width = 0.9)
pA <- ggplot(bars, aes(domain, pct, fill = score)) +
  geom_col(position = dodge, width = 0.88, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = dodge, vjust = -0.5, size = 2.7) +
  geom_text(data = ratio_lab, aes(x = domain, y = 51, label = lab), inherit.aes = FALSE,
            size = 2.6, colour = "#3F007D", lineheight = 0.9, vjust = 1) +
  scale_fill_manual(values = c("#4D4D4D", "#BDBDBD", "#8C96C6"), name = NULL) +
  scale_y_continuous(limits = c(0, 52), breaks = seq(0, 40, 10), expand = c(0, 0)) +
  labs(x = NULL, y = "Variance explained (%)", tag = "A") +
  theme_prism(base_size = BASE) +
  theme(legend.position = "top", legend.text = element_text(size = BASE - 1),
        axis.text.x = element_text(colour = DCOL[DOMAINS], face = "bold"),
        plot.tag = element_text(face = "bold", size = BASE + 4))

# ----------------------------------------------------------------------------- B
featB <- feat %>% mutate(row = factor(label, levels = rev(label)))
nB <- nrow(featB)
longB <- bind_rows(
  featB %>% transmute(row, domain, col = "TOR", r = corr_feature_vs_score_TOR),
  featB %>% transmute(row, domain, col = "EBB", r = corr_feature_vs_score_EBB)
) %>% mutate(col = factor(col, levels = c("TOR", "EBB")))
n_by_dom <- as.integer(table(featB$domain))
hlinesB <- nB - cumsum(n_by_dom)[1:2] + 0.5
# The TOR-weight column is drawn as text OUTSIDE the panel (clip = "off"), so the
# panel -- and therefore the colour bar, centred on the panel -- covers only the
# two coloured columns.
pB <- ggplot(longB, aes(col, row)) +
  geom_tile(aes(fill = r), colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r), colour = I(txt_col(r))), size = 2.6) +
  geom_text(data = featB, aes(x = 3, y = row, label = sprintf("%.2f", weight_TOR_used_for_both_cohorts)),
            inherit.aes = FALSE, size = 2.6) +
  annotate("text", x = 3, y = nB + 0.8, label = "TOR\nweight", fontface = "bold",
           size = BASE / .pt, vjust = 0, lineheight = 0.9) +
  annotate("segment", x = 0.5, xend = 3.45, y = hlinesB, yend = hlinesB, linewidth = 0.4) +
  scale_r(show = TRUE) + scale_x_discrete(position = "top") +
  scale_y_discrete(labels = function(x) parse(text = setNames(featB$label_expr, featB$label)[x])) +
  coord_cartesian(xlim = c(0.5, 2.5), ylim = c(0.5, nB + 0.5), expand = FALSE, clip = "off") +
  labs(x = NULL, y = NULL, tag = "B") +
  heat_theme() +
  theme(legend.position = "bottom", legend.location = "panel",
        legend.title = element_text(size = BASE - 1),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(colour = rev(DCOL[as.character(featB$domain)])),
        plot.margin = margin(5, 45, 5, 5))

# ----------------------------------------------------------------------------- C
lowerC <- inter %>% filter(match(row, DOMAINS) > match(col, DOMAINS)) %>%
  mutate(row = factor(row, levels = rev(DOMAINS[-1])), col = factor(col, levels = DOMAINS[-3]))
pC <- ggplot(lowerC, aes(col, row, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_r() + coord_equal() +
  labs(x = NULL, y = NULL, tag = "C") +
  heat_theme() +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1, colour = DCOL[DOMAINS[1:2]]),
        axis.text.y = element_text(colour = DCOL[rev(DOMAINS[-1])]))

# ----------------------------------------------------------------------------- D
genes <- unlist(MED_GENES, use.names = FALSE)
gdom  <- rep(names(MED_GENES), lengths(MED_GENES))
coefD <- expand_grid(gene = genes, model = DOMAINS) %>%
  left_join(sel %>% transmute(gene = symbol, model = domain,
                              # direction-coded scale as in B; identical model
                              coef = ifelse(symbol %in% REVERSED_GENES, -1, 1) * coefficient),
            by = c("gene", "model")) %>%
  mutate(glabel = paste0(gene, ifelse(gene %in% REVERSED_GENES, "†", "")),
         glabel = factor(glabel, levels = rev(paste0(genes, ifelse(genes %in% REVERSED_GENES, "†", "")))),
         model = factor(model, levels = DOMAINS))
n_ret <- coefD %>% group_by(model) %>% summarise(n = sum(!is.na(coef)), .groups = "drop")
pD <- ggplot(coefD, aes(model, glabel)) +
  geom_tile(aes(fill = coef), colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(coef), "", sprintf("%.2f", coef)),
                colour = I(txt_col(coef, 0.45))), size = 2.6) +
  geom_hline(yintercept = c(length(genes) - 4 + 0.5, length(genes) - 8 + 0.5), linewidth = 0.4) +
  scale_fill_gradient2(low = "#542788", mid = "white", high = "#B35806", midpoint = 0,
                       limits = c(-0.7, 0.7), na.value = "#F2F2F2", guide = "none") +
  scale_x_discrete(position = "top") +
  # gene symbols in italics (plotmath)
  scale_y_discrete(labels = function(x) parse(text = paste0('italic("', sub("†$", "", x), '")',
                                                            ifelse(grepl("†$", x), '*"†"', "")))) +
  annotate("text", x = 2, y = 0.2, vjust = 1, size = (BASE - 1) / .pt, lineheight = 0.95,
           label = paste0("Transcriptomic projection model coefficients\nGenes retained: ",
                          paste(sprintf("%d/14", n_ret$n), collapse = ", "))) +
  coord_cartesian(ylim = c(0.5, length(genes) + 0.5), expand = FALSE, clip = "off") +
  labs(x = NULL, y = NULL, tag = "D") +
  heat_theme() +
  theme(axis.text.x.top = element_text(face = "bold", colour = DCOL[DOMAINS], size = BASE - 1),
        plot.margin = margin(5, 5, 30, 5),   # room below the tiles for the note
        axis.text.y = element_text(colour = rev(DCOL[gdom])))

# ----------------------------------------------------------------------------- E
cols_E <- c("Transcriptomic\nprojection vs\nmulti-omic\nscore", "Transcriptomic\nprojection vs\ntranscriptome\ncomponent",
            "Transcriptomic\nprojection vs\nprotein\ncomponent", "Multi-omic\nscore vs\nprotein\ncomponent")
rowsE <- expand_grid(domain = DOMAINS, cs = c("TOR", "EBB")) %>%
  mutate(rlab = paste0(domain, ", ", cs, ifelse(cs == "TOR", " (OOF)", "")))
longE <- rowsE %>%
  left_join(perf %>% dplyr::select(domain, cs = cohort_short, pearson_r, r_with_rna_part,
                            r_with_protein_part, r_target_vs_protein_part),
            by = c("domain", "cs")) %>%
  pivot_longer(c(pearson_r, r_with_rna_part, r_with_protein_part, r_target_vs_protein_part),
               names_to = "metric", values_to = "r") %>%
  mutate(metric = factor(metric, levels = c("pearson_r", "r_with_rna_part",
                                            "r_with_protein_part", "r_target_vs_protein_part"),
                         labels = cols_E),
         rlab = factor(rlab, levels = rev(rowsE$rlab)))
pE <- ggplot(longE, aes(metric, rlab, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r), colour = I(txt_col(r))), size = 2.8) +
  geom_vline(xintercept = 3.5, linewidth = 0.4) +
  geom_hline(yintercept = c(2.5, 4.5), linewidth = 0.4) +
  scale_r() + scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, tag = "E") +
  heat_theme() +
  theme(axis.text.x.top = element_text(size = BASE - 3.5, lineheight = 0.85),   # smaller: longer header text
        axis.text.y = element_text(colour = rev(DCOL[rowsE$domain])))

# ----------------------------------------------------------------------------- layout & save
# cairo_pdf embeds the special characters (dagger) reliably; fall back to pdf.
PDF_DEVICE <- if (isTRUE(capabilities("cairo"))) cairo_pdf else "pdf"
# Each row is wrapped so that panels are aligned WITHIN a row only; otherwise
# patchwork aligns panel A's plotting area with B's (long y labels), pushing
# A's y-axis title away from its axis.
top    <- (pA | pC) + plot_layout(widths = c(2.4, 1.3))
bottom <- (pB | pD | pE) + plot_layout(widths = c(0.75, 0.85, 1.3))
fig1   <- wrap_elements(full = top) / wrap_elements(full = bottom) + plot_layout(heights = c(1, 1.95))

for (ext in c("pdf", "png"))
  ggsave(file.path(OUT, paste0("Figure1.", ext)), fig1, width = 11, height = 9.4,
         dpi = 300, bg = "white", device = if (ext == "pdf") PDF_DEVICE else "png")

panel_list <- list(A_variance_explained = list(pA, 7, 3.6),
                   B_feature_score_correlations = list(pB, 4.2, 6.2),
                   C_domain_intercorrelations = list(pC, 4.2, 2.8),
                   D_transcriptomic_projection_coefficients = list(pD, 3.4, 6.0),
                   E_transcriptomic_projection_vs_score_components = list(pE, 4.6, 4.8))
for (nm in names(panel_list)) {
  p <- panel_list[[nm]]
  for (ext in c("pdf", "png"))
    ggsave(file.path(PANELS, paste0("Figure1", nm, ".", ext)), p[[1]], width = p[[2]], height = p[[3]],
           dpi = 300, bg = "white", device = if (ext == "pdf") PDF_DEVICE else "png")
}
message("Saved Figure1 and panels to ", OUT)
