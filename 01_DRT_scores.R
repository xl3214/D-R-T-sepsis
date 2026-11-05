# =============================================================================
# Damage / Resistance / Tolerance (D/R/T) domain scores
# Script: 01_DRT_scores.R 
# =============================================================================
#
# WHAT THIS SCRIPT DOES
#   Computes patient- (or sample-) level D/R/T domain scores with the two
#   approaches described in the manuscript, from raw input files:
#
#   1. Multi-omic score ("multiomic")            [requires Olink + RNA-seq]
#      PC1 of the 14 core Chowdhury mediators (Olink protein + transcriptome), features
#      z-scored within each cohort. PC1 weights are estimated in the discovery
#      cohort (TOR) and applied unchanged to other cohorts.
#   2. Transcriptomic projection of the multi-omic score ("surrogate")  [RNA-seq]
#      LASSO model trained in TOR to predict the multi-omic score from the 14
#      mediator genes (lambda.1se); coefficients locked.
#
# TWO WAYS TO RUN
#   DERIVE = TRUE  : re-estimates every locked quantity (multi-omic weights and
#                    LASSO coefficients) from the discovery
#                    cohort, compares them with the published values embedded
#                    below (Section 2), and scores all cohorts with the
#                    re-derived values. Requires the discovery cohort's data.
#   DERIVE = FALSE : scores any cohort with the published locked values. Use
#                    this for new/external cohorts.
#   Either way, each cohort can be scored with any subset of the approaches
#   (setting SCORES); "multiomic" is only computed for cohorts with Olink data.
#
# STANDARDIZATION (identical to the manuscript)
#   Multi-omic: each feature z-scored within the cohort, among samples with
#     complete data for that domain.
#   Transcriptomic projection: each gene z-scored against the cohort's REFERENCE
#     samples (default: all scored samples; for longitudinal cohorts, the
#     baseline septic samples), and the same transformation applied to all
#     samples, so later timepoints and controls share one scale.
#   Every score is reported as score_raw and as score_perSD (mean 0, SD 1 in
#     the reference samples).
#
# RNA-SEQ PROCESSING (per cohort, setting rna_processing)
#   "vst"         raw counts with Ensembl IDs -> Entrez IDs (duplicates: keep
#                 the Ensembl ID with the highest median count) -> drop genes
#                 with <= 2 counts in every sample -> DESeq2 blind VST -> drop
#                 zero-variance genes. Used for TOR and EBB.
#   "tmm_log2cpm" raw counts (gene symbols or Entrez IDs; duplicate rows
#                 summed) -> TMM-normalized log2-CPM (edgeR, prior.count = 1).
#                 Used for ACESO-Ghana and GSE185263.
#   "none"        already-normalized log-scale matrix, used as provided.
#
# INPUT FILES
#   (CSV or tab-delimited .tsv/.txt files are both accepted.)
#   RNA-seq    genes x samples; first column = gene ID; column names = sample IDs.
#   Olink      samples x proteins; one ID column; NPX values in columns named
#              il6, il8, icam1, il10, ho1, angpt1, angpt2 (only these are used).
#   Metadata   one row per sample; ID column matching the RNA-seq column names.
#              Samples absent from the metadata are not scored. Rows are
#              processed in metadata order (this order determines the
#              cross-validation folds when DERIVE = TRUE).
#
# REPRODUCIBILITY
#   All steps are deterministic except the LASSO cross-validation (DERIVE =
#   TRUE), which uses fixed seeds (20261002 + repeat). Results were produced
#   with R 4.4.1, DESeq2 1.44.0, edgeR 4.2.2, glmnet 4.1-10 and
#   org.Hs.eg.db 3.19.1; other versions may give small numerical differences.
#   With DERIVE = TRUE the script reports the maximum difference between the
#   re-derived and the published values (derivation_vs_published.csv).
#
# OUTPUT  (<OUTPUT_DIR>/)
#   <cohort>/<cohort>_DRT_scores_long.csv        sample x method x domain
#   <cohort>/<cohort>_DRT_scores_wide_perSD.csv  sample x (method__domain)
#   <cohort>/gene_coverage.csv
#   derivation/  (DERIVE = TRUE) PC1 weights, variance explained and feature-score
#                correlations in the discovery cohort, coefficients,
#                transcriptomic projection cross-validated performance, comparison with the
#                published values
#   transfer_checks/  (multi-omic, non-discovery cohorts) variance captured,
#                PC1/PC2 separation (North's rule), feature-score correlations,
#                external performance of the transcriptomic projection
#   method_dictionary.csv, run_settings.csv, console_log.txt, sessionInfo.txt
#
# REQUIRED PACKAGES
#   tidyverse; DESeq2 ("vst"); edgeR ("tmm_log2cpm"); AnnotationDbi and
#   org.Hs.eg.db (Ensembl/symbol mapping); glmnet (DERIVE = TRUE).
# =============================================================================


# =============================================================================
# SECTION 1. CONFIGURATION -- edit paths and options here
# =============================================================================
OUTPUT_DIR <- "DRT_scores_output"
SCORES     <- c("multiomic", "surrogate")            # either or both
DERIVE     <- TRUE                                    # TRUE: re-derive from DISCOVERY; FALSE: use published values
DISCOVERY  <- "TOR"                                   # cohort used when DERIVE = TRUE
OOF_METHOD <- "first_repeat"   # TOR out-of-fold predictions used as the TOR projection scores: "first_repeat"
                               # (CV repeat 1; as published) or "average_10" (mean over the 10 CV repeats)

# One entry per cohort. Fields:
#   rna_file, rna_gene_id ("ensembl" | "entrez" | "symbol"),
#   rna_processing ("vst" | "tmm_log2cpm" | "none"),
#   olink_file (NULL if none), olink_id_col,
#   meta_file, meta_id_col, keep_cols (metadata columns copied to the output),
#   reference (NULL = all scored samples; or a named list of metadata
#              conditions, e.g. list(COHORT = "S", COLL_TIME = "00h"))
COHORTS <- list(
  TOR = list(
    rna_file = "data/tororo_all_274_raw_ensembl_only.csv", rna_gene_id = "ensembl", rna_processing = "vst",
    olink_file = "data/tororo_all_287_olink_proteomics.csv", olink_id_col = "pid",
    meta_file = "data/tororo_all_298_clinical_data.csv", meta_id_col = "pid", keep_cols = NULL,
    reference = NULL),
  EBB = list(
    rna_file = "data/cohort2_entebbe_all_128_raw_ensembl.csv", rna_gene_id = "ensembl", rna_processing = "vst",
    olink_file = "data/entebbe_all_278_olink_proteomics.csv", olink_id_col = "pid",
    meta_file = "data/entebbe_all_301_clinical_data.csv", meta_id_col = "pid", keep_cols = NULL,
    reference = NULL),
  Ghana = list(
    rna_file = "Analyses_Sept_2026/ACESO_Ghana/ACESO_GHANA_COGAPSPUB_PMID38816375_N526FUA_RAWCOUNTS_18DEC25.csv",
    rna_gene_id = "symbol", rna_processing = "tmm_log2cpm",
    olink_file = NULL, olink_id_col = NULL,
    meta_file = "Analyses_Sept_2026/ACESO_Ghana/ACESO_GHANA_COGAPSPUB_PMID38816375_N526FUA_METACLINICAL_18DEC25.csv",
    meta_id_col = "SAMPLE_ID", keep_cols = c("SUBJECT_ID", "COLL_TIME", "COHORT"),
    reference = list(COHORT = "S", COLL_TIME = "00h")), 
  GSE185263 = list(
    rna_file = "Multinational_sepsis_GSE185263/GSE185263_Raw_Counts_Entrez.csv",
    rna_gene_id = "entrez", rna_processing = "tmm_log2cpm",
    olink_file = NULL, olink_id_col = NULL,
    meta_file = "Multinational_sepsis_GSE185263/GSE185263_metadata.csv",
    meta_id_col = "geo_accession",
    keep_cols = c("title", "disease state:ch1", "collection site:ch1", "in hospital mortality:ch1"),
    reference = list("disease state:ch1" = "sepsis"))
)


# =============================================================================
# SECTION 2. PUBLISHED (LOCKED) DEFINITIONS -- do not edit
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble); library(purrr); library(readr) })

DOMAINS <- c("Damage", "Resistance", "Tolerance")
SEED_BASE <- 20261002
N_CV_REPEATS <- 10
N_FOLDS <- 10

# 2a. Core mediators (Chowdhury et al.). reverse = multiplied by -1 so that
#     higher values indicate greater domain activity.
CORE_MEDIATORS <- tribble(
  ~domain,      ~mediator,          ~gene,      ~entrez, ~protein,  ~reverse,
  "Resistance", "IL-6",             "IL6",      "3569",  "il6",     FALSE,
  "Resistance", "IL-8",             "CXCL8",    "3576",  "il8",     FALSE,
  "Resistance", "E-selectin",       "SELE",     "6401",  NA,        FALSE,
  "Resistance", "ICAM-1",           "ICAM1",    "3383",  "icam1",   FALSE,
  "Tolerance",  "IL-10",            "IL10",     "3586",  "il10",    FALSE,
  "Tolerance",  "HO-1",             "HMOX1",    "3162",  "ho1",     FALSE,
  "Tolerance",  "Angiopoietin-1",   "ANGPT1",   "284",   "angpt1",  TRUE,
  "Tolerance",  "TIMP-2",           "TIMP2",    "7077",  NA,        FALSE,
  "Tolerance",  "IGFBP7",           "IGFBP7",   "3490",  NA,        FALSE,
  "Tolerance",  "Antithrombin III", "SERPINC1", "462",   NA,        TRUE,
  "Damage",     "Angiopoietin-2",   "ANGPT2",   "285",   "angpt2",  FALSE,
  "Damage",     "Syndecan-1",       "SDC1",     "6382",  NA,        FALSE,
  "Damage",     "HMGB1",            "HMGB1",    "3146",  NA,        FALSE,
  "Damage",     "PAI-1",            "SERPINE1", "5054",  NA,        FALSE
)

# 2b. Published multi-omic PC1 weights (estimated in TOR). Feature names:
#     <protein>_prot (Olink NPX) or <GENE>_gene (RNA-seq), after reversal.
PUBLISHED_MULTIOMIC_WEIGHTS <- tribble(
  ~domain, ~feature, ~weight,
  "Resistance", "il6_prot", 0.5845562181070707,
  "Resistance", "IL6_gene", -0.067710891454596,
  "Resistance", "il8_prot", 0.4737835583917897,
  "Resistance", "CXCL8_gene", 0.1462595487386211,
  "Resistance", "SELE_gene", -0.1041296380915812,
  "Resistance", "icam1_prot", 0.4479955999382397,
  "Resistance", "ICAM1_gene", 0.4430615173517834,
  "Tolerance", "il10_prot", 0.3760706739744338,
  "Tolerance", "IL10_gene", 0.4742453207336217,
  "Tolerance", "ho1_prot", 0.1782179832955976,
  "Tolerance", "HMOX1_gene", 0.4781611286512052,
  "Tolerance", "angpt1_prot", 0.1473351478702846,
  "Tolerance", "ANGPT1_gene", 0.1206562581564478,
  "Tolerance", "TIMP2_gene", 0.3436443733718673,
  "Tolerance", "IGFBP7_gene", 0.4202063973349674,
  "Tolerance", "SERPINC1_gene", -0.2057475620124946,
  "Damage", "angpt2_prot", 0.3744409435883756,
  "Damage", "ANGPT2_gene", -0.5929587317568054,
  "Damage", "SDC1_gene", 0.2153696028692719,
  "Damage", "HMGB1_gene", 0.4067820528109011,
  "Damage", "SERPINE1_gene", 0.5443695600131311
)

# 2c. Published LASSO transcriptomic projection coefficients (14 mediator genes, lambda.1se).
#     Intercept = 0 (training outcome had mean 0; predictors z-scored).
#     NOTE: coefficients apply to UNREVERSED gene expression, because the
#     script applies them to z-scored expression as measured. In the manuscript
#     (Figure 1D), ANGPT1 and SERPINC1 are shown on the reversed (direction-
#     coded) scale, so their signs appear opposite to the values below
#     (ANGPT1 +0.06, SERPINC1 -0.18). Both describe the identical model.
PUBLISHED_SURROGATE <- tribble(
  ~domain, ~gene, ~entrez, ~coefficient,
  "Damage", "ANGPT2", "285", -0.6362445878148799,
  "Damage", "SERPINE1", "5054", 0.5522134525094885,
  "Damage", "HMGB1", "3146", 0.35026986193876,
  "Damage", "SDC1", "6382", 0.192660346984681,
  "Damage", "ICAM1", "3383", 0.0276372581643196,
  "Damage", "IL10", "3586", 0.0111121160953912,
  "Damage", "TIMP2", "7077", 0.0051025471416385,
  "Resistance", "ICAM1", "3383", 0.6857080150918311,
  "Resistance", "IL10", "3586", 0.1184398312090858,
  "Resistance", "ANGPT2", "285", -0.1053326471097108,
  "Resistance", "IGFBP7", "3490", 0.0676087319217101,
  "Tolerance", "IL10", "3586", 0.6480526188670981,
  "Tolerance", "HMOX1", "3162", 0.5778118648687163,
  "Tolerance", "IGFBP7", "3490", 0.4042543098568097,
  "Tolerance", "SERPINC1", "462", 0.1804281462753437,
  "Tolerance", "ICAM1", "3383", 0.1601615857412286,
  "Tolerance", "TIMP2", "7077", 0.0953976118207809,
  "Tolerance", "ANGPT1", "284", -0.0640955388245968,
  "Tolerance", "SDC1", "6382", 0.045635528489735,
  "Tolerance", "IL6", "3569", 0.0205085825007732
)
PUBLISHED_SURROGATE_LAMBDA <- c(Damage = 0.05672156121290934, Resistance = 0.2074664744862773, Tolerance = 0.057071996612860394)

METHOD_INFO <- tribble(
  ~method,                           ~approach,  ~role,
  "multiomic_PC1",                   "multiomic", "Primary multi-omic score",
  "rna_surrogate_lasso14_1se",       "surrogate", "Transcriptomic projection of the multi-omic score",
  "rna_surrogate_lasso14_1se_OOF",   "surrogate", paste0("Discovery cohort only: out-of-fold transcriptomic projection predictions (", OOF_METHOD, "; DERIVE = TRUE)")
)


# =============================================================================
# SECTION 3. SETUP
# =============================================================================
need <- c(if (any(map_chr(COHORTS, "rna_processing") == "vst")) c("DESeq2", "AnnotationDbi", "org.Hs.eg.db"),
          if (any(map_chr(COHORTS, "rna_processing") == "tmm_log2cpm")) "edgeR",
          if (DERIVE && "surrogate" %in% SCORES) "glmnet")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Please install: ", paste(miss, collapse = ", "))
stopifnot(all(SCORES %in% c("multiomic", "surrogate")))
if (DERIVE && !DISCOVERY %in% names(COHORTS)) stop("DISCOVERY cohort '", DISCOVERY, "' is not defined in COHORTS.")
if (DERIVE && is.null(COHORTS[[DISCOVERY]]$olink_file)) stop("DERIVE = TRUE requires Olink data for the discovery cohort.")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
log_con <- file(file.path(OUTPUT_DIR, "console_log.txt"), open = "wt")
sink(log_con, type = "output", split = TRUE); sink(log_con, type = "message")
close_log <- function() { sink(type = "message"); sink(type = "output"); close(log_con) }
# If the script stops with an error, run close_log() to restore console output.
message("D/R/T scoring -- ", format(Sys.time()), " -- DERIVE = ", DERIVE, " -- SCORES = ", paste(SCORES, collapse = ", "))

normalize_symbol <- function(x) toupper(gsub("\\s+", "", trimws(as.character(x))))
# Reads CSV, or tab-delimited files (.tsv / .txt / .tsv.gz), keeping column names as is
read_table_auto <- function(path, ...) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (grepl("\\.(tsv|txt)(\\.gz)?$", path, ignore.case = TRUE))
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
  else utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
}
read_matrix <- function(path) {
  dd <- read_table_auto(path, row.names = NULL)
  ids <- trimws(as.character(dd[[1]]))
  m <- as.matrix(data.frame(lapply(dd[, -1, drop = FALSE], function(x) suppressWarnings(as.numeric(x))), check.names = FALSE))
  colnames(m) <- trimws(colnames(dd)[-1]); rownames(m) <- NULL
  list(ids = ids, m = m)
}


# =============================================================================
# SECTION 4. RNA-SEQ PROCESSING
# =============================================================================
# "vst": project-standard pipeline used for TOR and EBB. The four functions
# below reproduce, line for line, the shared helper functions used in the
# original analysis (ensemblToNCBI, prepRNAseq, normalizeRNASEQwithVST,
# removeNoVarRows), so that VST values are identical.
ENSEMBL_KEYS <- NULL
ensemblToNCBI <- function(ensIds) {
  if (is.null(ENSEMBL_KEYS)) ENSEMBL_KEYS <<- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBL")
  if (sum(grepl("ENSG", ensIds)) != length(ensIds)) { warning("Non-Ensembl gene IDs present"); return(NULL) }
  ensIds <- sapply(strsplit(ensIds, ".", fixed = TRUE), function(t) t[1])      # remove version suffix
  ensIds <- ensIds[ensIds %in% ENSEMBL_KEYS]                                     # valid Ensembl IDs only
  if (length(ensIds) == 0) return(NULL)
  map <- as.matrix(suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = ensIds,
                                                          keytype = "ENSEMBL", columns = c("ENTREZID"))))
  t <- which(duplicated(map[, 1])); if (length(t) > 0) map <- map[-t, , drop = FALSE]
  res <- map[, 2]; names(res) <- map[, 1]
  res
}
prepRNAseq <- function(f_name, thresh = 2) {
  dd <- utils::read.csv(f_name, row.names = 1, check.names = FALSE)
  eids <- ensemblToNCBI(rownames(dd))
  # several Ensembl IDs mapping to one gene ID: keep the one with the highest median count
  dups <- duplicated(eids)
  if (sum(dups) > 0) {
    dups <- unique(eids[duplicated(eids)])
    drop <- unlist(lapply(dups, function(eid) {
      xx <- apply(dd[names(eids)[eids == eid], ], 1, stats::median)
      names(sort(xx, decreasing = TRUE))[2:length(xx)]
    }))
    eids <- eids[setdiff(names(eids), drop)]
  }
  ddn <- dd[names(eids), ]; rownames(ddn) <- eids
  keep <- apply(ddn, 1, function(r) sum(r > thresh) > 0)                        # drop genes with <= 2 counts in every sample
  ddn[keep, ]
}
normalizeRNASEQwithVST <- function(readCount) {
  condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = readCount, colData = condition, design = ~0))
  data.frame(SummarizedExperiment::assay(DESeq2::vst(dds)), check.names = FALSE)   # blind VST
}
removeNoVarRows <- function(dmat, min_var = 0) dmat[!(apply(dmat, 1, stats::var) <= min_var), ]
process_vst <- function(path) as.matrix(removeNoVarRows(normalizeRNASEQwithVST(prepRNAseq(path))))

process_tmm <- function(path, gene_id) {
  x <- read_matrix(path)
  ids <- x$ids
  if (gene_id == "symbol") { if (any(grepl("\\|", ids))) ids <- sub("^.*\\|", "", ids); ids <- normalize_symbol(ids) }
  if (gene_id == "entrez") ids <- sub("\\.0+$", "", ids)            # e.g. "3569.0" -> "3569"
  m <- x$m; m[is.na(m)] <- 0
  keep <- !is.na(ids) & ids != ""
  m <- rowsum(m[keep, , drop = FALSE], group = ids[keep], reorder = FALSE)
  edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts = m), method = "TMM"), log = TRUE, prior.count = 1)
}

process_none <- function(path, gene_id) {
  x <- read_matrix(path)
  ids <- if (gene_id == "symbol") normalize_symbol(x$ids) else sub("\\..*$", "", x$ids)
  m <- x$m[order(rowMeans(x$m, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
  ids <- ids[order(rowMeans(x$m, na.rm = TRUE), decreasing = TRUE)]
  m <- m[!duplicated(ids), , drop = FALSE]; rownames(m) <- ids[!duplicated(ids)]
  m
}

# Returns list(expr = genes x samples, key = "entrez" | "symbol")
load_rna <- function(cfg) {
  if (cfg$rna_processing == "vst") {
    if (cfg$rna_gene_id != "ensembl") stop("rna_processing = 'vst' expects Ensembl gene IDs.")
    list(expr = process_vst(cfg$rna_file), key = "entrez")
  } else if (cfg$rna_processing == "tmm_log2cpm") {
    if (cfg$rna_gene_id == "ensembl") stop("For 'tmm_log2cpm' use gene symbols or Entrez IDs.")
    list(expr = process_tmm(cfg$rna_file, cfg$rna_gene_id), key = cfg$rna_gene_id)
  } else {
    if (cfg$rna_gene_id == "ensembl") stop("For 'none' use gene symbols or Entrez IDs.")
    list(expr = process_none(cfg$rna_file, cfg$rna_gene_id), key = cfg$rna_gene_id)
  }
}


# =============================================================================
# SECTION 5. LOAD COHORTS
# =============================================================================
data <- imap(COHORTS, function(cfg, nm) {
  message("\n=== Loading ", nm, " ===")
  rna <- load_rna(cfg)
  if (!file.exists(cfg$meta_file)) stop(nm, ": metadata file not found: ", cfg$meta_file)
  meta <- read_table_auto(cfg$meta_file)
  need_cols <- c(cfg$meta_id_col, cfg$keep_cols, names(cfg$reference))
  absent <- setdiff(need_cols, names(meta))
  if (length(absent)) stop(nm, ": column(s) not found in the metadata file: ", paste(absent, collapse = ", "),
                           ". Available columns: ", paste(names(meta), collapse = ", "))
  meta[[cfg$meta_id_col]] <- trimws(as.character(meta[[cfg$meta_id_col]]))
  meta <- meta[!duplicated(meta[[cfg$meta_id_col]]), , drop = FALSE]
  ids <- meta[[cfg$meta_id_col]][meta[[cfg$meta_id_col]] %in% colnames(rna$expr)]   # metadata order
  if (!length(ids)) stop(nm, ": no metadata IDs ('", cfg$meta_id_col, "') match the RNA-seq column names (e.g. ",
                         paste(head(colnames(rna$expr), 3), collapse = ", "), ").")
  expr <- rna$expr[, ids, drop = FALSE]
  ref <- ids
  if (!is.null(cfg$reference)) {
    keep <- rep(TRUE, length(ids)); mm <- meta[match(ids, meta[[cfg$meta_id_col]]), , drop = FALSE]
    for (v in names(cfg$reference)) keep <- keep & trimws(as.character(mm[[v]])) %in% cfg$reference[[v]]
    ref <- ids[keep]
    if (!length(ref)) stop(nm, ": no samples match the reference definition. Values must match exactly (case-sensitive); ",
                           "check e.g. unique(meta[[\"", names(cfg$reference)[1], "\"]]).")
  }
  olink <- NULL
  if (!is.null(cfg$olink_file)) {
    olink <- read_table_auto(cfg$olink_file)
    olink[[cfg$olink_id_col]] <- trimws(as.character(olink[[cfg$olink_id_col]]))
    olink <- olink[!duplicated(olink[[cfg$olink_id_col]]), , drop = FALSE]
  }
  message(sprintf("  RNA-seq: %d genes; %d samples scored (with metadata); %d reference samples",
                  nrow(expr), length(ids), length(ref)))
  if (length(ref) < 20) warning(nm, ": fewer than 20 reference samples.")
  list(cfg = cfg, expr = expr, key = rna$key, ids = ids, ref = ref, meta = meta, olink = olink)
})


# =============================================================================
# SECTION 6. SCORING FUNCTIONS
# =============================================================================
# z-score rows of a genes x samples matrix against reference columns
z_against <- function(X, ref) {
  mu <- rowMeans(X[, ref, drop = FALSE]); s <- apply(X[, ref, drop = FALSE], 1, stats::sd)
  s[!is.finite(s) | s == 0] <- NA
  Z <- (X - mu) / s; Z[!is.finite(Z)] <- NA; Z
}

# Multi-omic features for one domain: samples x features (reversed), complete cases
multiomic_features <- function(d, dat, features) {
  med <- CORE_MEDIATORS %>% filter(domain == d)
  cols <- list(); mult <- c()
  for (i in seq_len(nrow(med))) {
    r <- med[i, ]; m <- if (r$reverse) -1 else 1
    if (!is.na(r$protein) && !is.null(dat$olink) && r$protein %in% names(dat$olink)) {
      v <- suppressWarnings(as.numeric(dat$olink[[r$protein]][match(dat$ids, dat$olink[[dat$cfg$olink_id_col]])]))
      cols[[paste0(r$protein, "_prot")]] <- v; mult[paste0(r$protein, "_prot")] <- m
    }
    key <- if (dat$key == "entrez") r$entrez else r$gene
    if (key %in% rownames(dat$expr)) { cols[[paste0(r$gene, "_gene")]] <- as.numeric(dat$expr[key, dat$ids]); mult[paste0(r$gene, "_gene")] <- m }
  }
  X <- as.data.frame(cols, check.names = FALSE); rownames(X) <- dat$ids
  if (!is.null(features)) {
    absent <- setdiff(features, names(X))
    if (length(absent)) { warning(d, ": features unavailable, multi-omic score not computed: ", paste(absent, collapse = ", ")); return(NULL) }
    X <- X[, features, drop = FALSE]; mult <- mult[features]
  }
  X <- X[stats::complete.cases(X), , drop = FALSE]
  Z <- sweep(scale(as.matrix(X)), 2, mult[colnames(X)], `*`)                # z-score within cohort, then reverse
  Z
}

# PC1 with the majority-positive sign rule ("mean" = orient PC1 to correlate positively
# with the mean of the standardized features; used only for reporting)
pc1 <- function(Z, orient = c("majority", "mean")) {
  orient <- match.arg(orient)
  p <- stats::prcomp(Z, center = FALSE, scale. = FALSE)
  w <- p$rotation[, 1]; s <- p$x[, 1]
  flip <- if (orient == "majority") sum(w < 0) > sum(w > 0) else isTRUE(stats::cor(s, rowMeans(Z)) < 0)
  if (flip) { w <- -w; s <- -s }
  list(weights = w, scores = s, eig = p$sdev^2)
}
tucker <- function(a, b) { k <- intersect(names(a), names(b)); sum(a[k] * b[k]) / sqrt(sum(a[k]^2) * sum(b[k]^2)) }
gene_key <- function(tbl, dat) if (dat$key == "entrez") tbl$entrez else tbl$gene


# =============================================================================
# SECTION 7. DERIVATION FROM THE DISCOVERY COHORT (DERIVE = TRUE)
# =============================================================================
W_MULTI <- split(PUBLISHED_MULTIOMIC_WEIGHTS, PUBLISHED_MULTIOMIC_WEIGHTS$domain) %>%
  map(~ stats::setNames(.x$weight, .x$feature))
SURR <- PUBLISHED_SURROGATE
oof_scores <- NULL

if (DERIVE) {
  message("\n=== Deriving locked quantities from ", DISCOVERY, " ===")
  dd <- data[[DISCOVERY]]
  der_dir <- file.path(OUTPUT_DIR, "derivation"); dir.create(der_dir, showWarnings = FALSE)
  shared_with <- function(d) {                       # features available in every cohort with Olink data
    sets <- map(Filter(function(x) !is.null(x$olink), data), function(x) colnames(multiomic_features(d, x, NULL)))
    Reduce(intersect, sets)
  }
  comp <- list()

  # 7a. Multi-omic PC1 weights
  disc_scores <- list(); W_new <- list(); disc_parts <- list(); disc_summ <- list(); disc_feat <- list()
  for (d in DOMAINS) {
    feats <- shared_with(d)
    Z <- multiomic_features(d, dd, feats)
    f <- pc1(Z, "majority")
    W_new[[d]] <- f$weights; disc_scores[[d]] <- f$scores
    disc_summ[[d]] <- tibble(cohort = DISCOVERY, domain = d, n_features = ncol(Z), n_samples = nrow(Z),
                             pct_variance_PC1 = 100 * f$eig[1] / sum(f$eig))
    disc_feat[[d]] <- tibble(cohort = DISCOVERY, domain = d, feature = colnames(Z), weight = unname(f$weights[colnames(Z)]),
                             r_feature_score = as.numeric(stats::cor(Z, f$scores)))
    # protein and transcriptome components of the multi-omic score (they sum to the score)
    prot <- grep("_prot$", colnames(Z), value = TRUE); rna <- setdiff(colnames(Z), prot)
    disc_parts[[d]] <- list(
      protein = if (length(prot)) as.numeric(Z[, prot, drop = FALSE] %*% f$weights[prot]) else rep(NA_real_, nrow(Z)),
      rna     = as.numeric(Z[, rna, drop = FALSE] %*% f$weights[rna]))
    pub <- W_MULTI[[d]]
    comp[[paste("multiomic", d)]] <- tibble(quantity = "multi-omic PC1 weights", domain = d,
      n_compared = length(intersect(names(pub), names(f$weights))),
      max_abs_difference = max(abs(pub[names(f$weights)] - f$weights), na.rm = TRUE),
      same_features = setequal(names(pub), names(f$weights)),
      note = sprintf("PC1 variance explained %.1f%%; sign rule based on the feature mean agrees: %s", 100 * f$eig[1] / sum(f$eig),
                     sign(stats::cor(f$scores, rowMeans(Z))) > 0))
  }
  W_MULTI <- W_new
  write_csv(imap_dfr(W_new, ~ tibble(domain = .y, feature = names(.x), weight = unname(.x))), file.path(der_dir, "multiomic_PC1_weights.csv"))
  write_csv(bind_rows(disc_summ), file.path(der_dir, "multiomic_PC1_variance_explained.csv"))
  write_csv(bind_rows(disc_feat), file.path(der_dir, "multiomic_PC1_feature_score_correlations.csv"))

  # 7b. LASSO transcriptomic projection (repeated 10-fold CV, median lambda.1se, refit on all)
  if ("surrogate" %in% SCORES) {
    surr_new <- list(); perf <- list(); oof <- list()
    for (d in DOMAINS) {
      y <- disc_scores[[d]]; pid <- names(y)
      keys <- gene_key(CORE_MEDIATORS, dd); keys <- keys[keys %in% rownames(dd$expr)]
      Xg <- dd$expr[keys, pid, drop = FALSE]
      sdv <- apply(Xg, 1, stats::sd); Xg <- Xg[is.finite(sdv) & sdv > 0, , drop = FALSE]
      X <- t((Xg - rowMeans(Xg)) / apply(Xg, 1, stats::sd))
      cvs <- map(seq_len(N_CV_REPEATS), function(r) { set.seed(SEED_BASE + r)
        glmnet::cv.glmnet(X, y, alpha = 1, nfolds = N_FOLDS, standardize = FALSE, keep = TRUE) })
      lam <- stats::median(map_dbl(cvs, "lambda.1se"))
      full <- glmnet::glmnet(X, y, alpha = 1, standardize = FALSE)
      b <- as.matrix(stats::coef(full, s = lam)); cf <- b[-1, 1]; cf <- cf[cf != 0]
      gsym <- if (dd$key == "entrez") CORE_MEDIATORS$gene[match(names(cf), CORE_MEDIATORS$entrez)] else names(cf)
      surr_new[[d]] <- tibble(domain = d, gene = gsym, entrez = CORE_MEDIATORS$entrez[match(gsym, CORE_MEDIATORS$gene)], coefficient = unname(cf))
      # Out-of-fold (OOF) predictions: for each of the 10 CV repeats, take the held-out
      # prediction at the lambda closest to the selected (median) lambda. The TOR projection
      # scores are those of repeat 1 (OOF_METHOD = "first_repeat") or the per-patient mean
      # across repeats ("average_10"). The range of OOF r across the 10 repeats is reported
      # in both cases to show sensitivity to fold assignment.
      oof_by_rep <- sapply(cvs, function(cv) cv$fit.preval[, which.min(abs(cv$lambda - lam))])
      oof_pred <- if (OOF_METHOD == "first_repeat") oof_by_rep[, 1] else rowMeans(oof_by_rep)
      r_by_rep <- apply(oof_by_rep, 2, function(p) stats::cor(p, y))
      oof[[d]] <- tibble(sample_id = pid, domain = d, score_raw = oof_pred)
      perf[[d]] <- tibble(domain = d, lambda = lam, intercept = b[1, 1], n_genes_selected = length(cf),
                          oof_pearson_r = stats::cor(oof_pred, y), oof_R2 = 1 - sum((y - oof_pred)^2) / sum((y - mean(y))^2),
                          oof_calibration_slope = unname(stats::coef(stats::lm(y ~ oof_pred))[2]),
                          oof_r_single_repeat_min = min(r_by_rep), oof_r_single_repeat_max = max(r_by_rep),
                          oof_r_first_repeat = r_by_rep[1],
                          oof_r_with_rna_component = stats::cor(oof_pred, disc_parts[[d]]$rna),
                          oof_r_with_protein_component = stats::cor(oof_pred, disc_parts[[d]]$protein),
                          score_r_with_protein_component = stats::cor(y, disc_parts[[d]]$protein))
      pub <- PUBLISHED_SURROGATE %>% filter(domain == d)
      j <- full_join(pub %>% dplyr::select(gene, pub = coefficient), surr_new[[d]] %>% dplyr::select(gene, new = coefficient), by = "gene")
      comp[[paste("surrogate", d)]] <- tibble(quantity = "transcriptomic projection coefficients", domain = d, n_compared = nrow(j),
        max_abs_difference = max(abs(coalesce(j$pub, 0) - coalesce(j$new, 0))), same_features = all(!is.na(j$pub) & !is.na(j$new)),
        note = sprintf("lambda %.4f (published %.4f); intercept %.2e; out-of-fold r = %.3f", lam, PUBLISHED_SURROGATE_LAMBDA[[d]], b[1, 1], stats::cor(oof_pred, y)))
    }
    SURR <- bind_rows(surr_new)
    oof_scores <- bind_rows(oof)
    write_csv(SURR, file.path(der_dir, "surrogate_coefficients.csv"))
    write_csv(bind_rows(perf), file.path(der_dir, "surrogate_crossvalidated_performance.csv"))
  }

  comp_tbl <- bind_rows(comp)
  write_csv(comp_tbl, file.path(der_dir, "derivation_vs_published.csv"))
  print(as.data.frame(comp_tbl %>% mutate(max_abs_difference = signif(max_abs_difference, 3))))
}


# =============================================================================
# SECTION 8. SCORE EACH COHORT
# =============================================================================
coverage_all <- list(); transfer <- list(); ext_perf <- list()
for (nm in names(data)) {
  dat <- data[[nm]]; message("\n=== Scoring ", nm, " ===")
  out <- list(); cov <- list(); mo_parts <- list()

  # 8a. Multi-omic score
  if ("multiomic" %in% SCORES && !is.null(dat$olink)) {
    for (d in DOMAINS) {
      w <- W_MULTI[[d]]
      Z <- multiomic_features(d, dat, names(w)); if (is.null(Z)) next
      s <- as.numeric(Z %*% w[colnames(Z)])
      out[[paste("mo", d)]] <- tibble(sample_id = rownames(Z), domain = d, method = "multiomic_PC1", score_raw = s, in_ref = TRUE)
      cov[[paste("mo", d)]] <- tibble(method = "multiomic_PC1", domain = d, n_features = length(w), n_found = ncol(Z), n_samples = nrow(Z))
      # protein and transcriptome components of the score (discovery-cohort weights; they sum to the score)
      prot <- grep("_prot$", colnames(Z), value = TRUE); rna <- setdiff(colnames(Z), prot)
      mo_parts[[d]] <- tibble(sample_id = rownames(Z),
        protein = if (length(prot)) as.numeric(Z[, prot, drop = FALSE] %*% w[prot]) else NA_real_,
        rna     = as.numeric(Z[, rna, drop = FALSE] %*% w[rna]))
      # transfer checks for cohorts other than the discovery cohort
      if (!(DERIVE && nm == DISCOVERY)) {
        own <- pc1(Z, "majority"); tv <- sum(apply(Z, 2, stats::var)); n <- nrow(Z)
        transfer[[paste(nm, d)]] <- list(
          variance = tibble(cohort = nm, domain = d, n_features = ncol(Z), n_samples = n,
            pct_variance_own_PC1 = 100 * stats::var(own$scores) / tv, pct_variance_discovery_axis = 100 * stats::var(s) / tv,
            ratio_discovery_axis_vs_own = stats::var(s) / stats::var(own$scores), tucker_phi = tucker(w, own$weights)),
          eigen = tibble(cohort = nm, domain = d, n_samples = n, lambda1 = own$eig[1], lambda2 = own$eig[2],
            ratio = own$eig[1] / own$eig[2], gap = own$eig[1] - own$eig[2], north_error = own$eig[1] * sqrt(2 / n),
            pc1_separable = (own$eig[1] - own$eig[2]) > own$eig[1] * sqrt(2 / n)),
          features = tibble(cohort = nm, domain = d, feature = colnames(Z), weight = unname(w[colnames(Z)]),
            r_feature_score = as.numeric(stats::cor(Z, s)), same_direction = sign(as.numeric(stats::cor(Z, s))) == sign(unname(w[colnames(Z)])),
            own_cohort_weight = unname(own$weights[colnames(Z)])))
      }
    }
  }

  # 8b. Transcriptomic projection
  if ("surrogate" %in% SCORES) {
    for (d in DOMAINS) {
      m <- SURR %>% filter(domain == d); k <- gene_key(m, dat); found <- k %in% rownames(dat$expr)
      Z <- z_against(dat$expr[k[found], dat$ids, drop = FALSE], dat$ref); Z[is.na(Z)] <- 0
      s <- as.numeric(crossprod(Z, m$coefficient[found]))
      out[[paste("su", d)]] <- tibble(sample_id = dat$ids, domain = d, method = "rna_surrogate_lasso14_1se", score_raw = s, in_ref = dat$ids %in% dat$ref)
      cov[[paste("su", d)]] <- tibble(method = "rna_surrogate_lasso14_1se", domain = d, n_features = nrow(m), n_found = sum(found),
                                      pct_weight_found = 100 * sum(abs(m$coefficient[found])) / sum(abs(m$coefficient)))
    }
    if (DERIVE && nm == DISCOVERY && !is.null(oof_scores))
      out[["oof"]] <- oof_scores %>% mutate(method = "rna_surrogate_lasso14_1se_OOF", in_ref = TRUE)
  }

  # 8c. External performance of the transcriptomic projection against the multi-omic score
  #     (cohorts with Olink data other than the discovery cohort; locked predictions,
  #     samples with both scores). Same metrics as the out-of-fold performance in 7b.
  if (all(c("multiomic", "surrogate") %in% SCORES) && !is.null(dat$olink) && !(DERIVE && nm == DISCOVERY)) {
    for (d in DOMAINS) {
      mo <- out[[paste("mo", d)]]; su <- out[[paste("su", d)]]
      if (is.null(mo) || is.null(su)) next
      j <- mo %>% dplyr::select(sample_id, y = score_raw) %>%
        inner_join(su %>% dplyr::select(sample_id, p = score_raw), by = "sample_id") %>%
        inner_join(mo_parts[[d]], by = "sample_id") %>% filter(is.finite(y), is.finite(p))
      ext_perf[[paste(nm, d)]] <- tibble(cohort = nm, domain = d, n = nrow(j),
        pearson_r = stats::cor(j$p, j$y), R2 = 1 - sum((j$y - j$p)^2) / sum((j$y - mean(j$y))^2),
        calibration_slope = unname(stats::coef(stats::lm(y ~ p, data = j))[2]),
        r_with_rna_component = stats::cor(j$p, j$rna),
        r_with_protein_component = if (all(is.finite(j$protein))) stats::cor(j$p, j$protein) else NA_real_,
        score_r_with_protein_component = if (all(is.finite(j$protein))) stats::cor(j$y, j$protein) else NA_real_)
    }
  }

  if (!length(out)) { message("  nothing to score"); next }

  long <- bind_rows(out) %>%
    group_by(method, domain) %>%
    mutate(score_perSD = (score_raw - mean(score_raw[in_ref])) / stats::sd(score_raw[in_ref])) %>%
    ungroup() %>%
    mutate(cohort = nm, reference_sample = in_ref, method = factor(method, levels = METHOD_INFO$method),
           domain = factor(domain, levels = DOMAINS)) %>%
    arrange(method, domain) %>% mutate(method = as.character(method), domain = as.character(domain)) %>%
    dplyr::select(cohort, sample_id, reference_sample, method, domain, score_raw, score_perSD)
  if (length(dat$cfg$keep_cols)) {
    mk <- dat$meta[, c(dat$cfg$meta_id_col, dat$cfg$keep_cols), drop = FALSE]; names(mk)[1] <- "sample_id"
    long <- long %>% left_join(mk, by = "sample_id") %>% relocate(all_of(dat$cfg$keep_cols), .after = sample_id)
  }
  wide <- long %>% mutate(col = paste(method, domain, sep = "__")) %>%
    dplyr::select(-reference_sample, -method, -domain, -score_raw) %>%
    pivot_wider(names_from = col, values_from = score_perSD)
  cdir <- file.path(OUTPUT_DIR, nm); dir.create(cdir, showWarnings = FALSE)
  write_csv(long, file.path(cdir, paste0(nm, "_DRT_scores_long.csv")))
  write_csv(wide, file.path(cdir, paste0(nm, "_DRT_scores_wide_perSD.csv")))
  cov_tbl <- bind_rows(cov) %>% mutate(cohort = nm, .before = 1)
  write_csv(cov_tbl, file.path(cdir, "gene_coverage.csv"))
  coverage_all[[nm]] <- cov_tbl
  message(sprintf("  %d samples; methods: %s", length(unique(long$sample_id)), paste(unique(long$method), collapse = ", ")))
}
print(as.data.frame(bind_rows(coverage_all)))

if (length(transfer)) {
  tdir <- file.path(OUTPUT_DIR, "transfer_checks"); dir.create(tdir, showWarnings = FALSE)
  write_csv(map_dfr(transfer, "variance"), file.path(tdir, "A_variance_captured_and_tucker_phi.csv"))
  write_csv(map_dfr(transfer, "eigen"),    file.path(tdir, "B_PC1_PC2_separation_North.csv"))
  write_csv(map_dfr(transfer, "features"), file.path(tdir, "C_feature_score_correlations.csv"))
  print(as.data.frame(map_dfr(transfer, "variance") %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
}
if (length(ext_perf)) {
  tdir <- file.path(OUTPUT_DIR, "transfer_checks"); dir.create(tdir, showWarnings = FALSE)
  write_csv(bind_rows(ext_perf), file.path(tdir, "D_transcriptomic_projection_external_performance.csv"))
  print(as.data.frame(bind_rows(ext_perf) %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
}

write_csv(METHOD_INFO %>% mutate(score_raw = "score as computed",
                                 score_perSD = "standardized to mean 0, SD 1 within the cohort's reference samples (multi-omic: all scored samples)"),
          file.path(OUTPUT_DIR, "method_dictionary.csv"))
format_setting <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.null(names(x))) return(paste(unlist(x), collapse = "; "))
  paste(names(x), unlist(x), sep = " = ", collapse = "; ")
}
write_csv(tibble(setting = c("DERIVE", "DISCOVERY", "SCORES", "run_time"),
                 value = c(DERIVE, DISCOVERY, paste(SCORES, collapse = ", "), format(Sys.time()))) %>%
            bind_rows(imap_dfr(COHORTS, function(cfg, nm) tibble(setting = paste0(nm, ": ", names(cfg)),
                       value = vapply(cfg, format_setting, character(1))))),
          file.path(OUTPUT_DIR, "run_settings.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
message("\nDone. Outputs in ", OUTPUT_DIR)
close_log()
