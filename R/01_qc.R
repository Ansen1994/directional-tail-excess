## =============================================================================
## 01_qc.R - per-sample quality control and filtering
##
## Input : data/raw_sc/<SAMPLE>/{count_matrix_sparse.mtx, count_matrix_genes.tsv,
##         count_matrix_barcodes.tsv, metadata.csv}
##         data/raw_sc is a junction to the original GEO download (read only)
## Output: results/tables/qc_per_sample.csv          per-sample counts + thresholds
##         results/tables/qc_cell_level.csv.gz       per-cell metrics + pass flag
##         results/figures/fig01_qc_overview.pdf/png
##         data/processed/01_counts_qc.rds           filtered merged count matrix
##
## Thresholding strategy (two constraints, the stricter one wins)
##   hard defaults (QC_MIN_FEATURES / QC_MAX_FEATURES / QC_MAX_MITO)
##   + per-sample adaptive MAD = median +/- QC_MAD_N * MAD
##
##   The nFeature UPPER bound is deliberately fixed and does not use MAD:
##   measured on this dataset, median + 3*MAD lands close to twice the median in
##   low-complexity samples (e.g. CID3946, median 421) and would remove nearly
##   half of the cells. nFeature is right-skewed, so MAD is a poor upper bound.
##   High-gene-count cells are left to doublet detection and downstream checks.
##   The nFeature LOWER bound and the mitochondrial upper bound do use MAD
##   (low-quality and dying cells), but a mitochondrial MAD bound below 10% is
##   rejected to avoid over-trimming tumour cells in intrinsically clean samples.
##   Effective per-sample thresholds are written to the table for review.
##
## NOTE: keep this file pure ASCII. R under LC_CTYPE=C on this machine fails to
## parse some multibyte string literals in source files.
## =============================================================================

source("config/00_config.R")
suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(ggplot2); library(patchwork)
})

set.seed(WB_SEED)
t0 <- Sys.time()

## -------------------------------------------------------------- sample scan
dirs <- list.dirs(RAW_DIR, recursive = FALSE, full.names = TRUE)
dirs <- dirs[file.exists(file.path(dirs, "count_matrix_sparse.mtx")) &
             file.exists(file.path(dirs, "metadata.csv"))]
dirs <- sort(dirs)
if (WB_SC_SAMPLE_LIMIT > 0) dirs <- head(dirs, WB_SC_SAMPLE_LIMIT)
if (!length(dirs)) stop("no sample folders found under ", RAW_DIR)
samples <- basename(dirs)
msg("%d samples: %s ... %s", length(dirs), samples[1], tail(samples, 1))
if (WB_SC_CELLS_PER_SAMPLE > 0)
  msg("[SMOKE TEST] at most %d cells per sample", WB_SC_CELLS_PER_SAMPLE)

mad_bounds <- function(x, n = QC_MAD_N) {
  m <- stats::median(x, na.rm = TRUE)
  s <- stats::mad(x, na.rm = TRUE)
  c(low = m - n * s, high = m + n * s)
}

mats <- list(); cell_qc <- list(); per_sample <- list()

for (d in dirs) {
  sid <- basename(d)
  gz  <- readLines(file.path(d, "count_matrix_genes.tsv"), warn = FALSE)
  bz  <- readLines(file.path(d, "count_matrix_barcodes.tsv"), warn = FALSE)
  m   <- Matrix::readMM(file.path(d, "count_matrix_sparse.mtx"))
  m   <- as(m, "CsparseMatrix")
  if (nrow(m) != length(gz))
    stop(sid, " gene dimension mismatch: ", nrow(m), " vs ", length(gz))
  if (ncol(m) != length(bz))
    stop(sid, " cell dimension mismatch: ", ncol(m), " vs ", length(bz))
  rownames(m) <- make.unique(gz); colnames(m) <- bz

  if (WB_SC_CELLS_PER_SAMPLE > 0 && ncol(m) > WB_SC_CELLS_PER_SAMPLE) {
    keep <- sort(sample.int(ncol(m), WB_SC_CELLS_PER_SAMPLE))
    m <- m[, keep, drop = FALSE]; bz <- bz[keep]
  }

  ## QC metrics on raw counts, before any normalisation
  n_count <- Matrix::colSums(m)
  n_feat  <- Matrix::colSums(m > 0)
  pat_pct <- function(pat) {
    i <- grep(pat, rownames(m), ignore.case = TRUE)
    if (!length(i)) return(rep(0, ncol(m)))
    as.numeric(Matrix::colSums(m[i, , drop = FALSE])) / pmax(n_count, 1) * 100
  }
  pct_mito <- pat_pct(GENE_PAT_MITO)
  pct_ribo <- pat_pct(GENE_PAT_RIBO)
  pct_hb   <- pat_pct(GENE_PAT_HB)

  ## effective thresholds
  fb <- mad_bounds(n_feat); mb <- mad_bounds(pct_mito)
  min_feat <- QC_MIN_FEATURES
  max_feat <- QC_MAX_FEATURES          # fixed on purpose, see header
  max_mito <- QC_MAX_MITO
  if (QC_USE_MAD) {
    min_feat <- max(min_feat, floor(fb[["low"]]))
    mh <- mb[["high"]]
    if (mh >= 10) max_mito <- min(max_mito, mh)
  }
  keep <- n_feat >= min_feat & n_feat <= max_feat &
          pct_mito <= max_mito & n_count >= QC_MIN_COUNTS

  ## attribution: how many cells each criterion alone would remove (overlapping)
  fail_low_feat   <- sum(n_feat < min_feat)
  fail_high_feat  <- sum(n_feat > max_feat)
  fail_mito       <- sum(pct_mito > max_mito)
  fail_low_count  <- sum(n_count < QC_MIN_COUNTS)

  per_sample[[sid]] <- data.frame(
    sample = sid, n_cells_raw = ncol(m), n_cells_keep = sum(keep),
    keep_pct = round(100 * sum(keep) / ncol(m), 2),
    min_features = min_feat, max_features = max_feat, max_mito = max_mito,
    med_features = round(stats::median(n_feat), 1),
    med_counts = round(stats::median(n_count), 1),
    med_mito = round(stats::median(pct_mito), 2),
    fail_low_feat = fail_low_feat, fail_high_feat = fail_high_feat,
    fail_mito = fail_mito, fail_low_count = fail_low_count,
    pct_feat_over_6000 = round(100 * mean(n_feat > 6000), 2),
    stringsAsFactors = FALSE)

  ## original curated labels (for downstream agreement checks, never for filtering)
  md <- tryCatch(read.csv(file.path(d, "metadata.csv"), row.names = 1,
                          stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) NULL)
  pick <- function(col) {
    if (is.null(md) || !col %in% colnames(md))
      return(rep(NA_character_, length(bz)))
    as.character(md[match(bz, rownames(md)), col])
  }
  ref_subtype <- pick("subtype")
  ref_major   <- pick("celltype_major")
  ref_minor   <- pick("celltype_minor")
  ref_subset  <- pick("celltype_subset")
  if (all(is.na(ref_major)))
    msg("  %-9s WARNING: curated labels not matched, agreement check will skip", sid)

  cell_qc[[sid]] <- data.frame(
    cell_id = bz, sample = sid, nCount_RNA = n_count, nFeature_RNA = n_feat,
    percent_mito = pct_mito, percent_ribo = pct_ribo, percent_hb = pct_hb,
    ref_subtype = ref_subtype, ref_major = ref_major,
    ref_minor = ref_minor, ref_subset = ref_subset,
    pass_qc = keep, stringsAsFactors = FALSE)

  mats[[sid]] <- m[, keep, drop = FALSE]
  msg("  %-9s raw %5d -> keep %5d (%.1f%%) | nFeature[%d,%d] mito<=%.1f",
      sid, ncol(m), sum(keep), 100 * sum(keep) / ncol(m), min_feat, max_feat, max_mito)
}

qc_tab <- do.call(rbind, per_sample); rownames(qc_tab) <- NULL
cell_df <- do.call(rbind, cell_qc);  rownames(cell_df) <- NULL
write_csv_safe(qc_tab, file.path(TAB_DIR, "qc_per_sample.csv"))
con <- gzfile(file.path(TAB_DIR, "qc_cell_level.csv.gz"), "wt")
utils::write.csv(cell_df, con, row.names = FALSE); close(con)

## -------------------------------------------------------- merge and save
counts <- do.call(cbind, mats)
rm(mats); invisible(gc(FALSE))
counts <- counts[Matrix::rowSums(counts > 0) >= QC_MIN_CELLS_GENE, ]
msg("filtered matrix: %d genes x %d cells, nnz=%.2e",
    nrow(counts), ncol(counts), length(counts@x))

meta <- cell_df[cell_df$pass_qc, c("cell_id", "sample", "nCount_RNA", "nFeature_RNA",
                                   "percent_mito", "percent_ribo", "percent_hb",
                                   "ref_subtype", "ref_major", "ref_minor", "ref_subset")]
rownames(meta) <- meta$cell_id
stopifnot(identical(colnames(counts), rownames(meta)))
saveRDS(list(counts = counts, meta = meta, per_sample = qc_tab,
             params = list(min_features = QC_MIN_FEATURES, max_features = QC_MAX_FEATURES,
                           max_mito = QC_MAX_MITO, min_counts = QC_MIN_COUNTS,
                           use_mad = QC_USE_MAD, mad_n = QC_MAD_N, seed = WB_SEED)),
        file.path(PROC_DIR, "01_counts_qc.rds"), compress = "gzip")

## ---------------------------------------------------------------- figures
p_metric <- function(var, ylab, hline = NULL, logy = FALSE) {
  ggplot(cell_df, aes(x = sample, y = .data[[var]])) +
    geom_violin(scale = "width", fill = "#7BA7BC", colour = NA, alpha = .55) +
    geom_boxplot(width = .12, outlier.shape = NA, fill = "white", linewidth = .25) +
    { if (!is.null(hline)) geom_hline(yintercept = hline, linetype = "dashed",
                                      colour = "#C0504D", linewidth = .35) } +
    labs(x = NULL, y = ylab) +
    { if (logy) scale_y_log10() } +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 6),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
}

p1 <- p_metric("nFeature_RNA", "Genes per cell", logy = TRUE)
p2 <- p_metric("nCount_RNA", "UMIs per cell", logy = TRUE)
p3 <- p_metric("percent_mito", "Mitochondrial %", hline = QC_MAX_MITO)

p4 <- ggplot(cell_df, aes(x = nCount_RNA, y = nFeature_RNA, colour = percent_mito)) +
  geom_point(size = .15, alpha = .35) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_gradientn(colours = c("#4A6FA5", "#9CAF88", "#E8A87C", "#C0504D")) +
  labs(x = "UMIs per cell (log10)", y = "Genes per cell (log10)", colour = "Mito %") +
  theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank())

p5 <- ggplot(qc_tab, aes(x = sample)) +
  geom_col(aes(y = n_cells_raw), fill = "#D9D9D9", width = .7) +
  geom_col(aes(y = n_cells_keep), fill = "#7BA7BC", width = .7) +
  labs(x = NULL, y = "Cells", title = "Before (grey) / after (blue) QC") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 6),
        plot.title = element_text(size = 8), panel.grid.minor = element_blank())

p6 <- ggplot(cell_df, aes(x = percent_mito,
                          colour = ifelse(pass_qc, "Passed", "Failed"),
                          fill = ifelse(pass_qc, "Passed", "Failed"))) +
  geom_density(alpha = .35, linewidth = .3) +
  coord_cartesian(xlim = c(0, max(QC_MAX_MITO * 1.25, 15))) +
  scale_colour_manual(values = c(Passed = "#7BA7BC", Failed = "#C0504D")) +
  scale_fill_manual(values = c(Passed = "#7BA7BC", Failed = "#C0504D")) +
  labs(x = "Mitochondrial %", y = "Density", colour = NULL, fill = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

fig <- (p1 | p2 | p3) / (p5 | p4 | p6) +
  plot_annotation(tag_levels = "A",
                  caption = sprintf("GSE176078 QC: %d samples, %d/%d cells retained (%.1f%%)",
                                    nrow(qc_tab), sum(qc_tab$n_cells_keep),
                                    sum(qc_tab$n_cells_raw),
                                    100 * sum(qc_tab$n_cells_keep) / sum(qc_tab$n_cells_raw)))
ggsave(file.path(FIG_DIR, "fig01_qc_overview.pdf"), fig, width = 13, height = 9)
ggsave(file.path(FIG_DIR, "fig01_qc_overview.png"), fig, width = 13, height = 9, dpi = 300)

## ---------------------------------------------------------------- summary
tot_raw <- sum(qc_tab$n_cells_raw); tot_keep <- sum(qc_tab$n_cells_keep)
con <- file(file.path(LOG_DIR, "01_qc_summary.txt"), "w")
writeLines(c(
  sprintf("run at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("samples: %d", nrow(qc_tab)),
  sprintf("cells: %d -> %d (%.2f%% retained)", tot_raw, tot_keep,
          100 * tot_keep / tot_raw),
  sprintf("genes after filtering (>=%d cells): %d", QC_MIN_CELLS_GENE, nrow(counts)),
  sprintf("failure attribution (overlapping): low nFeature %d / high nFeature %d / high mito %d / low UMI %d",
          sum(qc_tab$fail_low_feat), sum(qc_tab$fail_high_feat),
          sum(qc_tab$fail_mito), sum(qc_tab$fail_low_count)),
  sprintf("params: min_features=%d max_features=%d max_mito=%.1f min_counts=%d MAD=%s(n=%.1f)",
          QC_MIN_FEATURES, QC_MAX_FEATURES, QC_MAX_MITO, QC_MIN_COUNTS,
          QC_USE_MAD, QC_MAD_N),
  sprintf("elapsed: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
), con); close(con)

writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo_01_qc.txt"))
msg("done: %d -> %d cells in %.1f min", tot_raw, tot_keep,
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
