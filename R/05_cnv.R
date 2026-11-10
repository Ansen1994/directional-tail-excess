## =============================================================================
## 05_cnv.R -- per-cell copy-number score, malignant vs normal epithelium
##
## Why not copykat: it writes two ~200 MB bin-by-cell dumps plus a matching
## heatmap per sample. One 4k-cell sample needed ~50 min and >1 GB on this
## machine, i.e. hours for the cohort, for a call we only need at
## "malignant vs normal" resolution.
##
## Method (arm-level, inferCNV-style):
##   1. genes ordered by genomic position (chr, then start)
##   2. per-cell log normalised expression
##   3. baseline = per-gene mean of NORMAL EPITHELIAL cells from all OTHER
##      samples (leave-one-sample-out). This is the critical choice: an
##      immune-cell baseline leaves a lineage expression programme that is
##      unevenly distributed across the genome, which swamps the CNV signal.
##      Measured on CID44991, an immune baseline gave AUC 0.34 for cancer vs
##      normal epithelium; the epithelial LOSO baseline gives clean separation
##      (96% of cancer cells above the reference median, 0% of normal cells).
##   4. residual, clipped at +/-0.35, smoothed with a 101-gene running mean
##   5. score = SD of the cell's smoothed profile after centring on its own
##      median (without centring the score mostly measures expression breadth)
##   6. threshold = 95th percentile of the baseline cells' own scores
##
## Resolution caveat: arm-level. Detects whole-arm gains/losses, not focal
## amplifications. Not suitable for subclone reconstruction.
##
## Normal epithelial cells of sample X are excluded from X's baseline, so they
## remain a fair test set for specificity.
## =============================================================================

source("config/00_config.R")
suppressPackageStartupMessages({
  library(Matrix)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(patchwork)
})

set.seed(WB_SEED)
t0 <- Sys.time()

## ------------------------------------------------------- gene position table
MAP_FILE <- file.path(PROC_DIR, "gene_position.csv")
if (file.exists(MAP_FILE)) {
  gmap <- read.csv(MAP_FILE, stringsAsFactors = FALSE)
  msg("gene map loaded from cache: %d genes", nrow(gmap))
} else {
  obj0 <- readRDS(file.path(PROC_DIR, "04_doublets.rds"))
  genes <- rownames(obj0$counts)
  rm(obj0); invisible(gc(FALSE))
  msg("building gene position map for %d genes ...", length(genes))
  chr <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL",
                               columns = c("ENTREZID", "CHR"))
  chr <- chr[!is.na(chr$CHR), ]
  loc <- suppressWarnings(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(chr$ENTREZID), keytype = "ENTREZID",
    columns = "CHRLOC"))
  loc <- loc[!is.na(loc$CHRLOC), ]
  loc$start <- abs(loc$CHRLOC)
  loc <- loc[!duplicated(loc$ENTREZID), c("ENTREZID", "start")]
  gmap <- merge(chr, loc, by = "ENTREZID")
  gmap <- gmap[!duplicated(gmap$SYMBOL), c("SYMBOL", "CHR", "start")]
  gmap <- gmap[gmap$CHR %in% c(as.character(1:22), "X", "Y"), ]
  gmap$chr_num <- suppressWarnings(ifelse(gmap$CHR == "X", 23,
                                   ifelse(gmap$CHR == "Y", 24, as.integer(gmap$CHR))))
  gmap <- gmap[order(gmap$chr_num, gmap$start), ]
  write_csv_safe(gmap, MAP_FILE)
  msg("gene map built: %d genes with position", nrow(gmap))
}

## ------------------------------------------------------------------- prepare
obj <- readRDS(file.path(PROC_DIR, "04_doublets.rds"))
counts <- obj$counts
meta   <- obj$meta
ok <- meta$pass_doublet

is_norm <- ok & meta$ref_major == "Normal Epithelial"
is_epi  <- ok & meta$ref_major %in% c("Cancer Epithelial", "Normal Epithelial")
norm_by_sample <- split(which(is_norm), meta$sample[is_norm])
msg("normal epithelial cells available in %d samples (total %d)",
    length(norm_by_sample), sum(is_norm))

gm <- gmap[gmap$SYMBOL %in% rownames(counts), ]
msg("genes with genomic position: %d", nrow(gm))

samples <- unique(meta$sample[ok])
if (WB_SC_SAMPLE_LIMIT > 0) samples <- head(samples, WB_SC_SAMPLE_LIMIT)

calls_list <- vector("list", length(samples))
per_sample <- vector("list", length(samples))

for (i in seq_along(samples)) {
  sid <- samples[i]
  epi <- which(is_epi & meta$sample == sid)
  ref <- unlist(norm_by_sample[setdiff(names(norm_by_sample), sid)])
  if (length(ref) > CNV_MAX_REF) ref <- sample(ref, CNV_MAX_REF)

  if (length(epi) < CNV_MIN_EPI || length(ref) < CNV_MIN_REF) {
    msg("[%d/%d] %-9s skipped (epi=%d ref=%d)", i, length(samples), sid,
        length(epi), length(ref))
    per_sample[[i]] <- data.frame(sample = sid, n_epi = length(epi),
                                  n_ref = length(ref), threshold = NA_real_,
                                  n_malignant = 0L, pct_malignant = NA_real_,
                                  status = "skipped", stringsAsFactors = FALSE)
    next
  }

  sel <- c(epi, ref)
  sub <- counts[gm$SYMBOL, sel, drop = FALSE]
  lib <- Matrix::colSums(sub)
  nrm <- sub
  nrm@x <- nrm@x / rep(lib, diff(nrm@p)) * 1e4
  lg <- log1p(as.matrix(nrm))
  rm(nrm, sub); invisible(gc(FALSE))

  nref <- length(ref)
  idx_ref <- (length(sel) - nref + 1):length(sel)
  base <- rowMeans(lg[, idx_ref, drop = FALSE])
  resid <- lg - base
  resid[resid >  CNV_CLIP] <-  CNV_CLIP
  resid[resid < -CNV_CLIP] <- -CNV_CLIP

  sm <- apply(resid, 2, function(v) {
    s <- stats::filter(v, rep(1 / CNV_WIN, CNV_WIN), sides = 2)
    s[is.na(s)] <- v[is.na(s)]
    s
  })
  rm(resid, lg); invisible(gc(FALSE))

  score <- apply(sweep(sm, 2, apply(sm, 2, stats::median), "-"), 2, stats::sd)
  rm(sm); invisible(gc(FALSE))

  thr <- stats::quantile(score[idx_ref], CNV_REF_Q, names = FALSE)
  call <- ifelse(score > thr, "malignant", "normal")
  call[idx_ref] <- "reference"

  calls_list[[i]] <- data.frame(cell_id = meta$cell_id[sel], sample = sid,
                                cnv_score = round(score, 5), cnv_call = call,
                                stringsAsFactors = FALSE)
  n_mal <- sum(call == "malignant")
  per_sample[[i]] <- data.frame(
    sample = sid, n_epi = length(epi), n_ref = nref,
    threshold = round(thr, 4), n_malignant = n_mal,
    pct_malignant = round(100 * n_mal / length(epi), 2),
    status = "ok", stringsAsFactors = FALSE)
  msg("[%d/%d] %-9s epi=%4d ref=%4d thr=%.4f malignant=%4d (%.1f%%)",
      i, length(samples), sid, length(epi), nref, thr, n_mal,
      100 * n_mal / length(epi))
}

per_sample <- do.call(rbind, per_sample)
calls <- do.call(rbind, calls_list)
write_csv_safe(per_sample, file.path(TAB_DIR, "cnv_per_sample.csv"))
write_csv_safe(calls, file.path(TAB_DIR, "cnv_cell_calls.csv"))

## ------------------------------------------------------------- accuracy check
cc <- merge(calls, meta[, c("cell_id", "sample", "ref_major", "ref_subtype")],
            by = c("cell_id", "sample"))
cc <- cc[cc$cnv_call != "reference", ]
tb <- table(factor(cc$ref_major, levels = c("Cancer Epithelial", "Normal Epithelial")),
            factor(cc$cnv_call, levels = c("malignant", "normal")))
write_csv_safe(as.data.frame.matrix(tb), file.path(TAB_DIR, "cnv_vs_ref_labels.csv"))
print(tb)

denom_c <- sum(tb["Cancer Epithelial", ])
denom_n <- sum(tb["Normal Epithelial", ])
sens <- if (denom_c > 0) tb["Cancer Epithelial", "malignant"] / denom_c else NA_real_
spec <- if (denom_n > 0) tb["Normal Epithelial", "normal"] / denom_n else NA_real_
msg("sensitivity %.1f%% (n=%d) ; specificity %.1f%% (n=%d)",
    100 * sens, denom_c, 100 * spec, denom_n)
write_csv_safe(data.frame(metric = c("sensitivity", "specificity"),
                          value = round(100 * c(sens, spec), 2),
                          n = c(denom_c, denom_n)),
               file.path(TAB_DIR, "cnv_accuracy.csv"))

## -------------------------------------------------------------------- figure
ps <- per_sample[per_sample$status == "ok", ]
ps$subtype <- sapply(ps$sample, function(s) {
  v <- meta$ref_subtype[meta$sample == s]; names(sort(table(v), decreasing = TRUE))[1]
})
p1 <- ggplot(ps, aes(reorder(sample, pct_malignant), pct_malignant, fill = subtype)) +
  geom_col(width = .72) + coord_flip() +
  scale_fill_manual(values = c(`ER+` = "#7BA7BC", TNBC = "#C97B63", `HER2+` = "#9CAF88")) +
  labs(x = NULL, y = "Malignant among epithelial (%)", fill = NULL,
       title = "Malignant fraction per tumour") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

dd <- data.frame(score = cc$cnv_score, label = cc$ref_major)
p2 <- ggplot(dd, aes(score, fill = label, colour = label)) +
  geom_density(alpha = .3, linewidth = .3) +
  scale_fill_manual(values = c(`Cancer Epithelial` = "#C97B63",
                               `Normal Epithelial` = "#7BA7BC")) +
  scale_colour_manual(values = c(`Cancer Epithelial` = "#C97B63",
                                 `Normal Epithelial` = "#7BA7BC")) +
  labs(x = "CNV score", y = "Density", fill = NULL, colour = NULL,
       title = sprintf("Score separation (sens %.1f%%, spec %.1f%%)",
                       100 * sens, 100 * spec)) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

fig <- p1 | p2
ggsave(file.path(FIG_DIR, "fig10_cnv.png"), fig, width = 11, height = 5, dpi = 250)
ggsave(file.path(FIG_DIR, "fig10_cnv.pdf"), fig, width = 11, height = 5)

writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo_05_cnv.txt"))
msg("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
