## =============================================================================
## 04_doublets.R -- per-sample doublet detection (scDblFinder)
##
## Why per sample: doublets are formed within a single 10x lane, so the expected
## rate depends on that lane's loaded cell number. Pooling all 26 samples first
## would mix heterogeneous loading densities and bias the prior.
##
## Why it matters here: a doublet of epithelial + immune cell mimics an
## "epithelial-immune hybrid state", which is exactly the biology this project
## intends to describe. Undetected doublets would fabricate the headline result.
##
## Cells are NOT removed from the object. A flag column is added so that every
## downstream step can subset explicitly, and the atlas object stays traceable.
## =============================================================================

source("config/00_config.R")
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
  library(ggplot2)
  library(patchwork)
})

set.seed(WB_SEED)
t0 <- Sys.time()

obj <- readRDS(file.path(PROC_DIR, "01_counts_qc.rds"))
counts <- obj$counts
meta   <- obj$meta
stopifnot(identical(colnames(counts), meta$cell_id))

samples <- unique(meta$sample)
msg("loaded %d cells x %d genes from %d samples", ncol(counts), nrow(counts),
    length(samples))

dbl_score <- rep(NA_real_, ncol(counts))
dbl_class <- rep("singlet", ncol(counts))
dbl_notes <- rep("", ncol(counts))
per_sample <- vector("list", length(samples))

for (i in seq_along(samples)) {
  sid <- samples[i]
  idx <- which(meta$sample == sid)
  n <- length(idx)
  ## 10x Chromium: roughly 0.8% doublets per 1000 recovered cells
  dbr <- min(DBL_RATE_PER_1K * n / 1000, DBL_MAX_RATE)
  msg("[%d/%d] %-9s n=%5d prior_dbr=%.4f", i, length(samples), sid, n, dbr)

  sce <- SingleCellExperiment(list(counts = counts[, idx, drop = FALSE]))
  res <- tryCatch(
    scDblFinder(sce, dbr = dbr, BPPARAM = SerialParam(RNGseed = WB_SEED),
                verbose = FALSE),
    error = function(e) e
  )

  if (inherits(res, "error")) {
    msg("  scDblFinder failed: %s -- sample left unflagged", conditionMessage(res))
    per_sample[[i]] <- data.frame(
      sample = sid, n_cells = n, dbr_prior = dbr,
      n_doublet = 0L, doublet_pct = 0,
      med_score_singlet = NA_real_, med_score_doublet = NA_real_,
      status = "failed", stringsAsFactors = FALSE)
    dbl_notes[idx] <- "scDblFinder_failed"
    next
  }

  sc <- res$scDblFinder.score
  cl <- as.character(res$scDblFinder.class)
  dbl_score[idx] <- sc
  dbl_class[idx] <- cl
  nd <- sum(cl == "doublet")
  msg("  doublets=%d (%.2f%%)", nd, 100 * nd / n)

  per_sample[[i]] <- data.frame(
    sample = sid, n_cells = n, dbr_prior = dbr,
    n_doublet = nd, doublet_pct = round(100 * nd / n, 2),
    med_score_singlet  = round(median(sc[cl == "singlet"]), 3),
    med_score_doublet  = if (nd > 0) round(median(sc[cl == "doublet"]), 3) else NA_real_,
    status = "ok", stringsAsFactors = FALSE)
}

per_sample <- do.call(rbind, per_sample)
meta$doublet_score <- dbl_score
meta$doublet_class <- dbl_class
meta$doublet_note  <- dbl_notes
meta$pass_doublet  <- meta$doublet_class == "singlet"

write_csv_safe(per_sample, file.path(TAB_DIR, "doublets_per_sample.csv"))

tot <- sum(!meta$pass_doublet)
msg("total doublets: %d / %d (%.2f%%)", tot, nrow(meta), 100 * tot / nrow(meta))
msg("per-sample rate range: %.2f%% - %.2f%%",
    min(per_sample$doublet_pct), max(per_sample$doublet_pct))

## ---------------------------------------------------------------------- save
## Saved before plotting on purpose: figure code is the most fragile part and a
## plotting error must not throw away 4 minutes of doublet inference.
obj$meta <- meta
obj$per_sample_doublets <- per_sample
obj$params$doublet <- list(rate_per_1k = DBL_RATE_PER_1K, max_rate = DBL_MAX_RATE,
                           method = "scDblFinder", per_sample = TRUE)
saveRDS(obj, file.path(PROC_DIR, "04_doublets.rds"), compress = "gzip")
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo_04_doublets.txt"))
msg("saved 04_doublets.rds")

## ------------------------------------------------------------------- figures
ps <- per_sample
ps$sample <- factor(ps$sample, levels = ps$sample[order(ps$doublet_pct)])

p1 <- ggplot(ps, aes(x = sample, y = doublet_pct)) +
  geom_col(fill = "#7BA7BC", width = .7) +
  geom_point(aes(y = 100 * dbr_prior), colour = "#C0504D", size = 1.6,
             shape = 124, stroke = 1.1) +
  coord_flip() +
  labs(x = NULL, y = "Doublets (%)",
       title = "Detected vs expected doublet rate",
       subtitle = "bar = detected; red tick = prior from loaded cell number") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())

dd <- data.frame(score = dbl_score, class = dbl_class,
                 nFeature = meta$nFeature_RNA, nCount = meta$nCount_RNA)
dd <- dd[is.finite(dd$score), ]
set.seed(WB_SEED)
if (nrow(dd) > 40000) dd <- dd[sample(nrow(dd), 40000), ]

p2 <- ggplot(dd, aes(nFeature, nCount, colour = class)) +
  geom_point(size = .25, alpha = .35) +
  scale_colour_manual(values = c(singlet = "#B9C6CE", doublet = "#C0504D")) +
  scale_y_log10() +
  labs(x = "Genes per cell", y = "UMIs per cell (log10)", colour = NULL,
       title = "Doublets sit at high count / high gene corner") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p3 <- ggplot(dd, aes(score, fill = class, colour = class)) +
  geom_density(alpha = .35, linewidth = .3, adjust = 1.2) +
  scale_fill_manual(values = c(singlet = "#B9C6CE", doublet = "#C0504D")) +
  scale_colour_manual(values = c(singlet = "#B9C6CE", doublet = "#C0504D")) +
  labs(x = "scDblFinder score", y = "Density", fill = NULL, colour = NULL,
       title = "Score separation") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

fig <- (p1 | p2 | p3) +
  plot_annotation(tag_levels = "A",
                  caption = sprintf("GSE176078: %d samples, %d/%d cells flagged as doublets (%.1f%%)",
                                    nrow(ps), tot, nrow(meta), 100 * tot / nrow(meta)))
ggsave(file.path(FIG_DIR, "fig08_doublets.pdf"), fig, width = 13, height = 5.2)
ggsave(file.path(FIG_DIR, "fig08_doublets.png"), fig, width = 13, height = 5.2,
       dpi = 250)

## ---------------------------------------------------------------------- done
msg("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
