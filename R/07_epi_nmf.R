## =============================================================================
## 07_epi_nmf.R -- transcriptional programs of malignant epithelium (NMF)
##
## Why programs instead of clusters: discrete clusters shift with resolution and
## seed. A program is a continuous weight per cell, so "how much of program P does
## this cell express" survives re-clustering. Clusters are used only as a display
## and naming vehicle.
##
## Learning set: stratified cells per patient. Fitting on all 15k cells would let
## the largest tumours dominate the basis vectors.
##
## Projection: each cell is projected onto the learned gene-loading matrix, so
## every cell gets a weight vector even if it was not in the learning set.
##
## Implementation note: the NMF package (0.28) needs >280 s for a single brunet
## fit at 788 x 2647 on this machine, while RcppML does the same fit in ~0.1 s
## (both wrap the same multiplicative-update family). RcppML is therefore used
## for fitting and projection; multi-run consensus bookkeeping is done here
## rather than relying on NMF::nmf internals. Scale is identical because both
## return a non-negative gene-loading matrix.
## =============================================================================

source("config/00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(RcppML)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

set.seed(WB_SEED)
t0 <- Sys.time()

so <- readRDS(file.path(PROC_DIR, "06_epi_clustered.rds"))
msg("malignant epithelial cells: %d ; clusters: %d",
    ncol(so), length(unique(so$seurat_clusters)))

expr <- GetAssayData(so, layer = "data")
feat <- intersect(VariableFeatures(so), rownames(expr))
## drop genes that are mostly zero in this compartment: they add noise to NMF
pct <- Matrix::rowMeans(expr[feat, , drop = FALSE] > 0)
feat <- feat[pct >= NMF_MIN_PCT]
msg("HVG used for NMF: %d", length(feat))

V <- as(expr[feat, , drop = FALSE], "CsparseMatrix")
V <- as(V, "dgCMatrix")
rm(expr); invisible(gc(FALSE))

## ------------------------------------------------------ stratified learning set
cells_by_sample <- split(colnames(so), so$sample)
learn <- unlist(lapply(cells_by_sample, function(cs)
  if (length(cs) <= NMF_CELLS_PER_SAMPLE) cs else sample(cs, NMF_CELLS_PER_SAMPLE)))
msg("learning set: %d cells from %d samples", length(learn), length(cells_by_sample))

vle <- V[, learn, drop = FALSE]
## log-normalised data are already non-negative; guard anyway, then scale once and
## reuse the same factor for the projection so weights stay on one scale.
vmin <- min(vle)
if (vmin < 0) vle <- vle - vmin
vmax <- max(vle)
vle <- vle / vmax
msg("value range after scaling: [0, 1] (raw max %.2f)", vmax)

## adjusted Rand index from a contingency table.
## Defined before the fitting loop on purpose: R resolves the call at run time
## and the loop body executes before any definition placed later in the file.
.ari <- function(tb) {
  n <- sum(tb); a <- rowSums(tb); b <- colSums(tb)
  sij <- sum(choose(tb, 2))
  sa <- sum(choose(a, 2)); sb <- sum(choose(b, 2))
  exp <- sa * sb / choose(n, 2)
  mx <- (sa + sb) / 2
  if (mx == exp) return(1)
  (sij - exp) / (mx - exp)
}

## --------------------------------------------- multi-run fit per candidate rank
n_learn <- ncol(vle)
n_pair  <- NMF_NRUN * (NMF_NRUN - 1) / 2
msg("fitting NMF for ranks: %s (nrun = %d, tol = %g, maxit = %d)",
    paste(NMF_RANKS, collapse = ", "), NMF_NRUN, NMF_TOL, NMF_MAXIT)

res <- lapply(NMF_RANKS, function(k) {
  tk <- Sys.time()
  runs <- lapply(seq_len(NMF_NRUN), function(i)
    nmf(vle, k = k, seed = WB_SEED + i, tol = NMF_TOL,
        maxit = NMF_MAXIT, verbose = FALSE))

  ## hard assignment of each learning cell to its dominant program
  labs <- lapply(runs, function(m) max.col(t(m$h), ties.method = "first"))

  ## consensus matrix: fraction of run pairs that co-cluster a cell pair
  cons <- matrix(0, n_learn, n_learn)
  aris <- numeric(0)
  for (a in seq_len(NMF_NRUN - 1)) for (b in (a + 1):NMF_NRUN) {
    aris <- c(aris, .ari(table(labs[[a]], labs[[b]])))
    cons <- cons + outer(labs[[a]], labs[[a]], "==") *
                   outer(labs[[b]], labs[[b]], "==")
  }
  cons <- cons / n_pair

  ## cophenetic correlation of the consensus matrix (same definition as NMF::cophcor)
  d  <- as.dist(1 - cons)
  hc <- hclust(d, method = "average")
  cp <- cor(as.vector(d), as.vector(as.dist(cophenetic(hc))))
  rm(cons, d, hc); invisible(gc(FALSE))

  ## representative run: the medoid by ARI to all other runs
  ari_mat <- matrix(NA_real_, NMF_NRUN, NMF_NRUN); idx <- 1
  for (a in seq_len(NMF_NRUN - 1)) for (b in (a + 1):NMF_NRUN) {
    ari_mat[a, b] <- ari_mat[b, a] <- aris[idx]; idx <- idx + 1
  }
  medoid <- which.max(rowMeans(ari_mat, na.rm = TRUE))

  msg("  rank %2d: cophenetic=%.3f stability_ARI=%.3f medoid_run=%d (%.1f s)",
      k, cp, mean(aris), medoid, as.numeric(difftime(Sys.time(), tk, units = "secs")))

  list(k = k, cp = cp, stab = mean(aris), W = runs[[medoid]]$w,
       d = runs[[medoid]]$d, medoid = medoid)
})

sel <- data.frame(k = sapply(res, `[[`, "k"),
                  cophenetic = sapply(res, `[[`, "cp"),
                  stability = sapply(res, `[[`, "stab"))
write_csv_safe(sel, file.path(TAB_DIR, "epi_nmf_rank_selection.csv"))

## pick rank: highest cophenetic among ranks whose stability is within 90% of max
ok <- sel[sel$stability >= 0.9 * max(sel$stability, na.rm = TRUE), ]
best_k <- ok$k[which.max(ok$cophenetic)]
best <- res[[which(sapply(res, `[[`, "k") == best_k)]]
msg("selected rank: %d", best_k)

W <- best$W                                    # genes x programs
colnames(W) <- paste0("MP", seq_len(ncol(W)))
rownames(W) <- feat

prog_genes <- lapply(seq_len(ncol(W)), function(i) {
  g <- names(sort(W[, i], decreasing = TRUE))[1:NMF_TOP_GENES]
  data.frame(program = colnames(W)[i], rank = seq_along(g), gene = g,
             weight = round(W[g, i], 5), stringsAsFactors = FALSE)
})
prog_genes <- do.call(rbind, prog_genes)
write_csv_safe(prog_genes, file.path(TAB_DIR, "epi_mp_genes.csv"))

## ------------------------------------------------------ project all cells (NNLS)
## same scaling factor as the learning set, otherwise the weights live on a
## different scale than the basis vectors they are projected onto.
A <- V / vmax
h <- project(A, w = W)                         # programs x cells
## project() returns programs x cells and carries colnames(A) over, so transpose
## first (cells x programs), row-normalise to shares, then rebuild dimnames
## instead of overwriting them (overwriting fails on length mismatch).
share <- t(h)
share <- share / pmax(rowSums(share), 1e-12)   # per-cell shares, rows sum to 1
dimnames(share) <- list(colnames(so), colnames(W))
msg("projected %d cells onto %d programs", nrow(share), ncol(share))
rm(A, h, V); invisible(gc(FALSE))

usage <- data.frame(cell_id = rownames(share), round(share, 5), check.names = FALSE)
write_csv_safe(usage, file.path(TAB_DIR, "epi_mp_usage.csv"))

## dominant program per cell
dom <- factor(colnames(share)[max.col(share, ties.method = "first")],
              levels = colnames(share))
tab <- table(dom, so$seurat_clusters)
write_csv_safe(as.data.frame.matrix(tab), file.path(TAB_DIR, "epi_mp_by_cluster.csv"))

sub <- table(dom, so$ref_subtype)
sub_pct <- round(100 * sweep(sub, 2, colSums(sub), "/"), 1)
write_csv_safe(as.data.frame.matrix(sub_pct), file.path(TAB_DIR, "epi_mp_by_subtype.csv"))

## mean program share per cluster, for the heatmap (more informative than counts)
cl_mean <- sapply(levels(so$seurat_clusters), function(cl)
  colMeans(share[so$seurat_clusters == cl, , drop = FALSE]))
cl_pct <- round(100 * sweep(cl_mean, 2, colSums(cl_mean), "/"), 2)
write_csv_safe(data.frame(program = rownames(cl_pct), cl_pct, check.names = FALSE),
               file.path(TAB_DIR, "epi_mp_cluster_share.csv"))

## mean program share per sample: shows whether a program is patient-specific
## or recurrent, which decides whether it can carry a cohort-level claim.
smp_mean <- sapply(sort(unique(so$sample)), function(s)
  colMeans(share[so$sample == s, , drop = FALSE]))
out_smp <- data.frame(sample = colnames(smp_mean), t(round(smp_mean, 4)),
                      n_cells = as.integer(table(so$sample)[colnames(smp_mean)]),
                      check.names = FALSE)
write_csv_safe(out_smp, file.path(TAB_DIR, "epi_mp_by_sample.csv"))

## -------------------------------------------------------------- figures
sel_long <- data.frame(
  k = rep(sel$k, 2),
  metric = rep(c("cophenetic", "stability (ARI)"), each = nrow(sel)),
  value = c(sel$cophenetic, sel$stability))
p1 <- ggplot(sel_long, aes(k, value, colour = metric)) +
  geom_line(linewidth = .5) + geom_point(size = 1.6) +
  geom_vline(xintercept = best_k, linetype = 2, colour = "grey40") +
  scale_x_continuous(breaks = sel$k) +
  scale_colour_manual(values = c(cophenetic = "#4A6FA5", `stability (ARI)` = "#C0504D")) +
  labs(x = "Number of programs (k)", y = "Metric", colour = NULL,
       title = "Rank selection") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

hm <- as.data.frame(tab)
names(hm) <- c("program", "cluster", "n")
hm$frac <- ave(hm$n, hm$cluster, FUN = function(z) z / sum(z))
p2 <- ggplot(hm, aes(cluster, program, fill = frac)) +
  geom_tile(colour = "white", linewidth = .3) +
  scale_fill_gradient(low = "#F2F5F7", high = "#4A6FA5",
                      labels = function(z) sprintf("%.0f%%", 100 * z)) +
  scale_y_discrete(limits = rev(levels(factor(hm$program)))) +
  labs(x = "Cluster", y = NULL, fill = "share",
       title = "Dominant program by cluster") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 6))

hsub <- as.data.frame(sub_pct)
names(hsub) <- c("program", "subtype", "pct")
hsub$subtype <- factor(hsub$subtype, levels = names(colSums(sub)))
lvl_n <- paste0(levels(hsub$subtype), "\n(n=", format(colSums(sub), big.mark = ","), ")")
p3 <- ggplot(hsub, aes(subtype, program, fill = pct)) +
  geom_tile(colour = "white", linewidth = .3) +
  geom_text(aes(label = sprintf("%.0f", pct)), size = 2.2, colour = "grey15") +
  scale_fill_gradient(low = "#F2F5F7", high = "#C0504D") +
  scale_y_discrete(limits = rev(levels(factor(hsub$program)))) +
  scale_x_discrete(labels = lvl_n) +
  labs(x = NULL, y = NULL, fill = "% of subtype\ncells",
       title = "Dominant program by subtype") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank())

fig <- p1 | p2 | p3
ggsave(file.path(FIG_DIR, "fig11_epi_nmf.png"), fig, width = 15, height = 4.6, dpi = 250)
ggsave(file.path(FIG_DIR, "fig11_epi_nmf.pdf"), fig, width = 15, height = 4.6)

## program activity on the epithelial UMAP
for (p in colnames(share)) so[[p]] <- share[, p]
## combine = FALSE keeps the panels as plain ggplots, so per-panel theme edits
## use ggplot's own `+` instead of patchwork's `&` (unreliable across versions).
pl <- FeaturePlot(so, features = colnames(share), ncol = length(colnames(share)),
                  pt.size = .1, order = TRUE, cols = c("#F2F5F7", "#4A6FA5"),
                  combine = FALSE)
pl <- lapply(seq_along(pl), function(i) {
  p <- pl[[i]] + theme(plot.title = element_text(size = 8))
  if (i < length(pl)) p <- p + NoLegend()      # one shared colour bar at the end
  p
})
fig2 <- patchwork::wrap_plots(pl, nrow = 1)
ggsave(file.path(FIG_DIR, "fig12_epi_mp_umap.png"), fig2,
       width = 2.4 * length(colnames(share)), height = 2.6, dpi = 250)
ggsave(file.path(FIG_DIR, "fig12_epi_mp_umap.pdf"), fig2,
       width = 2.4 * length(colnames(share)), height = 2.6)

saveRDS(list(W = W, usage = share, rank = best_k, sel = sel,
             genes = feat, learning_cells = learn),
        file.path(PROC_DIR, "07_epi_nmf.rds"))
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo_07_nmf.txt"))
msg("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
