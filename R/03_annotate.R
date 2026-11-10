## =============================================================================
## 03_annotate.R - independent annotation and agreement with the original
##                 curated labels
##
## Input : data/processed/02_clustered.rds
##         results/tables/labels_reference.csv  (produced by 01b_attach_labels.R)
## Output: results/tables/cluster_markers.csv
##         results/tables/cluster_annotation.csv
##         results/tables/annotation_agreement.csv
##         results/tables/confusion_cluster_vs_ref.csv
##         results/figures/fig05_markers_dotplot.pdf/png
##         results/figures/fig06_annotation_umap.pdf/png
##         results/figures/fig07_agreement.pdf/png
##
## Method
##   Three independent lines of evidence are generated and then reconciled,
##   so that no single method is used to validate itself:
##   (1) data-driven cluster markers  (FindAllMarkers, Wilcoxon, presto)
##   (2) lineage module scores        (mean expression of 14 canonical marker
##                                    sets per cluster, z-scored across clusters)
##   (3) reference mapping            (SingleR + celldex HPCA)
##   Agreement is quantified with ARI and NMI: can an independent pipeline
##   reconstruct the published partitions? That is the reproducibility claim.
##
## NOTE: keep this file pure ASCII (see 02_integrate.R header).
## =============================================================================

source("config/00_config.R")
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
  library(SingleR); library(celldex); library(igraph); library(dplyr)
})

set.seed(WB_SEED)
t0 <- Sys.time()

obj <- readRDS(file.path(PROC_DIR, "02_clustered.rds"))
lab <- read.csv(file.path(TAB_DIR, "labels_reference.csv"), stringsAsFactors = FALSE)
msg("loaded %d cells, %d clusters", ncol(obj), length(unique(obj$seurat_clusters)))

## Free ~1.6 GB: scale.data is only needed by PCA, which is already done.
## Without this, SingleR runs out of memory on this machine (measured).
obj[["RNA"]]$scale.data <- NULL
obj@graphs <- list()
invisible(gc(FALSE))
msg("released scale.data and graphs before annotation")

## ------------------------------------------- attach original curated labels
idx <- match(colnames(obj), lab$cell_id)
if (all(is.na(idx))) stop("cell_id does not match labels_reference.csv")
for (cl in c("ref_subtype", "ref_major", "ref_minor", "ref_subset"))
  obj[[cl]] <- lab[[cl]][idx]
msg("label match rate: %.2f%%", 100 * mean(!is.na(idx)))

## -------------------------------------------------------- (1) cluster markers
Idents(obj) <- "seurat_clusters"
marker_csv <- file.path(TAB_DIR, "cluster_markers.csv")
if (file.exists(marker_csv) && Sys.getenv("WB_FORCE_MARKERS", "0") != "1") {
  markers <- read.csv(marker_csv, stringsAsFactors = FALSE)
  msg("markers loaded from cache (%d records)", nrow(markers))
} else {
  msg("computing cluster markers (Wilcoxon with presto) ...")
  markers <- FindAllMarkers(obj, only.pos = DE_ONLY_POS, min.pct = DE_MIN_PCT,
                            logfc.threshold = DE_LOGFC, verbose = FALSE)
  markers$cluster <- as.character(markers$cluster)
  write_csv_safe(markers, marker_csv)
  msg("%d marker records written", nrow(markers))
}
markers$cluster <- as.character(markers$cluster)
top_n <- 20
top <- markers %>% group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = top_n) %>% ungroup()

## ------------------------------------- (2) lineage scores per cluster
## Deliberately NOT AddModuleScore: its random control bins are noisy on sparse
## 10x data and produced lineages that contradicted the cluster's own markers
## (e.g. cluster 9, whose markers are PIGR/KRT15/GABRP, scored highest as Mast).
## Instead: mean log-normalised expression of each marker set per cluster, then
## z-scored across clusters. Transparent, and the audit table shows the values.
lin <- LINEAGE_MARKERS[vapply(LINEAGE_MARKERS,
                              function(g) sum(g %in% rownames(obj)) >= 3, logical(1))]
cl_levels <- levels(obj$seurat_clusters)
lin_avg <- AverageExpression(obj, features = unlist(lin), group.by = "seurat_clusters",
                             assays = "RNA", layer = "data", verbose = FALSE)$RNA
lin_mean <- sapply(lin, function(g) {
  gg <- intersect(g, rownames(lin_avg))
  colMeans(as.matrix(lin_avg[gg, , drop = FALSE]))
})
lin_mean <- matrix(lin_mean, nrow = length(cl_levels),
                   dimnames = list(cl_levels, names(lin)))
lin_z <- scale(lin_mean)                      # z across clusters: relative pattern only
lin_best <- colnames(lin_z)[max.col(lin_z, ties.method = "first")]
names(lin_best) <- rownames(lin_z)
lin_top_z <- apply(lin_z, 1, max)
rm(lin_avg); invisible(gc(FALSE))
msg("lineage scoring done over %d lineages (mean expression + z across clusters)",
    length(lin))

## ------------------------------------------------------------- (3) SingleR
## Memory is the binding constraint on this machine (31.8 GB total). The
## expression slice is taken first, then the reference is loaded, then SingleR
## runs; each step is followed by an explicit gc(). Feature set is restricted to
## HVGs + lineage markers, which is both faster and sufficient for label transfer.
msg("preparing SingleR input ...")
feat_use <- unique(c(VariableFeatures(obj), unlist(LINEAGE_MARKERS)))
feat_use <- intersect(feat_use, rownames(obj))
msg("SingleR feature set: %d genes", length(feat_use))
sr_cache <- file.path(TAB_DIR, "singler_calls.csv")
if (file.exists(sr_cache) && Sys.getenv("WB_FORCE_SINGLER", "0") != "1") {
  sr <- read.csv(sr_cache, stringsAsFactors = FALSE)
  stopifnot(identical(sr$cell_id, colnames(obj)))
  obj$singler_main <- sr$singler_main
  obj$singler_delta <- sr$singler_delta
  msg("SingleR calls loaded from cache")
} else {
  test_mat <- GetAssayData(obj, layer = "data")[feat_use, , drop = FALSE]
  invisible(gc(FALSE))
  msg("running SingleR against HumanPrimaryCellAtlasData ...")
  hpca <- celldex::HumanPrimaryCellAtlasData()
  pred <- SingleR(test = test_mat, ref = hpca, labels = hpca$label.main)
  rm(test_mat, hpca); invisible(gc(FALSE))

  ## delta = top1 - top2 score; vectorised instead of apply() over 95k rows
  sc2 <- pred$scores
  r_i <- seq_len(nrow(sc2))
  t1  <- max.col(sc2, ties.method = "first")
  s1  <- sc2[cbind(r_i, t1)]
  sc2[cbind(r_i, t1)] <- -Inf
  s2  <- sc2[cbind(r_i, max.col(sc2, ties.method = "first"))]
  obj$singler_main  <- pred$labels
  obj$singler_delta <- s1 - s2
  rm(sc2, pred); invisible(gc(FALSE))
  write_csv_safe(data.frame(cell_id = colnames(obj),
                            singler_main = obj$singler_main,
                            singler_delta = obj$singler_delta,
                            stringsAsFactors = FALSE), sr_cache)
}
n_sl <- length(unique(obj$singler_main))
msg("SingleR ready: %d reference labels", n_sl)

## Map the 36 HPCA labels onto the same 14 lineage names used above, so that
## agreement is judged on a common vocabulary rather than on string equality.
SINGLER_TO_LINEAGE <- c(
  Epithelial_cells = "Epithelial", Keratinocytes = "Epithelial",
  T_cells = "T_cell", NK_cell = "T_cell",
  "Pre.B_cell_CD34." = "B_cell", "Pro.B_cell_CD34." = "B_cell",
  B_cell = "B_cell",
  Plasma_cells = "Plasma", Macrophage = "Myeloid", Monocyte = "Myeloid",
  Neutrophils = "Myeloid", Myeloid_DCs = "Myeloid", DC = "Myeloid",
  Platelets = "Myeloid", Megakaryocytes = "Myeloid",
  Endothelial_cells = "Endothelial", Fibroblasts = "Fibroblast",
  Chondrocytes = "Fibroblast", Mesangial_cells = "Fibroblast",
  Smooth_muscle_cells = "Perivascular", Tissue_stem_cells = "Perivascular",
  MSC = "Perivascular", Adipocytes = "Adipocyte", Mast_cells = "Mast"
)
map_lin <- function(x) {
  k <- gsub("[^A-Za-z0-9_]", ".", x)
  out <- unname(SINGLER_TO_LINEAGE[k])
  ifelse(is.na(out), "Other", out)
}
obj$singler_lineage <- map_lin(obj$singler_main)

## ------------------------------------------------------ agreement metrics
ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  c2 <- function(v) sum(v * (v - 1) / 2)
  sum_ij <- c2(as.vector(tab)); a <- c2(rowSums(tab)); b <- c2(colSums(tab))
  exp_idx <- a * b / (n * (n - 1) / 2)
  max_idx <- (a + b) / 2
  (sum_ij - exp_idx) / (max_idx - exp_idx)
}
comm1 <- as.integer(factor(obj$seurat_clusters))
comm2 <- as.integer(factor(obj$ref_major))
comm3 <- as.integer(factor(obj$ref_minor))
comm4 <- as.integer(factor(obj$singler_main))
metrics <- data.frame(
  comparison = c("independent clusters vs original major",
                 "independent clusters vs original minor",
                 "SingleR vs original major"),
  ari = c(ari(comm1, comm2), ari(comm1, comm3), ari(comm4, comm2)),
  nmi = c(compare(comm1, comm2, method = "nmi"),
          compare(comm1, comm3, method = "nmi"),
          compare(comm4, comm2, method = "nmi")),
  stringsAsFactors = FALSE)
metrics[, 2:3] <- round(metrics[, 2:3], 4)
write_csv_safe(metrics, file.path(TAB_DIR, "annotation_agreement.csv"))
for (i in seq_len(nrow(metrics)))
  msg("agreement: %s -> ARI=%.3f NMI=%.3f", metrics$comparison[i],
      metrics$ari[i], metrics$nmi[i])

## ------------------------------------------------------------ per-cluster table
conf <- table(obj$seurat_clusters, obj$ref_major)
write_csv_safe(cbind(cluster = rownames(conf), as.data.frame.matrix(conf)),
               file.path(TAB_DIR, "confusion_cluster_vs_ref.csv"))

majority <- function(v) {
  t <- table(v, useNA = "no")
  if (!length(t)) return(c(lab = NA_character_, p = 0))
  c(lab = names(t)[which.max(t)], p = max(t) / sum(t))
}

anno <- do.call(rbind, lapply(cl_levels, function(cl) {
  cells <- obj$seurat_clusters == cl
  mj <- majority(obj$ref_major[cells])
  sj <- majority(obj$singler_main[cells])
  sl <- majority(obj$singler_lineage[cells])
  data.frame(
    cluster = cl, n_cells = sum(cells),
    top_markers = paste(head(top$gene[top$cluster == cl], 8), collapse = ", "),
    lineage_score = lin_best[cl],
    lineage_z = round(lin_top_z[cl], 2),
    singler = sj[["lab"]], singler_frac = round(as.numeric(sj[["p"]]), 3),
    singler_lineage = sl[["lab"]],
    ref_major = mj[["lab"]], ref_purity = round(as.numeric(mj[["p"]]), 3),
    stringsAsFactors = FALSE)
}))

## Confidence rules: the more the three sources agree, the more trustworthy.
## Agreement between the lineage score and SingleR is judged on the 14-name
## vocabulary (singler_lineage), not on raw string equality, because HPCA and
## the marker sets use different naming systems.
anno$concordant <- ifelse(
  !is.na(anno$ref_major) & anno$ref_purity >= 0.8, "ref-dominant",
  ifelse(anno$lineage_score == anno$singler_lineage, "lineage=SingleR", "conflict"))
anno$suggested <- with(anno, ifelse(concordant != "conflict",
                                    ifelse(!is.na(ref_major), ref_major, singler),
                                    paste0(lineage_score, " (unresolved)")))
anno$confidence <- with(anno, ifelse(
  concordant == "ref-dominant" & ref_purity >= 0.9, "high",
  ifelse(concordant == "conflict", "low", "moderate")))
write_csv_safe(anno, file.path(TAB_DIR, "cluster_annotation.csv"))
msg("annotation calls: high %d / moderate %d / low %d",
    sum(anno$confidence == "high"), sum(anno$confidence == "moderate"),
    sum(anno$confidence == "low"))

## ---------------------------------------------------------------- figures
## fig05 uses the canonical lineage marker panel rather than the top-6 markers
## of every cluster: with 34 clusters the latter gives 204 gene labels, which is
## unreadable at any sensible figure width. Full per-cluster markers remain in
## results/tables/cluster_markers.csv.
key_genes <- intersect(unique(unlist(LINEAGE_MARKERS)), rownames(obj))
p5 <- DotPlot(obj, features = key_genes, group.by = "seurat_clusters") +
  RotatedAxis() +
  labs(x = NULL, y = "Cluster",
       title = sprintf("Canonical lineage markers (%d genes)", length(key_genes))) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = .5),
        axis.text.y = element_text(size = 7),
        panel.grid.major = element_line(linewidth = .15, colour = "grey92"),
        plot.title = element_text(size = 9))
ggsave(file.path(FIG_DIR, "fig05_markers_dotplot.pdf"), p5, width = 15, height = 8.5)
ggsave(file.path(FIG_DIR, "fig05_markers_dotplot.png"), p5, width = 15, height = 8.5, dpi = 250)

emb <- as.data.frame(Embeddings(obj, "umap")); names(emb) <- c("UMAP_1", "UMAP_2")
emb$cluster <- obj$seurat_clusters
emb$ref_major <- obj$ref_major
emb$singler <- obj$singler_main
emb$cluster_num <- as.numeric(as.character(obj$seurat_clusters))
cent <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = emb, FUN = median)
p_umap <- function(col, ttl, ncol_leg = 1, leg_size = 5.5, manual = NULL) {
  p <- ggplot(emb, aes(UMAP_1, UMAP_2, colour = .data[[col]])) +
    geom_point(size = .1, alpha = .45) +
    guides(colour = guide_legend(override.aes = list(size = 2.2), ncol = ncol_leg)) +
    labs(title = ttl, x = "UMAP 1", y = "UMAP 2", colour = NULL) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), legend.key.height = unit(6.5, "pt"),
          legend.text = element_text(size = leg_size),
          plot.title = element_text(size = 8.5))
  if (!is.null(manual)) p <- p + scale_colour_manual(values = manual)
  p
}
## 34 discrete colours are not separable in one panel; use an ordered viridis
## ramp plus on-plot cluster numbers, and drop the legend (map is in fig07).
p6a <- ggplot(emb, aes(UMAP_1, UMAP_2, colour = cluster_num)) +
  geom_point(size = .1, alpha = .45) +
  geom_text(data = cent, aes(label = cluster), colour = "black", size = 2.6,
            fontface = "bold", show.legend = FALSE) +
  scale_colour_viridis_c(option = "D", direction = -1, guide = "none") +
  labs(title = sprintf("Independent clusters (n=%d)", length(cl_levels)),
       x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), plot.title = element_text(size = 8.5))

p6 <- p6a |
      p_umap("ref_major", "Original curated labels", ncol_leg = 1, leg_size = 6) |
      p_umap("singler", "SingleR (HPCA reference)", ncol_leg = 2, leg_size = 4.5)
ggsave(file.path(FIG_DIR, "fig06_annotation_umap.pdf"), p6, width = 18, height = 6)
ggsave(file.path(FIG_DIR, "fig06_annotation_umap.png"), p6, width = 18, height = 6, dpi = 250)

cdf <- as.data.frame(conf); names(cdf) <- c("cluster", "ref_major", "n")
cdf$frac <- cdf$n / ave(cdf$n, cdf$cluster, FUN = sum)
p7a <- ggplot(cdf, aes(cluster, ref_major, fill = frac)) +
  geom_tile(colour = "white", linewidth = .3) +
  scale_fill_gradientn(colours = c("#F5F5F5", "#A8C6D8", "#4A6FA5"),
                       labels = function(x) sprintf("%.0f%%", 100 * x)) +
  labs(x = "Independent cluster", y = "Original major label",
       fill = "Within-\ncluster\nshare") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 5.5,
                                   margin = margin(t = 2)),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        plot.margin = margin(t = 6, r = 6, b = 12, l = 6))

p7b <- ggplot(metrics, aes(x = reorder(comparison, ari), y = ari)) +
  geom_col(fill = "#7BA7BC", width = .6) +
  geom_text(aes(label = sprintf("%.3f", ari)), hjust = -.15, size = 2.6) +
  coord_flip() + ylim(0, max(1, max(metrics$ari) * 1.15)) +
  labs(x = NULL, y = "Adjusted Rand Index", title = "Clustering reproducibility") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 8.5),
        plot.margin = margin(t = 6, r = 12, b = 12, l = 6))

p7 <- p7a | p7b
ggsave(file.path(FIG_DIR, "fig07_agreement.pdf"), p7, width = 13, height = 7.5)
ggsave(file.path(FIG_DIR, "fig07_agreement.png"), p7, width = 13, height = 7.5, dpi = 250)

## ------------------------------------------------------------------- save
## Lightweight annotation record only: the full expression object stays in
## 02_clustered.rds. Saving a second full Seurat object would double disk use
## (~1 GB each) for no analytical gain, and downstream scripts can re-read 02.
out <- list(
  meta = obj@meta.data,
  umap = Embeddings(obj, "umap"),
  anno = anno,
  markers_top = top,
  agreement = metrics,
  cluster_levels = cl_levels,
  lineage_z = lin_z,
  params = obj@misc$params,
  note = "Lightweight annotation record. Full counts/data layers are in 02_clustered.rds."
)
saveRDS(out, file.path(PROC_DIR, "03_annotated.rds"), compress = "gzip")
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo_03_annotate.txt"))
msg("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
