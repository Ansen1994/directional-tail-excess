## =============================================================================
## 02_integrate.R - normalisation, feature selection, batch integration,
##                   multi-resolution clustering and UMAP
##
## Input : data/processed/01_counts_qc.rds
## Output: data/processed/02_clustered.rds
##         results/tables/cluster_sizes.csv
##         results/tables/cluster_vs_sample.csv
##         results/figures/fig02_umap_clusters.pdf/png
##         results/figures/fig03_umap_qc.pdf/png
##         results/figures/fig04_clustree.pdf/png
##
## Design notes
##   Integration uses Harmony on the sample variable, no external reference,
##   no cell removal, so the intrinsic structure is preserved. Whether the batch
##   effect is gone is checked visually in fig02 (colour by sample), not assumed.
##   Nine resolutions are computed in one pass and ALL are kept, so the main
##   resolution is not chosen post hoc by which one looks prettiest.
##   NOTE: keep this file pure ASCII. R under LC_CTYPE=C on this machine fails to
##   parse some multibyte string literals in source files.
## =============================================================================

source("config/00_config.R")
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(ggplot2)
  library(patchwork); library(clustree)
})

set.seed(WB_SEED)
t0 <- Sys.time()

in_rds <- file.path(PROC_DIR, "01_counts_qc.rds")
if (!file.exists(in_rds)) stop("missing ", in_rds, "; run R/01_qc.R first")
inp <- readRDS(in_rds)
msg("loaded %d genes x %d cells", nrow(inp$counts), ncol(inp$counts))

obj <- CreateSeuratObject(counts = inp$counts, meta.data = inp$meta,
                          min.cells = 0, min.features = 0, project = "GSE176078")
rm(inp); invisible(gc(FALSE))

## --------------------------------------------------- normalise + HVG + PCA
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = CL_N_HVG, verbose = FALSE)
hvgs <- VariableFeatures(obj)
msg("highly variable genes: %d", length(hvgs))

obj <- ScaleData(obj, features = hvgs, verbose = FALSE)
obj <- RunPCA(obj, features = hvgs, npcs = CL_N_PCS, verbose = FALSE)
msg("PCA done: %d components", CL_N_PCS)

## ------------------------------------------------------------- integration
if (!CL_SKIP_HARMONY && length(unique(obj[[CL_BATCH_VAR]][, 1])) > 1) {
  obj <- RunHarmony(obj, group.by.vars = CL_BATCH_VAR, reduction.use = "pca",
                    reduction.save = "harmony", verbose = FALSE)
  red <- "harmony"
  msg("harmony done, batch variable = %s (%d levels)", CL_BATCH_VAR,
      length(unique(obj[[CL_BATCH_VAR]][, 1])))
} else {
  red <- "pca"
  msg("harmony skipped, using PCA space")
}

## ------------------------------------------- neighbour graph + clustering
obj <- FindNeighbors(obj, reduction = red, dims = seq_len(CL_N_PCS), verbose = FALSE)
obj <- FindClusters(obj, resolution = CL_RESOLUTIONS, verbose = FALSE)

res_cols <- grep("^RNA_snn_res\\.", colnames(obj@meta.data), value = TRUE)
main_col <- paste0("RNA_snn_res.", CL_MAIN_RES)
if (!main_col %in% res_cols) {
  main_col <- res_cols[which.min(abs(as.numeric(sub("^RNA_snn_res\\.", "", res_cols)) - CL_MAIN_RES))]
  msg("requested resolution %.2f not found, using %s", CL_MAIN_RES, main_col)
}
obj$seurat_clusters <- obj[[main_col]][, 1]
Idents(obj) <- "seurat_clusters"
msg("main resolution %s: %d clusters", main_col, length(unique(obj$seurat_clusters)))

## ------------------------------------------------------------------- UMAP
obj <- RunUMAP(obj, reduction = red, dims = seq_len(CL_N_PCS),
               n.neighbors = CL_N_NEIGHBORS, min.dist = CL_MIN_DIST,
               seed.use = WB_SEED, verbose = FALSE)
obj@misc$reduction_used <- red
obj@misc$main_res_col   <- main_col
obj@misc$params <- list(seed = WB_SEED, n_hvg = CL_N_HVG, n_pcs = CL_N_PCS,
                        n_neighbors = CL_N_NEIGHBORS, min_dist = CL_MIN_DIST,
                        batch_var = CL_BATCH_VAR, resolutions = CL_RESOLUTIONS,
                        main_res = CL_MAIN_RES)

## -------------------------------------------------------------- statistics
sizes <- do.call(rbind, lapply(res_cols, function(cl) {
  tab <- table(obj[[cl]][, 1])
  data.frame(resolution = as.numeric(sub("^RNA_snn_res\\.", "", cl)),
             n_clusters = length(tab),
             min_size = min(tab), median_size = median(as.numeric(tab)),
             max_size = max(tab),
             singleton_clusters = sum(tab <= 3),
             stringsAsFactors = FALSE)
}))
write_csv_safe(sizes, file.path(TAB_DIR, "cluster_sizes.csv"))
msg("clusters per resolution: %s", paste(sizes$n_clusters, collapse = ", "))

cross <- as.data.frame.matrix(table(obj$seurat_clusters, obj$sample))
cross <- cbind(cluster = rownames(cross), cross)
write_csv_safe(cross, file.path(TAB_DIR, "cluster_vs_sample.csv"))

## ------------------------------------------------------------------ figures
umap_df <- as.data.frame(Embeddings(obj, "umap"))
names(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cluster <- obj$seurat_clusters
umap_df$sample  <- obj$sample
cent <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = umap_df, FUN = median)

p_cl <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = cluster)) +
  geom_point(size = .12, alpha = .5) +
  geom_text(data = cent, aes(label = cluster), colour = "black", size = 3.2,
            fontface = "bold") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1),
                               ncol = 2)) +
  labs(title = sprintf("%d clusters (resolution %.2f)",
                       length(unique(umap_df$cluster)),
                       as.numeric(sub("^RNA_snn_res\\.", "", main_col))),
       x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), legend.key.height = unit(9, "pt"),
        plot.title = element_text(size = 9))

p_sm <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = sample)) +
  geom_point(size = .12, alpha = .45) +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 2)) +
  labs(title = "By sample (batch check)", x = "UMAP 1", y = "UMAP 2", colour = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), legend.key.height = unit(8, "pt"),
        legend.text = element_text(size = 6), plot.title = element_text(size = 9))

fig2 <- p_cl | p_sm
ggsave(file.path(FIG_DIR, "fig02_umap_clusters.pdf"), fig2, width = 12, height = 5.6)
ggsave(file.path(FIG_DIR, "fig02_umap_clusters.png"), fig2, width = 12, height = 5.6, dpi = 300)

qc_plot <- function(var, lab, use_log = TRUE) {
  d <- umap_df; d$v <- obj[[var]][, 1]
  ggplot(d, aes(UMAP_1, UMAP_2, colour = v)) +
    geom_point(size = .12, alpha = .5) +
    scale_colour_gradientn(colours = c("#3B5B7A", "#7BA7BC", "#9CAF88", "#E8A87C"),
                           transform = if (use_log) "log10" else "identity") +
    labs(title = lab, x = "UMAP 1", y = "UMAP 2", colour = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 9),
          legend.key.height = unit(8, "pt"))
}
## percent_mito contains zeros; log10 would give -Inf, so keep it linear
fig3 <- qc_plot("nFeature_RNA", "Genes per cell") |
        qc_plot("nCount_RNA", "UMIs per cell") |
        qc_plot("percent_mito", "Mitochondrial %", use_log = FALSE)
ggsave(file.path(FIG_DIR, "fig03_umap_qc.pdf"), fig3, width = 14, height = 4.8)
ggsave(file.path(FIG_DIR, "fig03_umap_qc.png"), fig3, width = 14, height = 4.8, dpi = 300)

p_tree <- clustree(obj, prefix = "RNA_snn_res.", layout = "tree") +
  guides(edge_colour = "none", edge_alpha = "none", size = "none") +
  theme(legend.position = "right")
ggsave(file.path(FIG_DIR, "fig04_clustree.pdf"), p_tree, width = 11, height = 12)
ggsave(file.path(FIG_DIR, "fig04_clustree.png"), p_tree, width = 11, height = 12, dpi = 200)

## -------------------------------------------------------------------- save
saveRDS(obj, file.path(PROC_DIR, "02_clustered.rds"), compress = "gzip")

con <- file(file.path(LOG_DIR, "02_integrate_summary.txt"), "w")
writeLines(c(
  sprintf("run at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("cells: %d ; genes: %d", ncol(obj), nrow(obj)),
  sprintf("reduction: %s (batch var %s) ; HVG=%d ; PCs=%d",
          red, CL_BATCH_VAR, CL_N_HVG, CL_N_PCS),
  sprintf("main resolution: %s -> %d clusters", main_col,
          length(unique(obj$seurat_clusters))),
  sprintf("clusters per resolution: %s",
          paste(sprintf("%.1f:%d", sizes$resolution, sizes$n_clusters), collapse = "  ")),
  sprintf("smallest cluster: %d cells ; clusters with <=3 cells: %d",
          min(sizes$min_size), max(sizes$singleton_clusters)),
  sprintf("elapsed: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
), con); close(con)
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo_02_integrate.txt"))
msg("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
