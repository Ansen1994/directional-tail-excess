## =============================================================================
## 06_epi_clustered -- merge CNV shards, then subcluster malignant epithelium
##
## Design decisions
##   1. Joint clustering across subtypes, Harmony by sample only. Subtype is NOT
##      regressed out, because "which epithelial states are enriched in ER+ vs
##      TNBC" is one of the questions. Each cluster's subtype composition is
##      reported, and clusters with >90% single-subtype purity are flagged as
##      subtype-restricted rather than silently presented as shared states.
##   2. Only aneuploid (CNV-aberrant) epithelial cells enter the subclustering.
##      Normal epithelial cells are kept aside as a reference and projected
##      later, so they cannot pull the malignant manifold.
##   3. Doublets are excluded (pass_doublet).
## =============================================================================

source("config/00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

set.seed(WB_SEED)
t0 <- Sys.time()

obj <- readRDS(file.path(PROC_DIR, "04_doublets.rds"))
counts <- obj$counts
meta   <- obj$meta

cnv_file <- file.path(TAB_DIR, "cnv_cell_calls.csv")
if (!file.exists(cnv_file)) stop("run 05_cnv.R first: ", cnv_file)
cnv <- read.csv(cnv_file, stringsAsFactors = FALSE)
meta <- merge(meta, cnv[, c("cell_id", "sample", "cnv_call", "cnv_score")],
              by = c("cell_id", "sample"), all.x = TRUE)
rownames(meta) <- meta$cell_id
meta <- meta[colnames(counts), ]

EPI <- c("Cancer Epithelial", "Normal Epithelial")
is_epi   <- meta$ref_major %in% EPI
is_malig <- is_epi & !is.na(meta$cnv_call) & meta$cnv_call == "malignant" &
            meta$pass_doublet
msg("epithelial: %d ; malignant (CNV-positive): %d ; CNV-quiet epithelium: %d",
    sum(is_epi), sum(is_malig),
    sum(is_epi & meta$pass_doublet & !is.na(meta$cnv_call) & meta$cnv_call != "malignant"))

m <- meta[is_malig, ]
cnt <- counts[, rownames(m), drop = FALSE]

so <- CreateSeuratObject(counts = cnt, meta.data = m, project = "malignant_epi")
rm(cnt); invisible(gc(FALSE))

so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = CL_N_HVG, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 50, verbose = FALSE)

so <- RunHarmony(so, group.by.vars = CL_BATCH_VAR, reduction = "pca",
                 reduction.save = "harmony", verbose = FALSE)

so <- FindNeighbors(so, reduction = "harmony", dims = 1:CL_N_PCS,
                    k.param = CL_N_NEIGHBORS, verbose = FALSE)
so <- FindClusters(so, resolution = c(0.2, 0.4, 0.6, 0.8), verbose = FALSE)
so <- RunUMAP(so, reduction = "harmony", dims = 1:CL_N_PCS,
              min.dist = CL_MIN_DIST, verbose = FALSE)

so$seurat_clusters <- so$RNA_snn_res.0.4
Idents(so) <- "seurat_clusters"
msg("clusters at res 0.2/0.4/0.6/0.8: %s",
    paste(sapply(c("RNA_snn_res.0.2", "RNA_snn_res.0.4",
                   "RNA_snn_res.0.6", "RNA_snn_res.0.8"),
                 function(k) length(unique(so[[k]][, 1]))), collapse = " / "))

## ------------------------------------------------------- per-cluster summaries
emb <- as.data.frame(Embeddings(so, "umap"))
names(emb) <- c("UMAP_1", "UMAP_2")
emb$cluster <- so$seurat_clusters
emb$subtype <- so$ref_subtype
emb$sample  <- so$sample
write_csv_safe(emb, file.path(TAB_DIR, "epi_umap.csv"))

tab <- as.data.frame(table(so$seurat_clusters, so$ref_subtype))
names(tab) <- c("cluster", "subtype", "n")
tot <- tapply(tab$n, tab$cluster, sum)
tab$pct_within_cluster <- round(100 * tab$n / tot[tab$cluster], 2)

smry <- do.call(rbind, lapply(sort(unique(so$seurat_clusters)), function(cl) {
  cells <- so$seurat_clusters == cl
  st <- so$ref_subtype[cells]
  n_pat <- length(unique(so$sample[cells]))
  top <- names(sort(table(st), decreasing = TRUE))[1]
  data.frame(
    cluster = cl, n_cells = sum(cells), n_patients = n_pat,
    top_subtype = top,
    top_subtype_pct = round(100 * max(table(st)) / sum(cells), 1),
    median_genes = round(median(so$nFeature_RNA[cells]), 0),
    stringsAsFactors = FALSE)
}))
smry$subtype_restricted <- smry$top_subtype_pct >= 90
write_csv_safe(smry, file.path(TAB_DIR, "epi_cluster_summary.csv"))
write_csv_safe(tab,  file.path(TAB_DIR, "epi_cluster_by_subtype.csv"))
print(smry, row.names = FALSE)

## -------------------------------------------------------------------- figures
pal <- grDevices::hcl.colors(length(unique(so$seurat_clusters)), "Dynamic")
names(pal) <- levels(so$seurat_clusters)

p1 <- DimPlot(so, reduction = "umap", label = TRUE, repel = TRUE,
              pt.size = .25, cols = pal) +
  labs(title = sprintf("Malignant epithelium, %d clusters (res 0.4)",
                       length(pal))) +
  theme_bw(base_size = 9) + theme(panel.grid = element_blank(), legend.position = "none")

p2 <- DimPlot(so, reduction = "umap", group.by = "ref_subtype", pt.size = .25) +
  labs(title = "Subtype of origin", colour = NULL) +
  theme_bw(base_size = 9) + theme(panel.grid = element_blank(), legend.position = "top")

p3 <- DimPlot(so, reduction = "umap", group.by = "sample", pt.size = .2) +
  labs(title = "Sample of origin", colour = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), legend.position = "none")

fig <- p1 | p2 | p3
ggsave(file.path(FIG_DIR, "fig09_epi_umap.pdf"), fig, width = 16, height = 4.8)
ggsave(file.path(FIG_DIR, "fig09_epi_umap.png"), fig, width = 16, height = 4.8, dpi = 250)

saveRDS(so, file.path(PROC_DIR, "06_epi_clustered.rds"), compress = "gzip")
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo_06_epi.txt"))
msg("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
