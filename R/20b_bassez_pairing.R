# 统计 Bassez 数据的可用配对规模
# 关注：ER+ 亚型、Cancer_cell、Pre/On 配对、每例细胞数

suppressPackageStartupMessages({library(Seurat); library(data.table)})

obj <- readRDS("data/raw/bassez/bassez_downsampled.rds")
md <- as.data.table(obj@meta.data)
md[, cell := rownames(obj@meta.data)]

cat("=== layers ===\n"); print(Layers(obj))

cat("\n=== 总体 ===\n")
cat("细胞", nrow(md), " 患者", uniqueN(md$patient_id), "\n")

cat("\n=== 亚型 x 时点 ===\n")
print(dcast(md[, .N, by = .(BC_type, timepoint)], BC_type ~ timepoint, value.var = "N"))

cat("\n=== 亚型 x 细胞类型 ===\n")
print(dcast(md[, .N, by = .(BC_type, cellType)], BC_type ~ cellType, value.var = "N"))

cat("\n=== 队列 x 时点 ===\n")
print(dcast(md[, .N, by = .(cohort, timepoint)], cohort ~ timepoint, value.var = "N"))

# ---- 焦点：ER+ 的 Cancer_cell，按患者 x 时点 ----
cat("\n\n########## ER+ / Cancer_cell ##########\n")
er <- md[BC_type == "ER+" & cellType == "Cancer_cell"]
cat("细胞数", nrow(er), " 患者", uniqueN(er$patient_id), "\n")

tab <- dcast(er[, .N, by = .(patient_id, timepoint, cohort)],
             patient_id + cohort ~ timepoint, value.var = "N", fill = 0)
setnames(tab, c("patient_id", "cohort", "On", "Pre"), c("patient", "cohort", "On", "Pre"),
         skip_absent = TRUE)
if (!"On" %in% names(tab)) tab[, On := 0]
if (!"Pre" %in% names(tab)) tab[, Pre := 0]
tab[, min_n := pmin(Pre, On)]
setorder(tab, -min_n)
print(tab)

cat("\n配对（Pre>0 且 On>0）患者数：", tab[Pre > 0 & On > 0, .N], "\n")
cat("两侧均 >=30 的患者数：", tab[Pre >= 30 & On >= 30, .N], "\n")
cat("两侧均 >=50 的患者数：", tab[Pre >= 50 & On >= 50, .N], "\n")

# ---- TNBC 对照 ----
cat("\n\n########## TNBC / Cancer_cell ##########\n")
tn <- md[BC_type == "TNBC" & cellType == "Cancer_cell"]
tab2 <- dcast(tn[, .N, by = .(patient_id, timepoint)],
              patient_id ~ timepoint, value.var = "N", fill = 0)
if (!"On" %in% names(tab2)) tab2[, On := 0]
if (!"Pre" %in% names(tab2)) tab2[, Pre := 0]
tab2[, min_n := pmin(Pre, On)]
setorder(tab2, -min_n)
print(head(tab2, 25))
cat("TNBC 配对患者数：", tab2[Pre > 0 & On > 0, .N],
    " 两侧 >=30：", tab2[Pre >= 30 & On >= 30, .N], "\n")

# ---- 基因名检查：能否直接做程序评分 ----
cat("\n\n=== 关键基因是否存在 ===\n")
g <- rownames(obj)
probe <- c("ESR1", "PGR", "ERBB2", "MKI67", "EPCAM", "KRT8", "KRT18",
           "VIM", "CDH1", "CDH2", "TFF1", "AREG", "AGR2", "SCGB2A2",
           "IFIT1", "ISG15", "STAT1", "CXCL10", "COL1A1", "FN1",
           "TOP2A", "MCM2", "HLA-DRA", "CD274", "PDCD1LG2")
print(data.frame(gene = probe, present = probe %in% g,
                 row = match(probe, g)))

cat("\n=== 每个样本的细胞数上限（验证下采样）===\n")
print(summary(as.numeric(table(md$sample))))
