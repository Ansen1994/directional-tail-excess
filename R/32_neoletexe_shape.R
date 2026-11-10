#!/usr/bin/env Rscript
# ==============================================================================
# 32_neoletexe_shape.R
#
# NeoLetExe（内分泌治疗，ER+/HER2-，三时点）上复现形态判据。
#
# 要回答的问题：nc2025（化疗）与 Bassez（anti-PD1）上观察到的
# 「位移一旦发生就是整体平移，不是尾部扩张」是否在**第三种治疗类型**上仍成立。
#
# 时点：A = 基线，B = 换药（第一个 AI 用满 3 个月），C = 术前（第二个 AI 再用 3 个月）
#   A->B  第一个芳香化酶抑制剂
#   B->C  换药后
#   A->C  全程 6 个月
#
# 程序定义：与 Bassez 验证完全相同的五套文献基因集（ER 信号、增殖、IFN 应答、
#   EMT/基质、腔面分泌）。**不套用 nc2025 的 cNMF 元程序**——三者不同构，
#   混用会重犯第二轮已指出的错误。因此这里比的是"形状"，不是"幅度"。
#
# 患者编号加前缀 N（NP01...），避免与 nc2025 的 P01 混淆。
#
# 输出：
#   results/tables/neoletexe_shape.csv        每个 患者 x 程序 x 时点对 的位移与判定
#   results/tables/neoletexe_shape_summary.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(Seurat)
})

PROC <- "data/proc/neoletexe"; TAB <- "results/tables"
MIN_CELLS <- 50

PROG <- list(
  ER_signal = c("ESR1","PGR","TFF1","GREB1","PDZK1","CA12","MYB",
                "AREG","CTSD","SERPINA3"),
  Proliferation = c("MKI67","TOP2A","MCM2","PCNA","CCNB1","BUB1",
                    "AURKA","CDC20","NUSAP1","RRM2"),
  IFN_response = c("STAT1","IFIT1","IFIT3","ISG15","MX1","OAS1",
                   "CXCL10","IRF7","IFITM1","BST2"),
  EMT_matrix = c("VIM","FN1","CDH2","ZEB1","SNAI2","TWIST1",
                 "COL1A1","COL3A1","SPARC","TAGLN"),
  Luminal_secretory = c("AGR2","AGR3","SCGB2A2","SCGB1D2","PIP",
                        "KRT8","KRT18","KRT19","ELF3","FOXA1")
)

## ---------------- 1. 取恶性细胞 ----------------
cells <- fread(file.path(TAB, "neoletexe_cnv_cells.csv"))
mal <- cells[malignant == TRUE]
mal[, patient := paste0("N", patient)]     # 加前缀，与 nc2025 的同名患者区分
cat(sprintf("恶性细胞 %d 个 / %d 样本 / %d 例患者\n",
            nrow(mal), uniqueN(mal$sample), uniqueN(mal$patient)))

ep <- readRDS(file.path(PROC, "epi_combined.rds"))
cnt <- ep$counts[, mal$cell, drop = FALSE]
rm(ep); invisible(gc(FALSE))
stopifnot(identical(colnames(cnt), mal$cell))

## ---------------- 2. 打分 ----------------
obj <- CreateSeuratObject(counts = cnt, min.cells = 0, min.features = 0)
rm(cnt); invisible(gc(FALSE))
obj <- NormalizeData(obj, verbose = FALSE)
g <- rownames(obj)
PROG <- lapply(PROG, function(x) intersect(x, g))
cat("\n程序基因命中：\n")
for (n in names(PROG)) cat(sprintf("  %-18s %2d\n", n, length(PROG[[n]])))

for (n in names(PROG)) {
  obj <- AddModuleScore(obj, features = list(PROG[[n]]), name = n,
                        ctrl = 30, seed = 20261109)
}
md <- as.data.table(obj@meta.data)
md[, cell := rownames(obj@meta.data)]
rm(obj); invisible(gc(FALSE))
md <- merge(md, mal[, .(cell, patient, sample, tp)], by = "cell")
cat("打分完成：", nrow(md), " 细胞\n")

## 细胞级分数落盘，供 R/33 做克隆内分析复用
sc_cols <- paste0(names(PROG), "1")
cs <- md[, c("cell", "patient", "sample", "tp", sc_cols), with = FALSE]
setnames(cs, sc_cols, names(PROG))
fwrite(cs, file.path(TAB, "neoletexe_cell_scores.csv"))

## ---------------- 3. 分位位移与形状判定 ----------------
QS <- c(.10, .25, .50, .75, .90)
PAIRS <- list(c("A","B"), c("B","C"), c("A","C"))

## 判定规则与 R/28_fig5_shape.R 的 call_v2 完全一致：
## 展宽的门槛是"展宽幅度超过位移幅度"（|z50| < asym/2），而不是 |z50| < 0.5。
## asym 必须取 P75-P25（此时 asym = disp - 1，两者同尺度），
## 不能取 P90-P10：后者在分布单侧展宽时可以为负，abs(z50) < 负值 恒为假，
## 展宽永远判不出来（全表 129/220 个组合受影响）。另外 cohort 1 用的也是
## P75-P25，两队列必须同口径才能比构成。
classify <- function(z10, z25, z50, z75, z90, disp, asym) {
  if (!is.finite(z50) || !is.finite(disp)) return("insufficient")
  if (disp > 1.6 && abs(z50) < asym / 2)  return("dispersion up")
  if (disp < 0.63 && abs(z50) < 0.5)      return("dispersion down")
  if (z25 > 0.3 && z50 > 0.3 && z75 > 0.3) return("whole shift up")
  if (z25 < -0.3 && z50 < -0.3 && z75 < -0.3) return("whole shift down")
  if (z75 > 0.5 && abs(z25) < 0.35 && abs(z50) < 0.35) return("right-tail expansion")
  if (z25 > 0.5 && abs(z75) < 0.35) return("left truncation")
  "no change"
}

rows <- list()
for (pt in unique(md$patient)) {
  for (pr in PAIRS) {
    d0 <- md[patient == pt & tp == pr[1]]
    d1 <- md[patient == pt & tp == pr[2]]
    if (nrow(d0) < MIN_CELLS || nrow(d1) < MIN_CELLS) next
    for (n in names(PROG)) {
      col <- paste0(n, "1")
      a <- d0[[col]]; b <- d1[[col]]
      if (!is.finite(var(a)) || var(a) == 0) next
      qa <- quantile(a, QS, type = 7, names = FALSE)
      qb <- quantile(b, QS, type = 7, names = FALSE)
      iqr <- qa[4] - qa[2]          # 基线 IQR（Q75-Q25），QS 里 Q25/Q75 是第 2/4 位
      if (!is.finite(iqr) || iqr <= 0) next
      z <- (qb - qa) / iqr
      disp <- (qb[4] - qb[2]) / iqr
      asym <- z[4] - z[2]          # P75-P25：与 disp-1 同尺度，且与 cohort 1 一致
      rows[[length(rows) + 1]] <- data.table(
        cohort = "NeoLetExe", patient = pt, pair = paste(pr, collapse = ">"),
        from = pr[1], to = pr[2], program = n,
        n_from = nrow(d0), n_to = nrow(d1),
        q10_from = qa[1], q50_from = qa[3], q90_from = qa[5],
        q10_to = qb[1], q50_to = qb[3], q90_to = qb[5],
        z10 = z[1], z25 = z[2], z50 = z[3], z75 = z[4], z90 = z[5],
        disp = disp, asym = asym,
        call_v2 = classify(z[1], z[2], z[3], z[4], z[5], disp, asym))
    }
  }
}
res <- rbindlist(rows)
cat("\n组合数：", nrow(res), " 患者：", uniqueN(res$patient), "\n")

## 定向尾部超额：lead = 位移方向远端，trail = 另一端（与 R/28 同定义）
res[, `:=`(ex_lead  = fifelse(z50 > 0, z75 - z50, z25 - z50) * sign(z50),
           ex_trail = fifelse(z50 > 0, z25 - z50, z75 - z50) * sign(z50))]
fwrite(res, file.path(TAB, "neoletexe_shape.csv"))

## ---------------- 4. 汇总 ----------------
cat("\n=== 判定分布（按时点对）===\n")
print(dcast(res[, .N, by = .(pair, call_v2)], pair ~ call_v2, value.var = "N", fill = 0))
cat("\n=== 判定分布（全体）===\n")
print(res[, .N, by = call_v2][order(-N)])

cat("\n=== 位移幅度的分布 ===\n")
print(res[, .(n = .N, med_z50 = round(median(z50), 3),
              q90_abs = round(quantile(abs(z50), .9), 3),
              frac_abs_gt_.5 = round(mean(abs(z50) > .5), 3)), by = pair])

cat("\n=== 关键读数 ===\n")
big <- res[abs(z50) > 0.5]
cat(sprintf("实质位移（|z50|>0.5）：%d / %d (%.1f%%)\n",
            nrow(big), nrow(res), 100 * nrow(big) / nrow(res)))
cat(sprintf("其中尾部变化：%d 例（右尾 %d，左截 %d）\n",
            nrow(big[call_v2 %in% c("right-tail expansion", "left truncation")]),
            nrow(big[call_v2 == "right-tail expansion"]),
            nrow(big[call_v2 == "left truncation"])))
cat(sprintf("上移组合 trail 超额中位 %.3f，最大 %.3f\n",
            median(big[z50 > 0]$ex_trail), max(big[z50 > 0]$ex_trail)))
cat(sprintf("lead 超额中位 %.3f，最大 %.3f\n",
            median(big$ex_lead), max(big$ex_lead)))

summ <- res[, .(n = .N,
                n_shift = sum(call_v2 %in% c("whole shift up", "whole shift down")),
                n_tail = sum(call_v2 %in% c("right-tail expansion", "left truncation")),
                med_z50 = round(median(z50), 3),
                med_disp = round(median(disp), 2)), by = .(pair, program)]
fwrite(summ, file.path(TAB, "neoletexe_shape_summary.csv"))
cat("\n=== 按时点对 x 程序 ===\n"); print(summ)
cat("\n-> neoletexe_shape.csv / neoletexe_shape_summary.csv\n")
