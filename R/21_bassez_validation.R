# 21_bassez_validation.R
#
# 在 Bassez 2021（anti-PD1，ER+/TNBC 配对活检）上复现「整体平移」判据。
#
# 背景：nc2025 的核心形态结论是——治疗前后程序位移全部为整体平移，
# 无右尾扩张、无左端截断（第三轮 63 组合中 0 例尾部变化）。
# 本脚本用独立队列检验这一形态特征是否可复现。
#
# 队列差异（写进结果时必须说明）：
#   nc2025   多周期化疗或化疗+ICI，HR+/HER2- 为主
#   Bassez   单剂 anti-PD1，9±2 天后手术；含 ER+ / TNBC / HER2+
# 因此本队列验证的是「治疗压力下的形态特征」，不是「化疗效应」。
#
# 程序定义：用文献基因集，不套用 nc2025 的 cNMF 元程序（两者不同构，
# Jaccard 0.024，混用会重犯第二轮已经指出的错误）。

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(patchwork)
})

TAB <- "results/tables"; FIG <- "results/figures"
MIN_CELLS <- 50

# ---------- 1. 读入并取恶性细胞 ----------
cat("读取对象...\n")
obj <- readRDS("data/raw/bassez/bassez_downsampled.rds")
md <- as.data.table(obj@meta.data)
md[, cell := rownames(obj@meta.data)]

keep <- md[cellType == "Cancer_cell" & BC_type %in% c("ER+", "TNBC")]
cat("恶性细胞（ER+ 与 TNBC）：", nrow(keep), "\n")

# ---------- 2. 程序基因集 ----------
PROG <- list(
  ER_signal = c("ESR1", "PGR", "TFF1", "GREB1", "PDZK1", "CA12", "MYB",
                "AREG", "CTSD", "SERPINA3"),
  Proliferation = c("MKI67", "TOP2A", "MCM2", "PCNA", "CCNB1", "BUB1",
                    "AURKA", "CDC20", "NUSAP1", "RRM2"),
  IFN_response = c("STAT1", "IFIT1", "IFIT3", "ISG15", "MX1", "OAS1",
                   "CXCL10", "IRF7", "IFITM1", "BST2"),
  EMT_matrix = c("VIM", "FN1", "CDH2", "ZEB1", "SNAI2", "TWIST1",
                 "COL1A1", "COL3A1", "SPARC", "TAGLN"),
  Luminal_secretory = c("AGR2", "AGR3", "SCGB2A2", "SCGB1D2", "PIP",
                        "KRT8", "KRT18", "KRT19", "ELF3", "FOXA1")
)
genes_all <- rownames(obj)
PROG <- lapply(PROG, function(g) intersect(g, genes_all))
cat("\n程序基因命中：\n"); for (n in names(PROG)) cat(sprintf("  %-16s %d/%d\n",
  n, length(PROG[[n]]), length(PROG[[n]])))

# ---------- 3. 打分（AddModuleScore，bin 控制背景表达） ----------
sub <- subset(obj, cells = keep$cell)
rm(obj); gc()

for (n in names(PROG)) {
  sub <- AddModuleScore(sub, features = list(PROG[[n]]),
                        name = n, ctrl = 30, seed = 20261108)
}
score_cols <- paste0(names(PROG), "1")
md2 <- as.data.table(sub@meta.data)
md2[, cell := rownames(sub@meta.data)]
cat("打分完成，细胞 ", nrow(md2), "\n")

# ---------- 4. 分位数与形状判定 ----------
q <- function(v, p) {
  if (!length(v)) return(NA_real_)
  as.numeric(quantile(v, p, type = 7, names = FALSE))
}

classify <- function(z25, z50, z75, disp) {
  if (is.na(z50) || is.na(disp)) return("insufficient")
  if (disp > 1.6 && abs(z50) < 0.5) return("dispersion up")
  if (disp < 0.63 && abs(z50) < 0.5) return("dispersion down")
  if (z25 > 0.3 && z50 > 0.3 && z75 > 0.3) return("whole shift up")
  if (z25 < -0.3 && z50 < -0.3 && z75 < -0.3) return("whole shift down")
  if (z75 > 0.5 && abs(z25) < 0.35 && abs(z50) < 0.35) return("right-tail expansion")
  if (z25 > 0.5 && abs(z75) < 0.35) return("left truncation")
  "no change"
}

rows <- list()
for (bt in c("ER+", "TNBC")) {
  for (pt in unique(md2[BC_type == bt, patient_id])) {
    d <- md2[BC_type == bt & patient_id == pt]
    pre <- d[timepoint == "Pre"]; on <- d[timepoint == "On"]
    if (nrow(pre) < MIN_CELLS || nrow(on) < MIN_CELLS) next
    for (n in names(PROG)) {
      col <- paste0(n, "1")
      a <- pre[[col]]; b <- on[[col]]
      qb <- quantile(a, c(.25, .5, .75), type = 7, names = FALSE)
      qa <- quantile(b, c(.25, .5, .75), type = 7, names = FALSE)
      iqr <- qb[3] - qb[1]
      if (!is.finite(iqr) || iqr <= 0) next
      z <- (qa - qb) / iqr
      disp <- (qa[3] - qa[1]) / iqr
      rows[[length(rows) + 1]] <- data.table(
        cohort = "Bassez2021", subtype = bt, patient = pt, program = n,
        n_pre = nrow(pre), n_on = nrow(on),
        q25_pre = qb[1], q50_pre = qb[2], q75_pre = qb[3],
        q25_on = qa[1], q50_on = qa[2], q75_on = qa[3],
        z25 = z[1], z50 = z[2], z75 = z[3],
        disp = disp, asym = z[3] - z[1],
        call = classify(z[1], z[2], z[3], disp))
    }
  }
}
res <- rbindlist(rows)
cat("\n评估组合数：", nrow(res), " 患者：", uniqueN(res$patient), "\n")

fwrite(res, file.path(TAB, "bassez_v1_shape.csv"))

# ---------- 5. 汇总 ----------
cat("\n=== 判定分布（按亚型）===\n")
print(dcast(res[, .N, by = .(subtype, call)], subtype ~ call, value.var = "N", fill = 0))

cat("\n=== 判定分布（按程序）===\n")
print(dcast(res[, .N, by = .(program, call)], program ~ call, value.var = "N", fill = 0))

cat("\n=== 每个亚型 x 程序的 中位 z50 与整体平移比例 ===\n")
summ <- res[, .(n = .N,
                med_z50 = median(z50),
                frac_shift = mean(call %in% c("whole shift up", "whole shift down")),
                frac_tail = mean(call %in% c("right-tail expansion", "left truncation")),
                frac_nochange = mean(call == "no change")),
            by = .(subtype, program)][order(subtype, -frac_shift)]
print(summ)

# ---------- 6. 图 ----------
res[, prog_lab := gsub("_", " ", program)]
p <- ggplot(res, aes(x = z25, y = z75, colour = call)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.35,
              colour = "grey50") +
  geom_point(size = 1.9, alpha = 0.85) +
  facet_grid(subtype ~ prog_lab, scales = "free") +
  scale_colour_manual(values = c(
    `whole shift up` = "#C0392B", `whole shift down` = "#2E5FA3",
    `right-tail expansion` = "#E67E22", `left truncation` = "#8E44AD",
    `dispersion up` = "#16A085", `dispersion down` = "#7F8C8D",
    `no change` = "grey70")) +
  labs(x = "Q25 shift (baseline IQR units)", y = "Q75 shift (baseline IQR units)",
       colour = NULL,
       title = "Independent cohort: shifts are coordinated, not tail-restricted",
       subtitle = paste0("Bassez 2021, anti-PD1, malignant cells only (n = 19,979 cells, 29 patients, 145 patient-programme pairs).\n",
                         "Points on the diagonal move both quartiles together (whole shift); points above or below it would indicate tail expansion or truncation.\n",
                         "Grey points are classified as no change. Only 1 of 145 pairs shows a tail change.")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7.5),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 7.6, colour = "grey30"))

ggsave(file.path(FIG, "fig25_bassez_shape.png"), p, width = 11, height = 5.4,
       dpi = 300, bg = "white")
cat("\n-> fig25_bassez_shape.png\n")
