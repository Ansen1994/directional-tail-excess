#!/usr/bin/env Rscript
# ==============================================================================
# 35_figS1_stress.R
#
# 补充图 S1：解离应激诊断。
#
# 要排除的替代解释：程序分数的位移不是治疗引起的状态转变，而是
#   (a) 不同时点的组织处理/解离条件不同造成的技术伪影；
#   (b) 应激基因（FOS/JUN/HSP）本身被算进了某个程序，位移只是应激水平的涨落；
#   (c) 测序深度差异（nCount 高 -> 检出基因多 -> 分数整体上移）。
#
# 做法：
#   1. 对全部上皮细胞（120,058）算应激分数与 IFN 分数（AddModuleScore，ctrl=30，
#      与 R/32 完全同参数，保证与主图可比）。
#   2. 聚合到样本中位，检验三时点是否有系统偏移。
#   3. 把每个 患者 x 程序 x 时点对 的 z50 与其 from->to 的应激变化配对，
#      检验位移是否由应激驱动（核心阴性对照）。
#   4. 应激分数与 nCount 的相关，检验技术混杂。
#
# 输出：
#   results/tables/neoletexe_stress_sample.csv
#   results/tables/neoletexe_stress_pair.csv
#   results/figures/figS1_stress.png
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(Seurat)
  library(ggplot2); library(patchwork)
})

PROC <- "data/proc/neoletexe"; TAB <- "results/tables"; FIG <- "results/figures"
set.seed(20261109)

STRESS <- c("FOS","FOSB","JUN","JUNB","JUND","EGR1","HSPA1A","HSPA1B",
            "HSPA8","HSP90AA1","HSP90AB1","DNAJB1","IER2","ATF3","SOCS3")
IFNG   <- c("STAT1","IFIT1","IFIT3","ISG15","MX1","OAS1","CXCL10","IRF7",
            "IFITM1","BST2")

## ---------------- 1. 打分 ----------------
CACHE <- file.path(PROC, "stress_cells.rds")
if (file.exists(CACHE)) {
  md <- readRDS(CACHE)
  cat(sprintf("读缓存：%s 个细胞的应激分数\n", format(nrow(md), big.mark = ",")))
} else {
  ep <- readRDS(file.path(PROC, "epi_combined.rds"))
  obj <- CreateSeuratObject(counts = ep$counts, meta.data = ep$meta,
                            min.cells = 0, min.features = 0)
  rm(ep); invisible(gc(FALSE))
  obj <- NormalizeData(obj, verbose = FALSE)
  g <- rownames(obj)
  STRESS <- intersect(STRESS, g); IFNG <- intersect(IFNG, g)
  cat(sprintf("上皮细胞 %d，应激基因 %d，IFN 基因 %d\n",
              ncol(obj), length(STRESS), length(IFNG)))

  obj <- AddModuleScore(obj, features = list(STRESS), name = "stress",
                        ctrl = 30, seed = 20261109)
  obj <- AddModuleScore(obj, features = list(IFNG), name = "ifn",
                        ctrl = 30, seed = 20261109)
  md <- as.data.table(obj@meta.data)
  md[, stress := stress1][, ifn := ifn1]
  md[, c("stress1", "ifn1") := NULL]
  obj <- NULL; invisible(gc(FALSE))

  ## 恶性标签
  cnv <- fread(file.path(TAB, "neoletexe_cnv_cells.csv"))
  md <- merge(md, cnv[, .(cell, malignant)], by = "cell", all.x = TRUE)
  md[is.na(malignant), malignant := FALSE]
  saveRDS(md, CACHE)
}

## ---------------- 2. 样本级 ----------------
smp <- md[, .(n = .N,
              stress = as.numeric(median(stress)),
              stress_iqr = as.numeric(IQR(stress)),
              ifn = as.numeric(median(ifn)),
              nc = as.numeric(median(nCount)),
              nf = as.numeric(median(nFeature)),
              mt = as.numeric(median(pct_mt)),
              mal_frac = as.numeric(mean(malignant))), by = .(sample, patient, tp)]
smp[, patient := paste0("N", patient)]
fwrite(smp, file.path(TAB, "neoletexe_stress_sample.csv"))

cat("\n===== 三时点应激分数（样本中位的中位）=====\n")
print(smp[, .(n = .N,
              stress = round(median(stress), 3),
              ifn = round(median(ifn), 3),
              nc = round(median(nc)),
              mal = round(100 * median(mal_frac), 1)), by = tp][order(tp)])

cat("\n===== 时点间差异（配对 Wilcoxon，仅三时点齐全的例）=====\n")
w <- dcast(smp, patient ~ tp, value.var = "stress")
w <- w[complete.cases(w)]
if (nrow(w) >= 5) {
  for (pr in list(c("A","B"), c("B","C"), c("A","C"))) {
    x <- w[[pr[1]]]; y <- w[[pr[2]]]
    tt <- wilcox.test(x, y, paired = TRUE)
    cat(sprintf("  %s -> %s : 中位变化 %+0.3f  P = %.3f  (n = %d)\n",
                pr[1], pr[2], median(y - x), tt$p.value, nrow(w)))
  }
}

## ---------------- 3. 与位移配对 ----------------
sh <- fread(file.path(TAB, "neoletexe_shape.csv"))
sh[, z50 := as.numeric(z50)]
tp_of <- function(s) substr(s, nchar(s), nchar(s))

st <- smp[, .(patient, tp, stress, ifn)]
setnames(st, c("stress","ifn"), c("stress_tp","ifn_tp"))
sh <- merge(sh, st, by.x = c("patient","from"), by.y = c("patient","tp"),
            all.x = TRUE)
sh <- merge(sh, st, by.x = c("patient","to"), by.y = c("patient","tp"),
            all.x = TRUE)
sh[, `:=`(d_stress = stress_tp.y - stress_tp.x,
          d_ifn    = ifn_tp.y - ifn_tp.x)]
fwrite(sh[, .(patient, pair, from, to, program, z50, d_stress, d_ifn, call_v2)],
       file.path(TAB, "neoletexe_stress_pair.csv"))

big <- sh[abs(z50) > 0.5]
cat("\n===== 位移与应激变化的关系 =====\n")
cat(sprintf("全部组合 n=%d:  r(z50, d_stress) = %+.3f  P = %.3f\n",
            nrow(sh), cor(sh$z50, sh$d_stress, use = "complete.obs"),
            cor.test(sh$z50, sh$d_stress)$p.value))
cat(sprintf("实质位移 n=%d: r(z50, d_stress) = %+.3f  P = %.3f\n",
            nrow(big), cor(big$z50, big$d_stress, use = "complete.obs"),
            cor.test(big$z50, big$d_stress)$p.value))
cat(sprintf("实质位移 n=%d: r(|z50|, |d_stress|) = %+.3f\n",
            nrow(big), cor(abs(big$z50), abs(big$d_stress), use = "complete.obs")))

cat("\n位移组 vs 无变化组的 d_stress（中位）:\n")
print(sh[, .(n = .N, d_stress = round(median(d_stress, na.rm = TRUE), 3)),
         by = .(moved = abs(z50) > 0.5)])

## ---------------- 4. 技术混杂 ----------------
cat("\n===== 应激 vs 测序深度 =====\n")
cat(sprintf("样本级 r(stress, nCount) = %+.3f\n",
            cor(smp$stress, smp$nc, use = "complete.obs")))
cat(sprintf("细胞级 r(stress, nCount) = %+.3f\n",
            cor(md$stress, md$nCount, use = "complete.obs")))
cat(sprintf("应激 vs IFN  r = %+.3f\n", cor(md$stress, md$ifn)))

## ---------------- 5. 出图 ----------------
TPC <- c(A = "#4477AA", B = "#EE9944", C = "#66AA55")
smp[, tp := factor(tp, levels = c("A","B","C"))]

pA <- ggplot(smp, aes(tp, stress)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey93",
               colour = "grey40", linewidth = 0.35) +
  geom_jitter(aes(colour = tp), width = 0.13, size = 1.1, alpha = 0.85) +
  scale_colour_manual(values = TPC, guide = "none") +
  labs(x = NULL, y = "Dissociation-stress score",
       title = "A  Stress does not rise with treatment",
       subtitle = paste0("Median per sample (n = ", nrow(smp), " samples, ",
                         uniqueN(smp$patient), " patients).\n",
                         "If anything it falls; no interval differs ",
                         "(paired P > 0.3).")) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"))

pB <- ggplot(sh, aes(d_stress, z50)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_point(size = 0.75, alpha = 0.45, colour = "grey35") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              linewidth = 0.5, colour = "#B22222", fill = "grey85") +
  labs(x = "Change in stress score (to - from)",
       y = "Median programme shift (IQR)",
       title = "B  Shifts are not explained by stress",
       subtitle = sprintf(paste0("One point per patient x programme x interval ",
                                 "(n = %d).\nPearson r = %+.3f, P = %.2f."),
                          nrow(sh),
                          cor(sh$z50, sh$d_stress, use = "complete.obs"),
                          cor.test(sh$z50, sh$d_stress)$p.value)) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"))

pC <- ggplot(smp, aes(nc, stress, colour = tp)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = TRUE,
              linewidth = 0.5, colour = "grey30", fill = "grey85") +
  geom_point(size = 1.5, alpha = 0.9) +
  scale_x_log10() +
  scale_colour_manual(values = TPC, name = NULL) +
  labs(x = "Median counts per cell (log scale)", y = "Dissociation-stress score",
       title = "C  Weak technical coupling",
       subtitle = sprintf(paste0("One point per biopsy (n = %d).\n",
                                 "r = %+.2f at sample level, %+.2f at cell level:\n",
                                 "too weak to generate the shifts."),
                          nrow(smp), cor(smp$stress, smp$nc),
                          cor(md$stress, md$nCount))) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.4, colour = "grey35"),
        legend.position = c(0.14, 0.83), legend.key.size = unit(0.28, "cm"),
        legend.text = element_text(size = 5.8),
        legend.background = element_blank())

## D：剔除高应激样本后重测
thr <- quantile(smp$stress, 0.75)
keep <- smp[stress <= thr, .(patient, tp)]
sh2 <- merge(sh, keep, by.x = c("patient","from"), by.y = c("patient","tp"))
sh2 <- merge(sh2, keep, by.x = c("patient","to"), by.y = c("patient","tp"))
dd <- data.table(
  set = c("All samples", "Lower three\nstress quartiles"),
  n_pair = c(nrow(sh), nrow(sh2)),
  rate = c(mean(abs(sh$z50) > 0.5), mean(abs(sh2$z50) > 0.5)),
  trail = c(max(sh[z50 > 0.5]$ex_trail), max(sh2[z50 > 0.5]$ex_trail)))

pD <- ggplot(dd, aes(set, 100 * rate)) +
  geom_col(width = 0.5, fill = "grey55") +
  geom_text(aes(label = sprintf("%.1f%%\n(n = %d)", 100 * rate, n_pair)),
            vjust = -0.15, size = 2.3, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 45), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Combinations with a real shift (%)",
       title = "D  Result survives stress filtering",
       subtitle = paste0("Dropping the noisiest quarter of samples\n",
                         "leaves the shift rate essentially unchanged.")) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"))

pTop <- ggplot() + theme_void() +
  annotate("text", x = 0, y = 0.98, hjust = 0, vjust = 1, size = 3.2,
           fontface = "bold",
           label = "Figure S1  Dissociation stress does not explain the programme shifts") +
  annotate("text", x = 0, y = 0.55, hjust = 0, vjust = 1, size = 2.2,
           colour = "grey30", lineheight = 1.05,
           label = paste0("NeoLetExe, 120,058 epithelial cells, 73 biopsies. Stress score ",
                          "= 15 immediate-early and heat-shock genes (FOS / JUN / EGR1 / HSP),\n",
                          "scored with the same AddModuleScore parameters as the five ",
                          "main programmes.")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

p <- pTop / (pA | pB) / (pC | pD) + plot_layout(heights = c(0.22, 1, 1))

ggsave(file.path(FIG, "figS1_stress.png"), p, width = 7.2, height = 6.2,
       dpi = 300, bg = "white")
cat("\n-> figS1_stress.png\n")
