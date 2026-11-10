# fig31：Fig 6 外部验证（Bassez 2021，anti-PD1）
#
# 一句话：Fig 5 的"整体平移 + 上端动得更多"不是本队列的产物——在一个治疗类型、
#         时间尺度、亚型构成都不同的独立队列里，同一条梯度再现（r10 < 1 < r90）。
#
# A  尾部超额地图：两队列同坐标叠加；紫色参照 = Bassez 逐组合选择反事实（与 Fig 5D 同口径）
# B  归一化梯度：本队列实测 / Bassez 实测 / Bassez 选择反事实
# C  逐组合的 r10 与 r90，全部落在 r10 < 1 < r90 象限
# D  分层梯度：ER+ 与 TNBC 各自成线，不是亚型驱动的
# E  位移幅度 ECDF：两队列幅度差数倍，但形状判据不依赖幅度
#
# 队列不可比之处（写进图注，不能省）：
#   nc2025  多周期化疗 ± ICI，HR+/HER2-，活检随访以周计
#   Bassez  单剂 anti-PD1，9±2 天后手术，含 ER+/TNBC，随访以天计
#   两队列的程序定义也不同构（cNMF 元程序 vs 文献基因集，Jaccard 0.024），
#   因此本图比较的是**位移的形状**，不是位移的大小，也不是某个程序的生物学含义。
#
# 缓存：Bassez 细胞级打分写在 results/tables/bassez_v2_cell_scores.csv，
#       存在时跳过 Seurat 重算，便于反复调图。

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

tab <- "results/tables"; fig <- "results/figures"
CACHE <- file.path(tab, "bassez_v2_cell_scores.csv")
MIN_CELLS <- 50
QS <- c(.10, .25, .50, .75, .90)

theme_set(theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 10.2),
                  plot.subtitle = element_text(size = 7.3, colour = "grey30",
                                               lineheight = 1.06),
                  panel.grid.minor = element_blank(),
                  plot.title.position = "plot"))

COL_B <- "#7F8C8D"; COL_T <- "#C0392B"; COL_U <- "#2E5FA3"
COL_E <- "#B03A2E"; COL_N <- "#2874A6"
COL_CF <- "#6C3483"   # 选择反事实参照（与 Fig 5D 同方法；避开本图已用的蓝/红色点）

# ---------- 1. 细胞级打分（缓存） ----------
if (file.exists(CACHE)) {
  cat("读缓存：", CACHE, "\n")
  cs <- fread(CACHE)
} else {
  suppressPackageStartupMessages(library(Seurat))
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
  cat("读取 Bassez 对象（596 MB）...\n")
  obj <- readRDS("data/raw/bassez/bassez_downsampled.rds")
  md <- as.data.table(obj@meta.data); md[, cell := rownames(obj@meta.data)]
  keep <- md[cellType == "Cancer_cell" & BC_type %in% c("ER+", "TNBC")]
  cat("恶性细胞：", nrow(keep), "\n")
  PROG <- lapply(PROG, function(g) intersect(g, rownames(obj)))
  sub <- subset(obj, cells = keep$cell); rm(obj); gc()
  for (n in names(PROG))
    sub <- AddModuleScore(sub, features = list(PROG[[n]]), name = n,
                          ctrl = 30, seed = 20261108)
  md2 <- as.data.table(sub@meta.data); md2[, cell := rownames(sub@meta.data)]
  cs <- data.table(cell = md2$cell, patient = as.character(md2$patient_id),
                   subtype = md2$BC_type, timepoint = md2$timepoint)
  for (n in names(PROG)) cs[, (n) := md2[[paste0(n, "1")]]]
  rm(sub, md2); gc()
  fwrite(cs, CACHE)
  cat("-> ", CACHE, "\n")
}
PROGS <- c("ER_signal", "Proliferation", "IFN_response", "EMT_matrix",
           "Luminal_secretory")
cs[, patient := as.character(patient)]
cat("细胞 ", nrow(cs), " 患者 ", uniqueN(cs$patient), " 亚型 ",
    paste(unique(cs$subtype), collapse = "/"), "\n")

# ---------- 2. 分位剖面 + 选择反事实 ----------
rows <- list()
for (pt in unique(cs$patient)) {
  sub1 <- cs[patient == pt]
  for (sg in PROGS) {
    pre <- sub1[timepoint == "Pre"][[sg]]
    on <- sub1[timepoint == "On"][[sg]]
    if (length(pre) < MIN_CELLS || length(on) < MIN_CELLS) next
    qb <- quantile(pre, QS, type = 7, names = FALSE)
    qa <- quantile(on, QS, type = 7, names = FALSE)
    iqr <- qb[3] - qb[1]
    if (!is.finite(iqr) || iqr <= 0) next
    d50 <- (qa[3] - qb[3]) / iqr
    # 选择反事实：按位移方向只保留 20% 的基线细胞
    sel <- if (d50 > 0) pre[pre >= quantile(pre, 0.80)] else
                        pre[pre <= quantile(pre, 0.20)]
    qs <- quantile(sel, QS, type = 7, names = FALSE)
    r <- list(patient = pt, subtype = unique(sub1$subtype)[1], prog = sg,
              n_pre = length(pre), n_on = length(on), iqr = iqr,
              q25_base = qb[2], q75_base = qb[4],
              q25_on = qa[2], q75_on = qa[4])
    for (k in seq_along(QS)) {
      r[[paste0("d", QS[k] * 100)]] <- (qa[k] - qb[k]) / iqr
      r[[paste0("s", QS[k] * 100)]] <- (qs[k] - qb[k]) / iqr
    }
    rows[[length(rows) + 1]] <- as.data.table(r)
  }
}
P <- rbindlist(rows)
cat("组合：", nrow(P), " 患者：", uniqueN(P$patient), "\n")

P[, `:=`(z25 = d25, z50 = d50, z75 = d75)]
P[, disp := (q75_on - q25_on) / iqr]
P[, asym := z75 - z25]
P[, u := d75 - d50][, l := d25 - d50]

# 判定：与 fig30 同一套 call_v2 规则（展宽看幅度是否超过位移，不单看中位）
P[, call := fifelse(disp > 1.6 & abs(z50) < asym / 2, "dispersion up",
             fifelse(disp < 0.63 & abs(z50) < 0.5, "dispersion down",
             fifelse(z25 > 0.3 & z50 > 0.3 & z75 > 0.3, "whole shift up",
             fifelse(z25 < -0.3 & z50 < -0.3 & z75 < -0.3, "whole shift down",
             fifelse(z75 > 0.5 & abs(z25) < 0.35 & abs(z50) < 0.35,
                     "right-tail expansion",
             fifelse(z25 > 0.5 & abs(z75) < 0.35, "left truncation",
                     "no change"))))))]

# 归一化梯度
for (k in QS) {
  P[, (sprintf("r%d", k * 100)) := get(sprintf("d%d", k * 100)) / d50]
  P[, (sprintf("rs%d", k * 100)) := get(sprintf("s%d", k * 100)) / s50]
}

BIG <- P[d50 > 0.5]          # 只取上移：下移组合的分位比解读会翻转
cat("上移组合（d50>0.5）：", nrow(BIG), " 上移 ", sum(BIG$d50 > 0),
    " 下移 ", sum(BIG$d50 < 0), "\n")
cat("\n判定分布：\n"); print(P[, .N, by = call][order(-N)])
cat("\n尾部变化：", P[call %in% c("right-tail expansion", "left truncation"), .N],
    "/", nrow(P), "\n")
cat("下尾超额最大 ", round(max(P$l), 2), " 上尾最大 ", round(max(P$u), 2), "\n")

cat("\n== Bassez 归一化梯度（实质位移组合）==\n")
print(BIG[, .(p = QS * 100,
              obs = sapply(QS, function(k) median(get(sprintf("r%d", k * 100)))),
              sel = sapply(QS, function(k) median(get(sprintf("rs%d", k * 100)))))])
cat("r10<1<r90 的组合：", BIG[r10 < 1 & r90 > 1, .N], "/", nrow(BIG), "\n")
cat("逐患者：\n")
print(BIG[, .(n = .N, r10 = round(median(r10), 2), r90 = round(median(r90), 2)),
          by = .(subtype, patient)][order(subtype, -r90)])
cat("\n亚型分层：\n")
print(BIG[, .(n = .N, r10 = median(r10), r90 = median(r90)), by = subtype])

fwrite(P, file.path(tab, "bassez_v2_profile.csv"))
fwrite(BIG[, .(patient, subtype, prog, d10 = round(d10, 3), d50 = round(d50, 3),
               d90 = round(d90, 3), r10 = round(r10, 3), r90 = round(r90, 3),
               rs10 = round(rs10, 3), rs90 = round(rs90, 3))],
       file.path(tab, "bassez_v2_gradient.csv"))

# ---------- 3. 定向、与自队列对照 ----------
# 与 fig30 D 同一套定义：lead = 位移方向的远端，trail = 另一端，两者都乘 sign(d50)
P[, `:=`(ex_lead  = fifelse(d50 > 0, d75 - d50, d25 - d50) * sign(d50),
         ex_trail = fifelse(d50 > 0, d25 - d50, d75 - d50) * sign(d50))]
BIG  <- P[d50 > 0.5]        # 重新取一次，带上定向列
DOWN <- P[d50 < -0.5]
cat("下沉组合:", nrow(DOWN), " 后缘超额 >0.5 的:", sum(DOWN$ex_trail > 0.5), "\n")

g5  <- fread(file.path(tab, "nc2025_fig5_gradient.csv"))
sh5 <- fread(file.path(tab, "nc2025_fig5_shape_calls.csv"))
for (v in c("z50", "ex_lead", "ex_trail")) sh5[, (v) := as.numeric(get(v))]
sh5u <- sh5[z50 > 0.5]

nc_med <- sapply(QS, function(k) median(g5[[sprintf("r%d", k * 100)]]))
bs_med <- sapply(QS, function(k) median(BIG[[sprintf("r%d", k * 100)]]))
bs_sel <- sapply(QS, function(k) median(BIG[[sprintf("rs%d", k * 100)]]))
cat("\nnc2025 实测:", paste(round(nc_med, 2), collapse = " "), "\n")
cat("Bassez 实测:", paste(round(bs_med, 2), collapse = " "), "\n")
cat("Bassez 选择:", paste(round(bs_sel, 2), collapse = " "), "\n")

# ---------- 4. 图 ----------
# ---- 与 Fig 5D 同口径：删除旧的 "selection would land here" 红矩形 ----
# 旧矩形 x 起点 0.55 按 Q10 直觉画，与 ex_trail 的 Q25 口径不符（Bassez 反事实
# 中位只有 +0.19）。现直读 R/39 落盘的 21 个 Bassez 上移组合各自的 top-20%
# 反事实：淡紫带 = 10th-90th percentile，紫色虚线 = 中位。
# 下移侧 4 个越过 0.5 IQR 的点改用红色空心圈标出；0.5 是 S7 正文的**描述性
# 削顶判据**，不是选择反事实。这样保留"这条轴能记录削顶、只是上移侧不发生"
# 的论证，却不再依赖一个错误的固定矩形。
cf6 <- fread(file.path(tab, "selection_counterfactual_paired.csv"))
cf6 <- cf6[cohort == "Bassez" & z50 > 0.5]
cf6_med  <- median(cf6$cf)
cf6_band <- as.numeric(quantile(cf6$cf, c(.10, .90)))
n_cf6    <- nrow(cf6)
DOWN[, trunc_ev := ex_trail > 0.5]
n_trunc6 <- sum(DOWN$trunc_ev)
stopifnot(n_cf6 == nrow(BIG), sum(BIG$ex_trail < cf6_med) == n_cf6)
cat("\n[6A] Bassez 反事实中位 ", sprintf("%.3f", cf6_med),
    " 10-90% ", sprintf("%.3f", cf6_band[1]), " ~ ", sprintf("%.3f", cf6_band[2]),
    " 上移低于中位 ", n_cf6, "/", n_cf6,
    " 下移越过 0.5 判据 ", n_trunc6, "/", nrow(DOWN), "\n")

pA <- ggplot() +
  annotate("rect", xmin = cf6_band[1], xmax = cf6_band[2],
           ymin = -1.0, ymax = 4.1, fill = COL_CF, alpha = 0.16) +
  geom_vline(xintercept = cf6_med, colour = COL_CF,
             linetype = "22", linewidth = 0.55) +
  annotate("text", x = cf6_band[2] + 0.09, y = 2.45,
           label = paste0("selection\ncounterfactual\n+",
                          sprintf("%.2f", cf6_med), " IQR"),
           size = 2.05, colour = COL_CF, hjust = 0, vjust = 0.5,
           lineheight = 0.95) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(data = DOWN, aes(x = ex_trail, y = ex_lead), size = 1.3,
             alpha = 0.55, colour = "#9DB4C8", shape = 1, stroke = 0.55) +
  geom_point(data = BIG, aes(x = ex_trail, y = ex_lead), size = 1.55,
             alpha = 0.8, colour = COL_N) +
  geom_point(data = sh5u, aes(x = ex_trail, y = ex_lead), size = 1.9,
             alpha = 0.85, colour = COL_T, shape = 17) +
  geom_point(data = DOWN[trunc_ev == TRUE], aes(x = ex_trail, y = ex_lead),
             shape = 1, size = 3.3, colour = "#8B1A1A", stroke = 0.7,
             inherit.aes = FALSE) +
  annotate("text", x = 1.95, y = 3.92, label = "nc2025 upward (16)", size = 2.15,
           colour = COL_T, hjust = 1) +
  annotate("text", x = 1.95, y = 3.58, label = "Bassez upward (21)", size = 2.15,
           colour = COL_N, hjust = 1) +
  annotate("text", x = 1.95, y = 3.24, label = "Bassez downward (17)", size = 2.15,
           colour = "#7E96AC", hjust = 1) +
  scale_x_continuous(limits = c(-4.3, 2.05), breaks = seq(-4, 2, 1)) +
  scale_y_continuous(limits = c(-1.0, 4.1), breaks = c(0, 1, 2, 3, 4)) +
  labs(x = "Trailing-tail excess (IQR)", y = "Leading-tail excess (IQR)",
       title = "A  Both cohorts move at the leading side",
       subtitle = paste0(
         "Axes run along each pair's shift. Purple = ", n_cf6,
         " Bassez pair-specific\n",
         "selection counterfactuals (band = 10th-90th percentile, dashed\n",
         "line = median +", sprintf("%.2f", cf6_med),
         " IQR); no Bassez upward pair reaches it. Open red circles\n",
         "= ", n_trunc6, " of ", nrow(DOWN), " downward pairs above 0.5 IQR: ",
         "truncation\nshows on this axis, but not on the upward side.")) +
  theme(panel.grid.major = element_blank())

pB <- {
  grd <- data.table(p = rep(QS * 100, 3),
                    v = c(nc_med, bs_med, bs_sel),
                    set = rep(c("nc2025 observed", "Bassez observed",
                                "Bassez if selecting"), each = 5))
  grd[, set := factor(set, levels = c("nc2025 observed", "Bassez observed",
                                      "Bassez if selecting"))]
  ggplot(grd, aes(x = p, y = v, colour = set, linetype = set)) +
    geom_hline(yintercept = 1, linewidth = 0.3, linetype = "22", colour = "grey55") +
    geom_line(linewidth = 0.85) + geom_point(size = 1.7) +
    scale_colour_manual(values = c("nc2025 observed" = COL_T,
                                   "Bassez observed" = COL_N,
                                   "Bassez if selecting" = "grey45"), name = NULL) +
    scale_linetype_manual(values = c("22", "solid", "solid"), name = NULL) +
    scale_x_continuous(breaks = QS * 100,
                       labels = c("Q10", "Q25", "Q50", "Q75", "Q90")) +
    labs(x = "Quantile of the programme score",
         y = "Shift relative to\nthe median shift",
         title = "B  The same gradient, the same wrong direction",
         subtitle = paste0("Observed rises from 0.63 to 1.29 of the median in our ",
                           "cohort\nand from 0.27 to 1.44 in Bassez. Selecting the top ",
                           "fifth of baseline\ncells does the opposite in both.")) +
    theme(legend.position = "bottom", legend.key.size = unit(0.32, "cm"),
          legend.text = element_text(size = 6.6),
          legend.margin = margin(-4, 0, 0, 0)) +
    guides(colour = guide_legend(nrow = 2))
}

pC <- ggplot(BIG, aes(x = r10, y = r90)) +
  annotate("rect", xmin = -Inf, xmax = 1, ymin = 1, ymax = Inf,
           fill = "#C0392B", alpha = 0.07) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey60") +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey60") +
  geom_point(aes(colour = subtype), size = 2.0, alpha = 0.9) +
  annotate("text", x = -0.05, y = 2.62, label = sprintf("%d of %d pairs",
           BIG[r10 < 1 & r90 > 1, .N], nrow(BIG)),
           size = 2.3, colour = "#8B1A1A", hjust = 0, vjust = 1) +
  scale_colour_manual(values = c(`ER+` = COL_E, TNBC = COL_N), name = NULL) +
  annotate("text", x = 2.25, y = 2.95, hjust = 1, vjust = 1, size = 2.1,
           colour = "grey40",
           label = sprintf("%d pairs above\nthe axis (Q90 up to %.1f)",
                           nrow(BIG[r90 > 3]), max(BIG$r90))) +
  scale_x_continuous(limits = c(-0.15, 2.3)) +
  scale_y_continuous(limits = c(0, 3.05)) +
  labs(x = "Q10 shift / median shift", y = "Q90 shift / median shift",
       title = "C  One point, one upward pair",
       subtitle = paste0("Shaded quadrant = bottom moves less than the median while ",
                         "the top\nmoves more. ", BIG[r10 < 1 & r90 > 1, .N], " of ",
                         nrow(BIG), " pairs fall inside; the rest sit\nclose to the ",
                         "diagonal, none in the mirror quadrant.")) +
  theme(legend.position = "bottom", legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 6.6),
        legend.margin = margin(-4, 0, 0, 0)) +
  guides(colour = guide_legend(nrow = 1))

grd2 <- rbindlist(lapply(c("all", "ER+", "TNBC"), function(g) {
  d <- if (g == "all") BIG else BIG[subtype == g]
  if (!nrow(d)) return(NULL)
  data.table(p = QS * 100,
             v = sapply(QS, function(k) median(d[[sprintf("r%d", k * 100)]])),
             grp = sprintf("%s (%d)", if (g == "all") "all pairs" else g, nrow(d)))
}))
grd2[, grp := factor(grp, levels = unique(grp))]
pD <- ggplot(grd2, aes(x = p, y = v, colour = grp)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "22", colour = "grey55") +
  geom_line(aes(linetype = grp), linewidth = 0.8) + geom_point(size = 1.6) +
  scale_colour_manual(values = setNames(c("grey20", COL_E, COL_N), levels(grd2$grp)),
                      name = NULL) +
  scale_x_continuous(breaks = QS * 100,
                     labels = c("Q10", "Q25", "Q50", "Q75", "Q90")) +
  labs(x = "Quantile of the programme score",
       y = "Shift relative to\nthe median shift",
       title = "D  Not driven by one subtype",
       subtitle = paste0("Bassez mixes ER+ and TNBC lesions given the same anti-PD1. ",
                         "Both\nsubtypes give a rising gradient, so the shape is not an ",
                         "artefact of\nhormone-receptor biology.")) +
  theme(legend.position = "bottom", legend.key.size = unit(0.28, "cm"),
        legend.text = element_text(size = 6.1),
        legend.margin = margin(-4, 0, 0, 0)) +
  guides(colour = guide_legend(nrow = 2), linetype = "none")

XCAP <- 3.2
ecdf_d <- rbind(data.table(cohort = "nc2025", v = abs(sh5$z50)),
                data.table(cohort = "Bassez", v = abs(P$d50)))[is.finite(v)]
n_cap <- ecdf_d[v > XCAP, .N]
ecdf_d[, v := pmin(v, XCAP)]
pE <- ggplot(ecdf_d, aes(x = v, colour = cohort)) +
  stat_ecdf(linewidth = 0.8, geom = "step") +
  geom_vline(xintercept = 0.5, linewidth = 0.3, linetype = "22", colour = "grey50") +
  annotate("text", x = 0.55, y = 0.09, label = "effect threshold", size = 2.15,
           colour = "grey35", hjust = 0) +
  scale_colour_manual(values = c(nc2025 = COL_T, Bassez = COL_N), name = NULL) +
  scale_x_continuous(limits = c(0, XCAP + 0.05)) +
  labs(x = "|median shift| (baseline IQR)", y = "Cumulative share of pairs",
       title = "E  Sizes differ, shapes do not",
       subtitle = sprintf(paste0("Median |shift| %.2f IQR in our cohort versus %.2f in ",
                                 "Bassez - a different drug, schedule and subtype mix, ",
                                 "with %.0f%% of Bassez pairs\nbelow the 0.5 effect ",
                                 "threshold. The axis stops at 3.2 IQR, where ",
                                 "%d pairs are censored.\nWhat replicates is the shape ",
                                 "of the shift, not its size."),
                          median(abs(sh5$z50)), median(abs(P$d50)),
                          100 * mean(abs(P$d50) < 0.5), n_cap)) +
  theme(legend.position = c(0.80, 0.42), legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 6.8), legend.background = element_blank())

## 图级标题与图内解释段于 2026-11-10 移除（改由图注承担），与其他主图一致。
design <- "
AABB
CCDD
EEEE
"
f6 <- wrap_plots(A = pA, B = pB, C = pC, D = pD, E = pE,
                 design = design) +
  plot_layout(heights = c(1.00, 1.00, 0.72))

ggsave(file.path(fig, "fig31_fig6_external.png"), f6, width = 11.4, height = 8.9,
       dpi = 300, bg = "white")
cat("\n-> fig31_fig6_external.png\n")
