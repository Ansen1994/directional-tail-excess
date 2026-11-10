# fig30：Fig 5 分布形状判别（H5：不是挑出一个高分亚群）
#
# 一句话：位移贯穿整条分布，而且上端比下端动得多——这与"挑出预先存在的高分细胞"
#         所必然留下的痕迹方向相反（选择会让下端动得多、上端几乎不动）。
#
# A  21 个实质位移组合的分位位移剖面（Q10..Q90，基线 IQR 为单位）
# B  归一化梯度：实测 vs 选择反事实。实测 0.63 -> 1.29，选择 1.46 -> 0.73
# C  四个代表组合的单细胞密度（干净平移 / 巨额上移 / 分裂）
# D  尾部超额地图：横轴下尾超额、纵轴上尾超额；蓝色参照 = 逐组合选择反事实
# E  逐患者的上/下端梯度：7 例全部 r10 < 1 < r90
#
# 口径：位移一律除以该组合的基线 IQR（与 fig21/fig27/fig28 一致）。
#       细胞级分位来自 nc2025_v3_density.csv（31 个组合，两时点齐全）。
# 排版：图内一律 ASCII；半宽 panel 的 subtitle 每行 <= 84 字符。

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

tab <- "results/tables"; fig <- "results/figures"

theme_set(theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 10.2),
                  plot.subtitle = element_text(size = 7.3, colour = "grey30",
                                               lineheight = 1.06),
                  panel.grid.minor = element_blank(),
                  plot.title.position = "plot"))

COL_B <- "#7F8C8D"; COL_T <- "#C0392B"; COL_U <- "#2E5FA3"; COL_G <- "#7D8B99"

pat_key <- function(x) {          # 幂等：1 / "01" / "P01" 都吃
  y <- suppressWarnings(as.integer(x))
  ifelse(is.na(y), as.character(x), sprintf("P%02d", y))
}

sh <- fread(file.path(tab, "nc2025_v3_shape.csv"))
setnames(sh, "call", "call_v1")   # 否则 j 里裸写 call 会落到 base::call
dn <- fread(file.path(tab, "nc2025_v3_density.csv"))
sh[, patient := pat_key(patient)]
dn[, patient := pat_key(patient)]
for (v in c("z25", "z50", "z75", "asym", "disp", "iqr_base"))
  sh[, (v) := as.numeric(get(v))]
dn[, score := as.numeric(score)]

PATORD <- c("P01", "P09", "P12", "P13", "P16", "P18", "P23")
MPNAME <- c(MP2 = "MP2 Prolif", MP8 = "MP8 Lum-sec", MP7 = "MP7 ER-I",
            MP3 = "MP3 Stromal", MP11 = "MP11 EMT-III", MP6 = "MP6 TGFb",
            MP9 = "MP9 Neural", MP4 = "MP4 IFN", MP5 = "MP5 IEG")

QS <- c(.10, .25, .50, .75, .90)

# ---------- 细胞级分位（31 个组合，两时点） ----------
qu <- dn[, as.list(quantile(score, QS)), by = .(patient, mp, phase)]
setnames(qu, c("patient", "mp", "phase", paste0("q", QS * 100)))
qu <- dcast(qu, patient + mp ~ phase,
            value.var = paste0("q", QS * 100))
P <- merge(qu, sh[, .(patient, mp, z50, iqr_base)], by = c("patient", "mp"))

prof <- rbindlist(lapply(QS, function(p) {
  b <- P[[sprintf("q%d_baseline", round(p * 100))]]
  a <- P[[sprintf("q%d_on-treatment", round(p * 100))]]
  data.table(patient = P$patient, mp = P$mp, iqr = P$iqr_base, z50 = P$z50,
             p = p * 100, d = (a - b) / P$iqr_base)
}))

# ---------- 选择反事实：按位移方向只保留 20% 的基线细胞 ----------
cf <- rbindlist(lapply(seq_len(nrow(P)), function(i) {
  b <- dn[patient == P$patient[i] & mp == P$mp[i] & phase == "baseline", score]
  up <- P$z50[i] > 0
  sel <- if (up) b[b >= quantile(b, 0.80)] else b[b <= quantile(b, 0.20)]
  qb <- quantile(b, QS); qs <- quantile(sel, QS)
  data.table(patient = P$patient[i], mp = P$mp[i], iqr = P$iqr_base[i],
             z50 = P$z50[i], p = QS * 100, d = (qs - qb) / P$iqr_base[i])
}))

# 只保留**上移**组合。两点理由：
# ① 待检验的主张是"程序被诱导升高"，下移组合不在该主张范围内；
# ② 位移为负时，分位比的解读会翻转：下移 + 保留低分细胞得到的是
#    r10 < 1 < r90，与"诱导 + 前缘动得更多"的形态一模一样。把两者混进
#    同一个中位数，等于把相反的结论平均掉。
keep <- prof[, .(patient, mp, z50)][!duplicated(paste(patient, mp))][z50 > 0.5]
prof[, key := paste(patient, mp)]
cf[,   key := paste(patient, mp)]
prof_b <- prof[key %in% paste(keep$patient, keep$mp)]
cf_b   <- cf[key %in% paste(keep$patient, keep$mp)]
prof_b[, r := d / d[p == 50], by = key]
cf_b[,   r := d / d[p == 50], by = key]
cat("实质位移组合:", uniqueN(prof_b$key), " 上移:", uniqueN(prof_b[z50 > 0]$key),
    " 下移:", uniqueN(prof_b[z50 < 0]$key), "\n")

# ================= A 分位位移剖面 =================
prof_b[, dir := ifelse(z50 > 0, "up", "down")]
prof_b[, split := key %in% c("P16 MP8", "P16 MP11")]
medA <- prof_b[, .(d = median(d)), by = p]

pA <- ggplot(prof_b, aes(x = p, y = d, group = key)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_line(data = prof_b[split == FALSE], aes(colour = dir),
            linewidth = 0.42, alpha = 0.60) +
  geom_line(data = prof_b[split == TRUE], colour = "grey20",
            linewidth = 0.55, linetype = "22") +
  geom_line(data = medA, aes(x = p, y = d), colour = "white", linewidth = 2.6,
            inherit.aes = FALSE) +
  geom_line(data = medA, aes(x = p, y = d), colour = "black", linewidth = 1.25,
            inherit.aes = FALSE) +
  geom_point(data = medA, aes(x = p, y = d), colour = "black", size = 1.0,
             inherit.aes = FALSE) +
  scale_colour_manual(values = c(up = COL_T, down = COL_U), guide = "none") +
  scale_x_continuous(breaks = QS * 100,
                     labels = c("Q10", "Q25", "Q50", "Q75", "Q90")) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.6),
                     breaks = c(-4, -2, -1, 0, 1, 2, 4, 8, 16),
                     labels = c("-4", "-2", "-1", "0", "1", "2", "4", "8", "16")) +
  labs(x = "Quantile of the programme score", y = "Shift in quantile\n(baseline IQR)",
       title = "A  Every decile moves; the upper end leads the median in 13 of 16",
       subtitle = paste0("One line per upward-moving patient x programme pair (median ",
                         "shift > 0.5 IQR;\nn = 16). Thick black = median profile. No ",
                         "line is flat, so no part of the\ndistribution stays put. Upper-end ",
                         "lead is significant against label\npermutation in 7 pairs only. ",
                         "Dashed = P16-MP8 and P16-MP11 (see C, D).\nVertical axis is ",
                         "pseudo-log, so shifts up to 16 IQR stay visible.")) +
  theme(panel.grid.major.x = element_blank())

# ================= B 归一化梯度 =================
band <- function(dt) dt[, .(lo = quantile(r, .25), hi = quantile(r, .75),
                            md = median(r)), by = p]
bO <- band(prof_b); bS <- band(cf_b)
bO[, src := "observed"]; bS[, src := "selection"]
bB <- rbind(bO, bS)

pB <- ggplot(bB, aes(x = p, y = md, colour = src, fill = src)) +
  geom_hline(yintercept = 1, colour = "grey60", linetype = "22", linewidth = 0.35) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c(observed = COL_T, selection = "grey45"),
                      labels = c(observed = "Observed",
                                 selection = "If high-scoring cells were selected"),
                      name = NULL) +
  scale_fill_manual(values = c(observed = COL_T, selection = "grey45"),
                    labels = c(observed = "Observed",
                               selection = "If high-scoring cells were selected"),
                    name = NULL) +
  scale_x_continuous(breaks = QS * 100, labels = c("Q10", "Q25", "Q50", "Q75", "Q90")) +
  scale_y_continuous(breaks = c(0.5, 1.0, 1.5, 2.0)) +
  labs(x = "Quantile of the programme score", y = "Shift relative to\nthe median shift",
       title = "B  The gradient runs the wrong way for selection",
       subtitle = paste0("Each shift is divided by that pair's own median shift, so a ",
                         "pure\ntranslation would be flat at 1. Observed rises from 0.63 ",
                         "at Q10 to\n1.29 at Q90; taking only the top-scoring baseline ",
                         "cells predicts the\nopposite, 1.46 falling to 0.73. Upward ",
                         "pairs only (n = 16).")) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "bottom", legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 7.0), legend.margin = margin(-8, 0, 0, 0)) +
  guides(colour = guide_legend(nrow = 1))

# ================= C 密度曲线 =================
CASES <- data.table(
  patient = c("P23", "P18", "P13", "P16"),
  mp      = c("MP7", "MP9", "MP3", "MP8"),
  tag     = c("P23 / MP7  (clean translation, +2.2)",
              "P18 / MP9  (whole curve shifts, +1.3)",
              "P13 / MP3  (upper half outruns, +0.9)",
              "P16 / MP8  (splits apart, +0.7)"))
dc <- rbindlist(lapply(seq_len(nrow(CASES)), function(i) {
  x <- dn[patient == CASES$patient[i] & mp == CASES$mp[i]]
  if (!nrow(x)) return(NULL)
  x[, tag := CASES$tag[i]]
  x
}))
dc[, tag := factor(tag, levels = CASES$tag)]
dc[, phase := factor(phase, levels = c("baseline", "on-treatment"),
                     labels = c("Baseline", "On treatment"))]
dc <- dc[, .SD[score <= quantile(score, 0.995)], by = tag]

pC <- ggplot(dc, aes(x = score, fill = phase, colour = phase)) +
  geom_density(alpha = 0.26, linewidth = 0.45, adjust = 1.15) +
  facet_wrap(~ tag, nrow = 1, scales = "free") +
  scale_fill_manual(values = c(Baseline = COL_B, `On treatment` = COL_T), name = NULL) +
  scale_colour_manual(values = c(Baseline = "#5D6D7E", `On treatment` = COL_T), name = NULL) +
  labs(x = "Programme score per malignant cell", y = "Density",
       title = "C  Cell-level distributions behind the quantile shifts",
       subtitle = paste0("In the first three the whole curve slides right; in P16-MP8 the ",
                         "lower half falls while the\nupper tail climbs. Both combinations with ",
                         "that split shape come from P16, the only patient whose\ntwo samples ",
                         "were sequenced in different batches (Fig. 2).")) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7.6),
        legend.position = "bottom", legend.key.size = unit(0.36, "cm"),
        legend.text = element_text(size = 7.0))

# ================= D 尾部超额地图（按位移方向定向） =================
# 若直接画 u = z75 - z50 与 l = z25 - z50，下移组合的坐标会整体镜像：
# 对它们而言，"下端不动、上端抬起"其实对应的是选择信号的另一侧。
# 因此按位移方向定向：lead = 位移方向的远端（上移取 Q75，下移取 Q25），
# trail = 另一端；两个超额都乘 sign(z50)，正数 = 沿位移方向多走的量。
sh[, `:=`(ex_lead  = fifelse(z50 > 0, z75 - z50, z25 - z50) * sign(z50),
          ex_trail = fifelse(z50 > 0, z25 - z50, z75 - z50) * sign(z50))]
sh[, mat := abs(z50) > 0.5]
sh[, split := patient == "P16" & mp %in% c("MP8", "MP11")]
cat("\n[D] 前缘超额最大 ", round(max(sh$ex_lead), 2),
    " 后缘超额最大 ", round(max(sh$ex_trail), 2),
    " 实质位移组合里后缘超额最大 ",
    round(max(sh[mat == TRUE]$ex_trail), 2), "\n")

# ---- 方案 3：删除旧的 "selection would land here" 红矩形（xmin = 0.55） ----
# 旧矩形按 Q10 直觉画，ex_trail 实际用 Q25/Q75 口径：top-20% 反事实中位仅 +0.24 IQR，
# 16 个组合无一进入矩形。现直接读 R/39 落盘的逐组合反事实，图上只标定 trailing 轴：
#   淡蓝带 = 16 个 pair-specific 反事实的 10th–90th percentile
#   蓝色竖虚线 = 中位。不画二维"选择象限"，避免把幅度当成判据。
cfp <- fread(file.path(tab, "selection_counterfactual_paired.csv"))
cfp <- cfp[cohort == "nc2025" & z50 > 0.5]
cf_med  <- median(cfp$cf)
cf_band <- as.numeric(quantile(cfp$cf, c(.10, .90)))
n_cf    <- nrow(cfp)
n_below <- sum(sh[z50 > 0.5]$ex_trail < cf_med)
stopifnot(n_cf == sum(sh$z50 > 0.5), n_below == n_cf)
cat("\n[D] 反事实参照 n = ", n_cf, "  中位 ", sprintf("%.3f", cf_med),
    "  10-90% ", sprintf("%.3f", cf_band[1]), " ~ ", sprintf("%.3f", cf_band[2]),
    "  观测点低于中位线 ", n_below, "/", n_cf, "\n")

pD <- ggplot(sh, aes(x = ex_trail, y = ex_lead)) +
  annotate("rect", xmin = cf_band[1], xmax = cf_band[2],
           ymin = -0.6, ymax = 4.15, fill = COL_U, alpha = 0.16) +
  geom_vline(xintercept = cf_med, colour = COL_U, linetype = "22", linewidth = 0.55) +
  annotate("text", x = cf_band[2] + 0.09, y = 2.05,
           label = paste0("selection\ncounterfactual\n+", sprintf("%.2f", cf_med)),
           size = 2.05, colour = COL_U, hjust = 0, vjust = 0.5, lineheight = 0.95) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey86", linetype = "22",
              linewidth = 0.3) +
  geom_point(aes(size = abs(z50), colour = mat), alpha = 0.8, shape = 16) +
  geom_point(data = sh[split == TRUE], shape = 1, size = 3.0,
             colour = "grey15", stroke = 0.6) +
  scale_size_continuous(range = c(0.7, 3.4), name = "|median shift|",
                        breaks = c(1, 3, 6),
                        guide = guide_legend(direction = "horizontal")) +
  scale_colour_manual(values = c(`TRUE` = COL_T, `FALSE` = "grey72"), guide = "none") +
  scale_x_continuous(limits = c(-4.3, 2.0), breaks = seq(-4, 1, 1)) +
  scale_y_continuous(limits = c(-0.6, 4.15), breaks = c(0, 1, 2, 3, 4)) +
  labs(x = "Trailing-tail excess (IQR away from the median shift)",
       y = "Leading-tail excess (IQR)",
       title = "D  The extra movement sits on the leading side only",
       subtitle = paste0(
         "Both axes run along each pair's own shift, so the 5 downward\n",
         "pairs are not mirrored. Blue = ", n_cf, " pair-specific selection\n",
         "counterfactuals: band = 10th-90th percentile, dashed line =\n",
         "median +", sprintf("%.2f", cf_med), " IQR; none of the ", n_cf,
         " upward pairs reaches it.\n",
         "Trailing excess tops out at ", sprintf("%.2f", max(sh[z50 > 0.5]$ex_trail)),
         " IQR, against ", round(max(sh$ex_lead), 2), " on the leading side.")) +
  theme(panel.grid.major = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(0.28, "cm"), legend.text = element_text(size = 6.6),
        legend.title = element_text(size = 6.6),
        legend.margin = margin(-8, 0, 0, 0))


# ================= E 逐患者梯度 =================
e <- prof_b[p %in% c(10, 90)]
e[, ratio := r]
e[, q := factor(p, levels = c(10, 90), labels = c("Q10", "Q90"))]
emed <- e[, .(ratio = median(ratio)), by = .(patient, q)]

pE <- ggplot(e, aes(x = ratio, y = factor(patient, levels = rev(PATORD)))) +
  geom_vline(xintercept = 1, colour = "grey40", linetype = "22", linewidth = 0.4) +
  geom_point(aes(colour = q), size = 1.2, alpha = 0.35,
             position = position_jitter(height = 0.10, width = 0, seed = 3)) +
  geom_segment(data = dcast(emed, patient ~ q, value.var = "ratio"),
               aes(x = Q10, xend = Q90,
                   y = factor(patient, levels = rev(PATORD)),
                   yend = factor(patient, levels = rev(PATORD))),
               colour = "grey45", linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(data = emed, aes(colour = q), size = 3.1) +
  scale_colour_manual(values = c(Q10 = COL_U, Q90 = COL_T), name = NULL) +
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2)) +
  labs(x = "Shift relative to the median shift", y = NULL,
       title = "E  The same gradient in every patient",
       subtitle = paste0("Small dots = individual pairs, large dots = the per-patient ",
                         "median.\nAll seven patients sit below 1 at Q10 and above 1 at ",
                         "Q90, so the\npattern is not carried by a single case.")) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom", legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7.2)) +
  guides(colour = guide_legend(nrow = 1))

# ================= 拼图 =================
## 图级标题与图内解释段于 2026-11-10 移除（改由图注承担），与其后主图一致。
design <- "
AABB
CCCC
DDEE
"
f5 <- wrap_plots(A = pA, B = pB, C = pC, D = pD, E = pE, design = design) +
  plot_layout(heights = c(1.00, 0.78, 1.00))

ggsave(file.path(fig, "fig30_fig5_shape.png"), f5,
       width = 11.4, height = 9.6, dpi = 300, bg = "white")
cat("-> fig30_fig5_shape.png\n")

# ================= 导出核对表 =================
out <- dcast(prof_b[, .(key, p, r)], key ~ p, value.var = "r")
setnames(out, c("key", "r10", "r25", "r50", "r75", "r90"))
out[, `:=`(r10 = round(r10, 3), r25 = round(r25, 3), r50 = round(r50, 3),
           r75 = round(r75, 3), r90 = round(r90, 3))]
out[, c("patient", "mp") := tstrsplit(key, " ")]
fwrite(out[, .(patient, mp, r10, r25, r50, r75, r90)],
       file.path(tab, "nc2025_fig5_gradient.csv"))
cat("-> nc2025_fig5_gradient.csv\n")

# ---- 形态判定：修正"展宽"分支的门槛 ----
# 原判据要求 disp > 1.6 且 |z50| < 0.5 才算展宽，结果把 P16-MP7（disp 4.12）
# 与 P16-MP8（disp 4.25）挤进"无变化"。展宽是否成立与中位动不动无逻辑关系，
# 改为"展宽幅度超过位移幅度"即判展宽：|z50| < asym/2。
sh[, call_v2 := fifelse(disp > 1.6 & abs(z50) < asym / 2, "dispersion up",
                fifelse(disp < 0.63 & abs(z50) < 0.5, "dispersion down",
                fifelse(z25 > 0.3 & z50 > 0.3 & z75 > 0.3, "whole shift up",
                fifelse(z25 < -0.3 & z50 < -0.3 & z75 < -0.3, "whole shift down",
                fifelse(z75 > 0.5 & abs(z25) < 0.35 & abs(z50) < 0.35, "right-tail expansion",
                fifelse(z25 > 0.5 & abs(z75) < 0.35, "left truncation", "no change"))))))]
sh[, changed := call_v1 != call_v2]
cat("\n形态判定修正：改名", sum(sh$changed), "例 ->",
    paste(sh[changed == TRUE, paste0(patient, "-", mp,
          " (", call_v1, " -> ", call_v2, ")")], collapse = "; "), "\n")
print(sh[, .N, by = call_v2][order(-N)])
fwrite(sh[, .(patient, mp, z50 = round(z50, 3), asym = round(asym, 3),
              disp = round(disp, 2),
              ex_lead = round(ex_lead, 3), ex_trail = round(ex_trail, 3),
              call_v1, call_v2, changed)],
       file.path(tab, "nc2025_fig5_shape_calls.csv"))
cat("-> nc2025_fig5_shape_calls.csv\n")

cat("\n-- 关键数字 --\n")
cat("实测梯度中位: ", paste(sprintf("Q%s=%.2f", QS * 100, bO$md), collapse = "  "), "\n")
cat("选择反事实:   ", paste(sprintf("Q%s=%.2f", QS * 100, bS$md), collapse = "  "), "\n")
cat("后缘超额最大:", round(max(sh$ex_trail), 3),
    " 前缘超额最大:", round(max(sh$ex_lead), 3), "\n")
cat("上移组合后缘超额最大:", round(max(sh[z50 > 0.5]$ex_trail), 3),
    " 下移组合后缘超额最大:", round(max(sh[z50 < -0.5]$ex_trail), 3), "\n")
cat("逐患者 r10 中位:", paste(sprintf("%s=%.2f", emed[q == "Q10"]$patient,
    emed[q == "Q10"]$ratio), collapse = " "), "\n")
cat("逐患者 r90 中位:", paste(sprintf("%s=%.2f", emed[q == "Q90"]$patient,
    emed[q == "Q90"]$ratio), collapse = " "), "\n")
