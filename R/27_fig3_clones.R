# fig29：Fig 3 排除预存耐药克隆——六个 panel
#   一句话：每个留存的克隆都各自上移，位移越出基线克隆之间的差异，
#           所以这次变化不是"换了细胞群体"，而是"同一群细胞变了状态"。
#
#   A  六例主导程序的留存克隆斜率图（基线 -> 末次随访），20/20 上移
#   B  分离度 sep_z：随访最差克隆 - 基线最好克隆（IQR 单位）
#   C  组成反事实：观测位移 vs 仅改变克隆占比所能产生的位移
#   D  克隆间 ICC：克隆身份解释了多少随访期的程序方差
#   E  克隆基线位置 vs 位移：排除"均值回归"这一替代解释
#   F  消失克隆的基线偏移：P09 例外的来源
#
# 数据：nc2025_v4_clone_{layer,summary,verdict}.csv + composition_counterfactual.csv

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

tab <- "results/tables"; fig <- "results/figures"

theme_set(theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 10.0),
                  plot.subtitle = element_text(size = 7.3, colour = "grey30",
                                               lineheight = 1.06),
                  panel.grid.minor = element_blank(),
                  plot.title.position = "plot"))

COL_B <- "#7F8C8D"; COL_T <- "#C0392B"; COL_C <- "#2C3E50"; COL_U <- "#2E5FA3"
COL_W <- "#B7950B"   # P09 阴性对照

MPNAME <- c(MP2 = "MP2 Prolif", MP8 = "MP8 Lum-sec", MP7 = "MP7 ER-I",
            MP3 = "MP3 Stromal", MP11 = "MP11 EMT-III", MP6 = "MP6 TGFb",
            MP9 = "MP9 Neural", MP4 = "MP4 IFN", MP5 = "MP5 IEG")
PATORD <- c("P01", "P09", "P13", "P16", "P18", "P23")   # 与 fig28 一致，去 P12（无亚克隆）

key <- function(x) {               # 幂等：整数 1 / 字符 "01" / 已带前缀 "P01" 都能处理
  y <- suppressWarnings(as.integer(x))
  ifelse(is.na(y), as.character(x), sprintf("P%02d", y))
}

ly <- fread(file.path(tab, "nc2025_v4_clone_layer.csv"))
sm <- fread(file.path(tab, "nc2025_v4_clone_summary.csv"))
vr <- fread(file.path(tab, "nc2025_v4_clone_verdict.csv"))
cf <- fread(file.path(tab, "nc2025_v4_composition_counterfactual.csv"))
shp <- fread(file.path(tab, "nc2025_v3_shape.csv"))
shp[, patient := key(patient)]
IQRB <- setNames(as.numeric(shp$iqr_base), paste(shp$patient, shp$mp))
for (d in list(ly, sm, cf)) d[, patient := key(patient)]
vr[, patient := key(patient)]

dom <- sm[is_dominant == "yes"][order(match(patient, key(PATORD)))]
SEP <- setNames(dom$sep_z, dom$patient)
ICC <- setNames(dom$icc_post, dom$patient)
NSEP <- setNames(grepl("not_separated", vr$verdict), vr$patient)

# 只保留主导程序下的克隆层（20 个留存克隆的全部证据都在这里）
dl <- merge(ly[, .(patient, mp, clone, status, n_base, n_post,
                   mb = as.numeric(med_base), mpst = as.numeric(med_post))],
            dom[, .(patient, mp_dom = mp, sep_z, n_retained)],
            by = "patient")[mp == mp_dom]
dl[, patient := factor(patient, levels = PATORD)]
dl[, is_p09 := patient == "P09"]

# ================= A：留存克隆的斜率 =================
dA <- dl[status == "retained" & !is.na(mpst)]
dA[, xb := 1L][, xp := 2L]
dA[, clone_lab := sprintf("%s-c%02d", patient, as.integer(clone))]
dA[, facet := factor(sprintf("%s  %s  (sep %+.2f)", patient, MPNAME[mp], sep_z),
                     levels = {
                       tt <- unique(dA[, .(patient, mp, sep_z)])
                       tt[order(match(patient, PATORD))][, sprintf("%s  %s  (sep %+.2f)",
                                                                  patient, MPNAME[mp], sep_z)]
                     })]

pA <- ggplot(dA) +
  geom_segment(aes(x = xb, xend = xp, y = mb, yend = mpst, group = clone_lab),
               colour = COL_T, linewidth = 0.38, alpha = 0.72,
               arrow = arrow(length = unit(0.018, "in"), type = "open")) +
  geom_point(aes(x = xb, y = mb), colour = COL_B, size = 1.2) +
  geom_point(aes(x = xp, y = mpst), colour = COL_T, size = 1.5) +
  facet_wrap(~ facet, nrow = 2, scales = "free_y") +
  scale_x_continuous(breaks = 1:2, labels = c("baseline", "last on-tx"),
                     limits = c(0.70, 2.30)) +
  labs(title = "A  Every retained clone moves up, on its own",
       subtitle = paste0("One line per clone of that patient's dominant programme ",
                         "(n = 20 clones, six patients).\n20 of 20 move in the same direction; ",
                         "no clone moves against the shift. Facet label gives that patient's ",
                         "separation."),
       x = NULL, y = "Programme score (clone median)") +
  theme(strip.text = element_text(size = 6.8),
        axis.text.x = element_text(size = 6.6),
        axis.text.y = element_text(size = 6.6),
        panel.spacing = unit(0.5, "cm"))

# ================= B：分离度 =================
dB <- copy(dom)
dB[, patient := factor(patient, levels = rev(PATORD))]
dB[, lab := sprintf("%s  %s", patient, MPNAME[mp])]
dB[, lab := factor(lab, levels = rev(sprintf("%s  %s", PATORD, MPNAME[mp[match(PATORD, patient)]])))]
dB[, neg := sep_z < 0]

pB <- ggplot(dB, aes(x = sep_z, y = lab)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.45) +
  geom_col(aes(fill = neg), width = 0.62) +
  geom_text(aes(label = sprintf("%+.2f", sep_z), hjust = ifelse(sep_z < 0, 1.12, -0.12)),
            size = 2.4, colour = "grey15") +
  scale_fill_manual(values = c(`FALSE` = COL_T, `TRUE` = COL_W), guide = "none") +
  scale_x_continuous(limits = c(-0.95, 6.85), breaks = c(0, 2, 4, 6)) +
  labs(title = "B  After treatment, the worst clone beats the best one",
       subtitle = paste0("Separation = lowest clone median after treatment\nminus highest ",
                         "before it, in baseline IQR.\nPositive in five of six."),
       x = "Separation (baseline IQR)", y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 7.0))

# ================= C：组成反事实 =================
dC <- cf[!is.na(obs_delta)][order(match(patient, PATORD))]
dC[, patient := factor(patient, levels = rev(PATORD))]
dC[, obs := as.numeric(obs_delta)][, cfv := as.numeric(cf_delta)]
dC[, share := as.numeric(cf_share)]
dC[, lab := sprintf("%s  %s", patient, MPNAME[mp])]
dC[, lab := factor(lab, levels = rev(dC$lab))]
dC[, is_p09 := patient == "P09"]

pC <- ggplot(dC) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_segment(aes(x = cfv, xend = obs, y = lab, yend = lab), colour = "grey72",
               linewidth = 0.9) +
  geom_point(aes(x = cfv, y = lab), colour = COL_B, size = 2.3) +
  geom_point(aes(x = obs, y = lab, colour = is_p09), size = 2.3) +
  geom_text(aes(x = pmax(obs, cfv, 0) + 0.05, y = lab,
                label = sprintf("%d%%", round(share * 100))),
            hjust = 0, size = 2.2, colour = "grey25") +
  scale_colour_manual(values = c(`FALSE` = COL_T, `TRUE` = COL_W), guide = "none") +
  scale_x_continuous(limits = c(-0.16, 2.55)) +
  labs(title = "C  Composition explains almost none of it",
       subtitle = paste0("Red = observed shift in median; grey = what\nreshuffling clone ",
                         "proportions alone would produce.\nPercentage = share attributable ",
                         "to composition."),
       x = "Shift in median programme score", y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 7.0))

# ================= D：克隆间 ICC =================
dD <- copy(dom)
dD[, patient := factor(patient, levels = rev(PATORD))]
dD[, lab := factor(sprintf("%s  %s", patient, MPNAME[mp]),
                   levels = rev(sprintf("%s  %s", PATORD, MPNAME[mp[match(PATORD, patient)]])))]
icc_med <- median(dom$icc_post)

pD <- ggplot(dD, aes(x = icc_post, y = lab)) +
  geom_vline(xintercept = icc_med, colour = "grey45", linetype = "22", linewidth = 0.35) +
  geom_col(fill = COL_C, width = 0.62) +
  geom_text(aes(label = sprintf("%.3f", icc_post), hjust = -0.15), size = 2.4,
            colour = "grey15") +
  scale_x_continuous(limits = c(0, 0.36), breaks = c(0, 0.1, 0.2, 0.3)) +
  labs(title = "D  Clone identity explains a minority of the variance",
       subtitle = sprintf(paste0("Share of on-treatment variance attributable\nto clone ",
                                 "identity. Median %.2f, range %.3f-%.3f;\nthe largest still ",
                                 "leaves 71%% unexplained."),
                          icc_med, min(dom$icc_post), max(dom$icc_post)),
       x = "Intraclass correlation (clone)", y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 7.0))

# ================= E：基线位置 vs 位移（排除均值回归）=================
dE <- dl[status == "retained" & !is.na(mpst)]
dE[, iqr_b := IQRB[paste(patient, mp)]]
dE[, xb_z := (mb - mean(mb)) / iqr_b, by = patient]
dE[, dy_z := (mpst - mb) / iqr_b]
r_E <- cor(dE$xb_z, dE$dy_z)
n_up <- sum(dE$dy_z > 0)

pE <- ggplot(dE, aes(x = xb_z, y = dy_z)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = is_p09), size = 2.2, alpha = 0.9) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey30",
              linewidth = 0.45, linetype = "21") +
  scale_colour_manual(values = c(`FALSE` = COL_T, `TRUE` = COL_W), guide = "none") +
  labs(title = "E  Not regression to the mean",
       subtitle = sprintf(paste0("Each clone by where it started (baseline,\nin baseline ",
                                 "IQR from the patient's own mean) and how far\nit moved. If ",
                                 "clones merely regressed, high starters\nwould fall. r = %+.2f; ",
                                 "%d of %d moved up."), r_E, n_up, nrow(dE)),
       x = "Baseline position of the clone (IQR)", y = "Shift (IQR)") +
  theme(axis.text = element_text(size = 7.0))

# ================= F：消失克隆的基线偏移 =================
dF <- copy(dom)
dF[, patient := factor(patient, levels = rev(PATORD))]
dF[, lab := factor(sprintf("%s  (%d lost of %d)", patient, n_dropped, n_clone_base),
                   levels = rev(sprintf("%s  (%d lost of %d)", PATORD,
                                        n_dropped[match(PATORD, patient)],
                                        n_clone_base[match(PATORD, patient)])))]
dF[, neg := drop_bias < -0.15]

pF <- ggplot(dF, aes(x = drop_bias, y = lab)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.4) +
  geom_col(aes(fill = neg), width = 0.62) +
  geom_text(aes(label = sprintf("%+.2f", drop_bias),
                hjust = ifelse(drop_bias < 0, 1.12, -0.12)),
            size = 2.4, colour = "grey15") +
  scale_fill_manual(values = c(`FALSE` = COL_B, `TRUE` = COL_W), guide = "none") +
  scale_x_continuous(limits = c(-0.42, 0.20)) +
  labs(title = "F  Only in P09 were the lost clones the low ones",
       subtitle = paste0("Baseline score of clones that disappear minus those\nthat persist ",
                         "(IQR). Near zero in five cases, so loss looks\nsampling-driven. P09 is ",
                         "the exception that keeps the\ncriterion falsifiable."),
       x = "Baseline offset of lost clones (IQR)", y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 7.0))

# ================= 拼图 =================
## 图级标题与图内解释段于 2026-11-10 移除（改由图注承担）。保留面板字母与
## 每个面板的一行标题，与其后五张主图一致。
design <- "
AAB
AAC
DEF
"
f3 <- wrap_plots(A = pA, B = pB, C = pC, D = pD, E = pE, F = pF,
                 design = design) +
  plot_layout(heights = c(1.05, 1.05, 0.92))

ggsave(file.path(fig, "fig29_fig3_clones.png"), f3,
       width = 11.4, height = 9.4, dpi = 300, bg = "white")
cat("-> fig29_fig3_clones.png\n")

# ================= 核对输出 =================
cat("\n-- A：留存克隆同向 --\n")
print(dE[, .(n = .N, up = sum(dy_z > 0), down = sum(dy_z <= 0)), by = patient][
  order(match(patient, PATORD))])
cat(sprintf("合计 n=%d，上移 %d\n", nrow(dE), n_up))

cat("\n-- B 分离度 --\n"); print(dom[, .(patient, mp, sep_z = round(sep_z, 2))])
cat("\n-- C 组成解释率 --\n")
print(cf[!is.na(cf_share), .(patient, mp, obs = round(obs_delta, 3),
                             cf = round(cf_delta, 3), share = cf_share)])
cat(sprintf("\n-- D ICC：中位 %.3f，范围 %.3f-%.3f（全 54 行中位 %.3f）--\n",
            icc_med, min(dom$icc_post), max(dom$icc_post), median(sm$icc_post)))
cat(sprintf("-- E r(基线位置, 位移) = %+.3f --\n", r_E))
cat("\n-- F drop_bias --\n")
print(dom[, .(patient, n_dropped, n_clone_base, drop_bias = round(drop_bias, 3))])
