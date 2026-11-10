#!/usr/bin/env Rscript
# ==============================================================================
# 36_figS2_sensitivity.R
#
# 补充图 S2：判据的稳健性。
#
# 主图的所有结论都依赖三个可调参数：
#   (1) 位移阈值 |z50| > 0.5 IQR
#   (2) 每个 患者 x 程序 x 时点对 的最少细胞数 50
#   (3) 时点对的选择（A>B / B>C / A>C 三种，互相不独立）
# 还要排除"结论由个别病人驱动"。
#
# 注意：阈值扫描是**后验重分类**——只改变"是否算实质位移"，
# 不重新判定 call_v2（形状类别依赖原阈值算出的 asym）。文中须注明。
#
# 输出：
#   results/tables/neoletexe_sensitivity.csv
#   results/figures/figS2_sensitivity.png
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

TAB <- "results/tables"; FIG <- "results/figures"

sh <- fread(file.path(TAB, "neoletexe_shape.csv"))
sh[, z50 := as.numeric(z50)]
sh[, `:=`(ex_trail = as.numeric(ex_trail), ex_lead = as.numeric(ex_lead))]
sh[, nmin := pmin(n_from, n_to)]
N0 <- nrow(sh)

## ---------------- A. 阈值扫描 ----------------
thrs <- seq(0.30, 0.90, by = 0.05)
sc <- rbindlist(lapply(thrs, function(t) {
  s <- sh[abs(z50) > t]
  data.table(thr = t, n_shift = nrow(s), rate = nrow(s) / N0,
             trail_max = if (nrow(s[z50 > 0])) max(s[z50 > 0]$ex_trail) else NA_real_,
             trail_n = sum(s$ex_trail > 0.5, na.rm = TRUE))
}))
fwrite(sc, file.path(TAB, "neoletexe_sensitivity.csv"))

## ---------------- B. 细胞数门槛 ----------------
ncs <- c(50, 80, 120, 200, 300, 500)
nb <- rbindlist(lapply(ncs, function(k) {
  s <- sh[nmin >= k]
  data.table(nmin = k, n_pair = nrow(s),
             n_shift = sum(abs(s$z50) > 0.5),
             rate = mean(abs(s$z50) > 0.5),
             med_abs = median(abs(s$z50)),
             trail_max = max(s[z50 > 0.5]$ex_trail, na.rm = TRUE))
}))
cat("\n===== 最小细胞数门槛 =====\n"); print(nb)

## ---------------- C. 时点对 ----------------
pc <- sh[, .(n = .N, shift = sum(abs(z50) > 0.5),
             rate = mean(abs(z50) > 0.5),
             med_z = median(z50), med_abs = median(abs(z50))), by = "pair"]
pc[, pair := factor(pair, levels = c("A>B", "B>C", "A>C"))]
pc <- pc[order(pair)]
cat("\n===== 三个时点对 =====\n"); print(pc)

## ---------------- D. 留一患者 ----------------
pts <- sort(unique(sh$patient))
lo <- rbindlist(lapply(pts, function(p) {
  s <- sh[patient != p]
  data.table(drop = p, rate = mean(abs(s$z50) > 0.5))
}))
base_rate <- mean(abs(sh$z50) > 0.5)
cat(sprintf("\n===== 留一患者：位移率 %.3f，范围 %.3f-%.3f =====\n",
            base_rate, min(lo$rate), max(lo$rate)))
print(lo[order(rate)][c(1:3, .N - 2:0)])

## ---------------- 出图 ----------------
COL <- "#4477AA"; RED <- "#B22222"

pA <- ggplot(sc, aes(thr, 100 * rate)) +
  geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.35, colour = RED) +
  geom_line(linewidth = 0.6, colour = COL) +
  geom_point(size = 1.3, colour = COL) +
  annotate("text", x = 0.52, y = 42, label = "main text\n0.5 IQR",
           hjust = 0, size = 2.0, colour = RED, lineheight = 0.9) +
  scale_x_continuous(breaks = thrs) +
  labs(x = "Shift threshold (IQR)", y = "Combinations called shifted (%)",
       title = "A  Threshold dependence is smooth",
       subtitle = paste0("No cliff at the chosen cut-off: the distribution of ",
                         "|z50|\nis continuous, so conclusions do not hinge on ",
                         "0.5 exactly.")) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.4, colour = "grey35"),
        axis.text.x = element_text(angle = 0, size = 6))

pB <- ggplot(nb, aes(factor(nmin), 100 * rate)) +
  geom_col(width = 0.55, fill = "grey58") +
  geom_text(aes(label = sprintf("%d/%d", n_shift, n_pair)),
            vjust = -0.25, size = 2.0) +
  scale_y_continuous(limits = c(0, 55), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Minimum cells per timepoint", y = "Shift rate (%)",
       title = "B  Not driven by sparse timepoints",
       subtitle = paste0("Raising the floor to 300 cells keeps the\n",
                         "shift rate near one third.")) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.4, colour = "grey35"))

pC <- ggplot(pc, aes(pair, 100 * rate)) +
  geom_col(width = 0.55, fill = COL) +
  geom_text(aes(label = sprintf("%.1f%%\n%d/%d", 100 * rate, shift, n)),
            vjust = -0.2, size = 2.1, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Interval", y = "Shift rate (%)",
       title = "C  Movement concentrates in the first interval",
       subtitle = paste0("A>B = first aromatase inhibitor; B>C = after switch;\n",
                         "A>C = whole course (not independent of the two).")) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.4, colour = "grey35"))

lo[, drop := factor(drop, levels = lo[order(rate)]$drop)]
lo_s <- lo[order(rate)]
brk <- c(head(lo_s$drop, 3), tail(lo_s$drop, 3))
pD <- ggplot(lo, aes(drop, 100 * rate)) +
  geom_hline(yintercept = 100 * base_rate, linetype = "22",
             linewidth = 0.35, colour = RED) +
  geom_point(size = 1.1, colour = "grey30") +
  coord_flip() +
  scale_x_discrete(breaks = brk) +
  labs(x = NULL, y = "Shift rate after dropping that patient (%)",
       title = "D  No single patient carries the result",
       subtitle = sprintf("Leave-one-patient-out over all %d patients; only the\nextremes are labelled. Dashed line = %.1f%%.",
                          length(pts), 100 * base_rate)) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.4, colour = "grey35"),
        axis.text.y = element_text(size = 6.2))

pTop <- ggplot() + theme_void() +
  annotate("text", x = 0, y = 0.98, hjust = 0, vjust = 1, size = 3.3,
           fontface = "bold", colour = "black",
           label = "Figure S2  Robustness of the shape call") +
  annotate("text", x = 0, y = 0.55, hjust = 0, vjust = 1, size = 2.2,
           colour = "grey30", lineheight = 1.05,
           label = paste0("NeoLetExe, 220 combinations. Threshold scans only relabel ",
                          "the same profiles as shifted / not shifted;\n",
                          "the shape category itself is never re-derived.")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

p <- pTop / (pA | pB) / (pC | pD) + plot_layout(heights = c(0.22, 1, 1))

ggsave(file.path(FIG, "figS2_sensitivity.png"), p, width = 7.2, height = 5.6,
       dpi = 300, bg = "white")
cat("\n-> figS2_sensitivity.png\n")
