# fig20：把机制结论接到临床终点（RCB）
# A 基线 MP7 与应答（横断面，n=2 vs 9）
# B 治疗后 MP7 变化与应答（配对，n=1 vs 6）
# C 时点敏感性：MP7 Δ 在全部治疗后时点上的方向

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

tab <- "results/tables"; fig <- "results/figures"
base <- fread(file.path(tab, "nc2025_baseline_mp_rcb.csv"))
call <- fread(file.path(tab, "nc2025_v2_call.csv"))
alltp <- fread(file.path(tab, "nc2025_v3_call_all_tp.csv"))
lab <- fread(file.path(tab, "nc2025_rcb_labels.csv"))

# 统一患者编号为 "P01" 形式
key <- function(x) sprintf("P%02d", as.integer(x))
base[, patient := key(gsub("P", "", patient))]
call[, patient := key(patient)]
alltp[, patient := key(patient)]

resp <- lab[, .(patient, response, rcb)]
call <- merge(call, resp, by = "patient", all.x = TRUE)
alltp <- merge(alltp, resp, by = "patient", all.x = TRUE)

lv <- c("Favorable (RCB 0-I)", "Unfavorable (RCB II-III)")
cols <- c("Favorable (RCB 0-I)" = "#2E86C1", "Unfavorable (RCB II-III)" = "#C0392B")
for (d in list(base, call, alltp)) {
  d[, response := factor(response, levels = c("favorable", "unfavorable"), labels = lv)]
}

# ---------- A 基线 MP7 ----------
pa <- ggplot(base, aes(x = response, y = MP7)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.25,
               colour = "grey35", fill = "grey85") +
  geom_jitter(aes(colour = response), width = 0.06, height = 0, size = 2.6, alpha = 0.9) +
  geom_text_repel(aes(label = patient), size = 2.9, colour = "grey25",
                  min.segment.length = 0, box.padding = 0.28,
                  segment.colour = "grey65", segment.size = 0.25,
                  max.overlaps = 20, seed = 42) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = NULL, y = "Baseline MP7 score (ER-I), patient median",
       title = "A  Baseline ER programme is lower in responders",
       subtitle = "Cross-sectional, 11 patients with baseline cells (\u226550); 2 vs 9") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey35"))

# ---------- B 治疗后变化 ----------
pb_dat <- call[mp == "MP7"]
pb <- ggplot(pb_dat, aes(x = response, y = d_med)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.25,
               colour = "grey35", fill = "grey85") +
  geom_jitter(aes(colour = response), width = 0.06, height = 0, size = 2.6, alpha = 0.9) +
  geom_text_repel(aes(label = patient), size = 2.9, colour = "grey25",
                  min.segment.length = 0, box.padding = 0.28,
                  segment.colour = "grey65", segment.size = 0.25,
                  max.overlaps = 20, seed = 42) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_y_continuous(breaks = seq(-0.4, 1.0, 0.2)) +
  labs(x = NULL, y = "Change in MP7 (last timepoint \u2212 baseline), median",
       title = "B  The only responder falls; non-responders mostly rise",
       subtitle = "Paired, 7 patients; 1 vs 6 \u2014 descriptive only, no test") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey35"))

# ---------- C 时点敏感性 ----------
cc <- alltp[mp == "MP7"]
cc[, tp := factor(tp, levels = c("T2", "T3", "T4", "T5"))]
pc <- ggplot(cc, aes(x = tp, y = d_med, group = patient)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dotted", colour = "grey70", linewidth = 0.3) +
  geom_line(aes(colour = response), linewidth = 0.6, alpha = 0.75) +
  geom_point(aes(colour = response), size = 2.2) +
  geom_text_repel(data = cc[tp %in% c("T3", "T5")],
                  aes(label = patient), size = 2.8, colour = "grey25",
                  nudge_x = 0.18, direction = "y", hjust = 0,
                  segment.size = 0.2, segment.colour = "grey70",
                  max.overlaps = 20, seed = 42) +
  facet_wrap(~ response) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = NULL, y = "Change in MP7 median vs baseline",
       title = "C  Direction depends on which timepoint you pick",
       subtitle = "Dotted lines = \u00b10.2 effect threshold; P16 flips sign between T4 and T5") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey35"))

fig20 <- (pa | pb) / pc +
  plot_layout(heights = c(1, 0.85)) +
  plot_annotation(
    title = "Clinical anchoring of the ER-programme induction",
    subtitle = "RCB read from Fig 1B; verified against Supplementary Data 4 (12 = 3 favorable + 9 unfavorable)"
  )

ggsave(file.path(fig, "fig20_a5_clinical_link.png"), fig20,
       width = 11, height = 8.4, dpi = 300)
cat("saved fig20\n")
