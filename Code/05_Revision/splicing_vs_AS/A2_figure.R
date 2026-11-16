#!/usr/bin/env Rscript
# A2_figure.R -- Supplementary Figure S10: DCRs vs ordinary alternative splicing (length, frame
# preservation, NMD, conservation, inclusion-matched conservation).
# Inputs  (A2 output dir, default $WORK/splicing_vs_AS): A2_events_annotated.csv, A2_table_NMD_tests.csv,
#         A2_table_conservation.csv, A2_table_frame_preservation.csv,
#         A2_table_conservation_inclusion_matched.csv   (from A2_stats.py)
# Outputs (figdir): A2_FigureS_DCR_vs_regularAS.pdf/.png ; A2_figure_panelE_data.csv in the A2 dir
# Usage:  Rscript A2_figure.R [a2dir] [figdir]
args <- commandArgs(trailingOnly = TRUE)
a2dir  <- if (length(args) >= 1) args[1] else file.path(Sys.getenv("WORK", "work"), "splicing_vs_AS")
figdir <- if (length(args) >= 2) args[2] else file.path(Sys.getenv("WORK", "work"), "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
# A2 supplementary figure: DCRs vs regular alternative splicing (house style: classic theme, no panel titles,
# bold capital panel letters, legends at bottom, DCR orange vs greys; brackets with stars)
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(cowplot); library(scales)})
E <- read.csv(file.path(a2dir, "A2_events_annotated.csv"), stringsAsFactors = FALSE)
nmd_t  <- read.csv(file.path(a2dir, "A2_table_NMD_tests.csv"))
cons_t <- read.csv(file.path(a2dir, "A2_table_conservation.csv"))
fp_t   <- read.csv(file.path(a2dir, "A2_table_frame_preservation.csv"))
inc_t  <- read.csv(file.path(a2dir, "A2_table_conservation_inclusion_matched.csv"))

lev  <- c("DCR","regASfp_DCRgenes","regASfp_nonDCRgenes","fsAS_DCRgenes","fsAS_nonDCRgenes","const_DCRgenes","const_nonDCRgenes")
labs <- c("D","fp\nD*","fp\nND","fs\nD*","fs\nND","c\nD*","c\nND")
grp  <- c("DCR","DCR gene","non-DCR gene","DCR gene","non-DCR gene","DCR gene","non-DCR gene")
names(labs) <- lev; names(grp) <- lev
E$arm2 <- factor(E$arm2, levels = lev)
E$grp  <- factor(grp[as.character(E$arm2)], levels = c("DCR","DCR gene","non-DCR gene"))
fills <- c("DCR" = "#FFA800", "DCR gene" = "#FFF0BD", "non-DCR gene" = "gray80")
lines <- c("DCR" = "#FFA800", "DCR gene" = "#E3B95A", "non-DCR gene" = "gray60")
boxfills <- c("DCR" = "#FFD98C", "DCR gene" = "#FFF0BD", "non-DCR gene" = "gray80")  # lighter DCR fill so the median line is visible (as in Fig 2C)
stars <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", "ns")))
th <- theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(legend.position = "bottom", axis.text = element_text(colour = "black"), legend.title = element_blank(),
        legend.key.size = unit(3, "mm"))
xsc <- scale_x_discrete(labels = labs)
bracket <- function(xend, y, lab, x0 = 1, tick = 0.02, size = 2.6) {
  list(annotate("segment", x = x0, xend = xend, y = y, yend = y, linewidth = 0.3),
       annotate("segment", x = c(x0, xend), xend = c(x0, xend), y = y, yend = y - tick, linewidth = 0.3),
       annotate("text", x = (x0 + xend) / 2, y = y, label = lab, vjust = -0.2, size = size))
}
AS5 <- lev[1:5]

# A - length
dA <- E %>% filter(arm2 %in% AS5)
pA <- ggplot(dA, aes(arm2, length, fill = grp, colour = grp)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, linewidth = 0.4, width = 0.65) +
  scale_y_log10(breaks = c(10, 30, 100, 300, 1000), labels = comma) + xsc +
  scale_fill_manual(values = boxfills, drop = FALSE) + scale_colour_manual(values = lines, drop = FALSE) +
  labs(x = NULL, y = "Segment length (nt)") + th + theme(legend.position = "none")

# B - multiple of 3 (all internal regular AS, before the frame-preservation split)
dB <- data.frame(arm = factor(c("DCR","regAS_DCRgenes","regAS_nonDCRgenes"), levels = c("DCR","regAS_DCRgenes","regAS_nonDCRgenes")),
                 grp = factor(c("DCR","DCR gene","non-DCR gene"), levels = names(fills)),
                 frac = c(fp_t$DCR_frac_mult3[1], fp_t$ref_frac_mult3))
pB <- ggplot(dB, aes(arm, 100 * frac, fill = grp, colour = grp)) + geom_col(width = 0.65, linewidth = 0.4) +
  scale_x_discrete(labels = c("D", "AS\nD*", "AS\nND")) +
  scale_fill_manual(values = fills) + scale_colour_manual(values = lines) +
  scale_y_continuous(limits = c(0, 85), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Length multiple of 3 (%)") + th + theme(legend.position = "none") +
  bracket(2, 66, stars(fp_t$p[1]), tick = 1.5) + bracket(3, 76, stars(fp_t$p[2]), tick = 1.5)

# C - NMD (pipeline 50-nt rule, event level)
dC <- E %>% filter(arm2 %in% AS5, !(single_model_exon_gene %in% c("True", TRUE)), nmd_evaluable %in% c("True", TRUE)) %>%
  group_by(arm2, grp) %>% summarise(n = n(), k = sum(nmd_arm_specific %in% c("True", TRUE)), .groups = "drop") %>%
  mutate(p = k / n, lo = mapply(function(k, n) binom.test(k, n)$conf.int[1], k, n),
         hi = mapply(function(k, n) binom.test(k, n)$conf.int[2], k, n))
ref_order <- c("regASfp_DCRgenes","regASfp_nonDCRgenes","fsAS_DCRgenes","fsAS_nonDCRgenes")
pC <- ggplot(dC, aes(arm2, 100 * p, fill = grp, colour = grp)) + geom_col(width = 0.65, linewidth = 0.4) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = 0.2, colour = "black", linewidth = 0.3) + xsc +
  scale_fill_manual(values = fills) + scale_colour_manual(values = lines) +
  scale_y_continuous(limits = c(0, 62), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "NMD-susceptible events (%)") + th + theme(legend.position = "none")
for (i in seq_along(ref_order)) {
  pv <- nmd_t$p_adj_GEE[nmd_t$reference == ref_order[i]]
  pC <- pC + bracket(i + 1, 38 + 5.5 * (i - 1), stars(pv), tick = 1)
}

# D - phyloP (length/GC-matched Wilcoxon for brackets)
dD <- E %>% filter(phyloP_bp_frac >= 0.9)
pD <- ggplot(dD, aes(arm2, phyloP_mean, fill = grp, colour = grp)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4, width = 0.65) + xsc +
  scale_fill_manual(values = boxfills, drop = FALSE) + scale_colour_manual(values = lines, drop = FALSE) +
  coord_cartesian(ylim = c(-2.2, 13.2)) + scale_y_continuous(breaks = seq(-2, 8, 2)) +
  labs(x = NULL, y = "Mean PhyloP (100-way)") + th + theme(legend.position = "none")
for (i in 2:7) {
  pv <- cons_t$p_matched_wilcoxon[cons_t$reference == lev[i]]
  pD <- pD + bracket(i, 7.4 + 0.95 * (i - 2), stars(pv), tick = 0.25)
}

# E - phyloP by inclusion level (median, bootstrap 95% CI)
set.seed(1)
bootmed <- function(x) { b <- replicate(1000, median(sample(x, replace = TRUE))); c(median(x), quantile(b, c(.025, .975))) }
dE <- dD %>% filter(arm2 %in% lev) %>%
  mutate(stratum = cut(incl_frac, c(0, 0.5, 0.75, 0.999, 1.01), right = FALSE, labels = c("<0.50", "0.50-0.75", "0.75-<1", "1")),
         stratum = ifelse(grepl("^const", arm2), "1", as.character(stratum)),
         cls = ifelse(arm2 == "DCR", "DCR", ifelse(grepl("^const", arm2), "constitutive", "regular AS")),
         grp2 = as.character(grp)) %>%
  filter(!(arm2 %in% c("fsAS_DCRgenes", "fsAS_nonDCRgenes"))) %>%
  group_by(stratum, arm2, grp2, cls) %>% summarise(n = n(), m = bootmed(phyloP_mean)[1], lo = bootmed(phyloP_mean)[2],
                                                    hi = bootmed(phyloP_mean)[3], .groups = "drop") %>%
  mutate(stratum = factor(stratum, levels = c("<0.50", "0.50-0.75", "0.75-<1", "1")),
         grp2 = factor(grp2, levels = names(fills)))
write.csv(dE, file.path(a2dir, "A2_figure_panelE_data.csv"), row.names = FALSE)
pd <- position_dodge(width = 0.55)
pE <- ggplot(dE, aes(stratum, m, colour = grp2, fill = grp2, shape = cls, group = arm2)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, position = pd, linewidth = 0.4) +
  geom_point(size = 2, position = pd, stroke = 0.6) +
  scale_colour_manual(values = c("DCR" = "#FFA800", "DCR gene" = "#E3B95A", "non-DCR gene" = "gray55"),
                      labels = c("DCR (D)", "DCR gene (D*)", "non-DCR gene (ND)")) +
  scale_fill_manual(values = fills, guide = "none") +
  scale_shape_manual(values = c("DCR" = 21, "regular AS" = 22, "constitutive" = 24),
                     breaks = c("DCR", "regular AS", "constitutive"),
                     labels = c("DCR", "frame-preserving AS (fp)", "constitutive (c)")) +
  labs(x = "Fraction of the gene's isoforms coding the segment", y = "Median PhyloP (100-way)") + th +
  guides(colour = guide_legend(nrow = 1, override.aes = list(shape = 21, fill = fills)),
         shape = guide_legend(nrow = 1, override.aes = list(colour = "black", fill = "white"))) +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0))

top <- plot_grid(pA, pB, pC, nrow = 1, rel_widths = c(1.15, 0.75, 1.15), labels = c("A", "B", "C"),
                 label_size = 12, label_fontface = "bold", label_fontfamily = "Helvetica")
bot <- plot_grid(pD, pE, nrow = 1, rel_widths = c(1.2, 1.25), labels = c("D", "E"),
                 label_size = 12, label_fontface = "bold", label_fontfamily = "Helvetica", align = "h", axis = "tb")
fig <- plot_grid(top, bot, ncol = 1, rel_heights = c(1, 1.22))
ggsave(file.path(figdir, "A2_FigureS_DCR_vs_regularAS.pdf"), fig, width = 180, height = 150, units = "mm", device = "pdf", useDingbats = FALSE, family = "Helvetica")
ggsave(file.path(figdir, "A2_FigureS_DCR_vs_regularAS.png"), fig, width = 180, height = 150, units = "mm", dpi = 300)
