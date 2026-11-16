#!/usr/bin/env Rscript
# supp_gc_content_controls.R -- Supplementary figure: GC-content controls (R2-M1).
#   A  per-exon GC, human DCR exons (D) vs non-DCR exons of DCR genes (D*) and of
#      non-DCR genes (ND)                       [data/gc_per_exon_human.tsv]
#   B  per-region GC, mouse orthologs of human DCRs vs random mouse control
#                                               [data/gc_per_region_mouse.tsv]
#      tests for A-B: gc_audit_tests.csv; summaries checked against gc_audit_summary.csv
#   C  dual-coding proportion, DCR vs length-only and GC+length-matched controls
#                                               [gc_matched_fig2b.csv]
#   D  DCR self-shuffle (own codons) null vs observed count   [self_shuffle_control.csv]
#   E  mean phyloP, baseline vs GC+length-matched controls     [phyloP_gc_control.csv]
# Per-exon GC tables are regenerated from the DOOVER beds (bedtools nuc / fasta) by
# figures/make_gc_tables.py.
# Usage: Rscript Code/05_Revision/figures/supp_gc_content_controls.R [indir] [datadir] [outdir]

args <- commandArgs(trailingOnly = TRUE)
indir   <- if (length(args) >= 1) args[1] else "inputs"
datadir <- if (length(args) >= 2) args[2] else "data"
outdir  <- if (length(args) >= 3) args[3] else "figures"
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "dcr_style.R"))
suppressPackageStartupMessages(library(data.table))

summ  <- read_first(file.path(indir, "gc_audit_summary.csv"))
tests <- read_first(file.path(indir, "gc_audit_tests.csv"))
gcm   <- read_first(file.path(indir, "gc_matched_fig2b.csv"))
shuf  <- read_first(file.path(indir, "self_shuffle_control.csv"))
phy   <- read_first(file.path(indir, "phyloP_gc_control.csv"))
hum <- fread(file.path(datadir, "gc_per_exon_human.tsv"), header = FALSE,
             col.names = c("cls", "chrom", "start", "end", "gc", "len"))
mou <- fread(file.path(datadir, "gc_per_region_mouse.tsv"))

# regenerated per-exon values must reproduce the audited summary
chk <- rbind(hum[, .(n = .N, med = round(median(gc), 2)), by = cls],
             mou[, .(n = .N, med = round(median(gc, na.rm = TRUE), 2)), by = cls])
chk <- merge(chk, summ[, c("cls", "n", "gc_median")], by = "cls")
stopifnot(all(chk$n.x == chk$n.y), all(abs(chk$med - chk$gc_median) <= 0.01))

tv <- function(a, b) tests[tests$group_A == a & tests$group_B == b, ]

# ---- A: human ----
cls_map <- c(DCR_exon_human = "D", nonDCR_exon_DCR_gene = "D*", nonDCR_exon_nonDCR_gene = "ND")
hum[, grp := factor(cls_map[cls], levels = c("D", "D*", "ND"))]
nlab <- hum[, .N, by = grp]
brA <- data.frame(x1 = c(1, 1), x2 = c(2, 3), y = c(100, 110),
                  label = c(p_stars(tv("DCR_exon_human", "nonDCR_exon_DCR_gene")$MWU_p),
                            p_stars(tv("DCR_exon_human", "nonDCR_exon_nonDCR_gene")$MWU_p)))
violin_box <- function(d, fills, lines) {
  list(geom_violin(aes(fill = grp, colour = grp), width = 0.85, linewidth = 0.35,
                   scale = "width", trim = TRUE, alpha = 0.55),
       geom_boxplot(aes(fill = grp, colour = grp), width = 0.22, outlier.shape = NA,
                    linewidth = 0.4, fatten = 1.6),
       scale_fill_manual(values = fills, guide = "none"),
       scale_colour_manual(values = lines, guide = "none"))
}
med_col <- function(p) p   # boxplot median line keeps the class line colour (Fig 2C)
gA <- ggplot(hum, aes(grp, gc)) +
  violin_box(hum, EXON_CLASS_COL, c(D = "#E09400", "D*" = "#D9C36A", ND = "grey50")) +
  sig_bracket(brA, tip = 2.5) +
  geom_text(data = nlab, aes(grp, -4, label = paste0("n = ", fmt_n(N))), vjust = 1,
            size = DCR_SIZE$small / PT, family = DCR_FONT, colour = "grey30") +
  scale_y_continuous(limits = c(-11, 120), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = "Human exon class", y = "GC content (%)") +
  theme_dcr()

# ---- B: mouse ----
mou[, grp := factor(ifelse(cls == "mouse_DCR_ortholog", "DCR\northolog", "random\ncontrol"),
                    levels = c("DCR\northolog", "random\ncontrol"))]
mou <- mou[!is.na(gc)]
nlabB <- mou[, .N, by = grp]
brB <- data.frame(x1 = 1, x2 = 2, y = 104,
                  label = p_stars(tv("mouse_DCR_ortholog", "random_ctrl_mouse")$MWU_p))
fB <- c("DCR\northolog" = DCR_COL[["dcr"]], "random\ncontrol" = DCR_COL[["ctrl_light"]])
lB <- c("DCR\northolog" = "#E09400", "random\ncontrol" = "grey50")
gB <- ggplot(mou, aes(grp, gc)) + violin_box(mou, fB, lB) +
  sig_bracket(brB, tip = 2.5) +
  geom_text(data = nlabB, aes(grp, -4, label = paste0("n = ", fmt_n(N))), vjust = 1,
            size = DCR_SIZE$small / PT, family = DCR_FONT, colour = "grey30") +
  scale_y_continuous(limits = c(-11, 120), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = "Mouse region", y = "GC content (%)") +
  theme_dcr()

# ---- C: GC + length matched dual-coding proportion ----
lvC <- c(DCR_observed = "DCR\northolog", random_ctrl_length_only = "control,\nlength-\nmatched",
         random_ctrl_GClength_matched = "control,\nGC + length-\nmatched")
dC <- data.frame(grp = factor(lvC[gcm$comparison], levels = lvC), pct = gcm$dualcoding_pct,
                 lo = suppressWarnings(as.numeric(gcm$ci_low)),
                 hi = suppressWarnings(as.numeric(gcm$ci_high)))
dC$fill <- c(DCR_COL[["dcr"]], DCR_COL[["ctrl_light"]], DCR_COL[["ctrl_dark"]])[match(gcm$comparison, names(lvC))]
gC <- ggplot(dC, aes(grp, pct)) +
  geom_col(aes(fill = fill), width = 0.62) +
  geom_errorbar(data = dC[!is.na(dC$lo), ], aes(ymin = lo, ymax = hi), width = 0.18,
                linewidth = 0.4) +
  geom_text(aes(y = pmax(pct, hi, na.rm = TRUE) + 2, label = sprintf("%.1f%%", pct)),
            vjust = 0, size = DCR_SIZE$annot / PT, family = DCR_FONT) +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "Mouse regions dual-coding (%)") +
  theme_dcr()

# ---- D: self-shuffle null ----
sv <- setNames(shuf$value, shuf$metric)
dD <- data.frame(grp = factor(c("observed", "self-shuffle\nnull"),
                              levels = c("self-shuffle\nnull", "observed")),
                 v = c(sv[["observed_count"]], sv[["null_mean"]]),
                 lo = c(NA, sv[["null_ci_low"]]), hi = c(NA, sv[["null_ci_high"]]),
                 fill = c(DCR_COL[["dcr"]], DCR_COL[["ctrl_dark"]]))
gD <- ggplot(dD, aes(v, grp)) +
  geom_errorbarh(data = dD[!is.na(dD$lo), ], aes(xmin = lo, xmax = hi), height = 0.18,
                 linewidth = 0.45) +
  geom_point(aes(fill = fill), shape = 21, size = 3, stroke = 0.4) +
  annotate("text", x = mean(dD$v), y = 1.5, vjust = 0.5,
           label = sprintf("Z = %.1f", sv[["Z"]]), size = DCR_SIZE$annot / PT,
           family = DCR_FONT) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(1040, 1345), breaks = seq(1100, 1300, 100)) +
  labs(x = "Dual-coding orthologs (count)", y = NULL) +
  theme_dcr()

# ---- E: phyloP with GC control ----
pv <- setNames(phy$value, phy$test)
dE <- data.frame(set = factor(rep(c("baseline", "GC + length-\nmatched"), each = 3),
                              levels = c("baseline", "GC + length-\nmatched")),
                 grp = factor(rep(c("D", "D*", "ND"), 2), levels = c("D", "D*", "ND")),
                 m = c(pv[["baseline_D_mean"]], pv[["baseline_Dg_mean"]], pv[["baseline_Ng_mean"]],
                       pv[["GCmatched_D_mean"]], pv[["GCmatched_Dg_mean"]], pv[["GCmatched_Ng_mean"]]))
dodge <- 0.75
xpos <- function(s, g) as.numeric(s) + (as.numeric(g) - 2) * dodge / 3
brE <- data.frame(x1 = xpos(2, 1), x2 = c(xpos(2, 2), xpos(2, 3)),
                  y = c(2.93, 3.03),
                  label = p_stars(c(pv[["D_vs_matchedDg_MWU_p"]], pv[["D_vs_matchedNg_MWU_p"]])))
gE <- ggplot(dE, aes(set, m, fill = grp)) +
  geom_point(position = position_dodge(dodge), shape = 21, size = 3, stroke = 0.4,
             colour = "black") +
  sig_bracket(brE, tip = 0.025) +
  scale_fill_manual(values = EXON_CLASS_COL, name = NULL) +
  scale_y_continuous(limits = c(2.5, 3.1), breaks = seq(2.5, 3.1, 0.1)) +
  labs(x = NULL, y = "Mean phyloP (100-way)") +
  theme_dcr() + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
                      legend.box.margin = margin(-6, 0, 0, 0))

top <- dcr_grid(list(gA, gB), labels = c("A", "B"), nrow = 1, rel_widths = c(1.35, 1),
                align = "h", axis = "tb")
bot <- dcr_grid(list(gC, gD, gE), labels = c("C", "D", "E"), nrow = 1,
                rel_widths = c(1.05, 1.1, 1), align = "h", axis = "b")
fig <- cowplot::plot_grid(top, bot, ncol = 1, rel_heights = c(1.05, 1))
save_dcr(fig, "Supp_GC_content_controls", 7.2, 6.0, outdir)
write.csv(chk, file.path(outdir, "Supp_GC_content_controls_check.csv"), row.names = FALSE)
