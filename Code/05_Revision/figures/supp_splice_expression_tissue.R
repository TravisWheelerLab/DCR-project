#!/usr/bin/env Rscript
# supp_splice_expression_tissue.R -- Supplementary figure: DCR regions vs controls in
# splice-site conservation, expression class and tissue specificity (R2-M4b/M5/M6).
#   A  splice-site phyloP, DCR vs control, junction x site
#        per-site values: splice_site_phylop.tsv ; tests: splice_site_conservation.csv
#   B  % DCR genes per GTEx expression quintile, 95% CI   [A4_dcr_by_expression.csv]
#   C  tissue-specificity tau, DCR vs non-DCR genes within quintiles
#        per-gene: A4_gene_expression_tau.csv ; tests: A4_tau_by_expression.csv
# Usage: Rscript Code/05_Revision/figures/supp_splice_expression_tissue.R [indir] [outdir]

args <- commandArgs(trailingOnly = TRUE)
indir  <- if (length(args) >= 1) args[1] else "inputs"
outdir <- if (length(args) >= 2) args[2] else "figures"
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "dcr_style.R"))
suppressPackageStartupMessages(library(data.table))

ss   <- fread(file.path(indir, "splice_site_phylop.tsv"))
sst  <- read_first(file.path(indir, "splice_site_conservation.csv"))
q    <- read_first(file.path(indir, "A4_dcr_by_expression.csv"))
tau  <- fread(file.path(indir, "A4_gene_expression_tau.csv"))
taut <- read_first(file.path(indir, "A4_tau_by_expression.csv"))

# ---- A ----
ss <- ss[!is.na(phyloP)]
cnt <- dcast(ss[, .N, by = .(junction, site, cls)], junction + site ~ cls, value.var = "N")
chkA <- merge(cnt, sst, by = c("junction", "site"))
stopifnot(all(chkA$DCR == chkA$n_DCR), all(abs(chkA$control - chkA$n_ctrl) <= 10))
ss[, js := factor(paste0(junction, "\n", site),
                  levels = c("upstream\ndonor", "upstream\nacceptor",
                             "downstream\nacceptor", "downstream\ndonor"))]
ss[, cls := factor(cls, levels = c("DCR", "control"))]
whisk <- ss[, {b <- boxplot.stats(phyloP)$stats; .(lo = b[1], hi = b[5])}, by = .(js, cls)]
ylimA <- c(min(whisk$lo), max(whisk$hi))
sstA <- sst[sst$junction != "ALL", ]
sstA$js <- factor(paste0(sstA$junction, "\n", sstA$site), levels = levels(ss$js))
topA <- whisk[, .(top = max(hi)), by = js]
brA <- merge(data.frame(js = sstA$js, label = p_stars(sstA$p_adj)), topA, by = "js")
xA <- as.numeric(brA$js)
brA <- data.frame(x1 = xA - 0.2, x2 = xA + 0.2, y = brA$top + 0.04 * diff(ylimA),
                  label = brA$label)
gA <- ggplot(ss, aes(js, phyloP, fill = cls, colour = cls)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.4, fatten = 1.6,
               position = position_dodge(0.75), coef = 1.5) +
  sig_bracket(brA, tip = 0.02 * diff(ylimA)) +
  scale_fill_manual(values = c(DCR = DCR_COL[["dcr"]], control = DCR_COL[["ctrl_light"]]),
                    labels = c("DCR exon", "control exon"), name = NULL) +
  scale_colour_manual(values = c(DCR = "#B87900", control = "grey40"),
                      labels = c("DCR exon", "control exon"), name = NULL) +
  coord_cartesian(ylim = c(ylimA[1], ylimA[2] + 0.14 * diff(ylimA))) +
  labs(x = "Junction and splice site", y = "Splice-site phyloP (100-way)") +
  theme_dcr()

# ---- B ----
q$bin <- factor(sub(" \\(.*\\)", "", q$expr_bin), levels = sub(" \\(.*\\)", "", q$expr_bin))
overall <- 100 * sum(q$n_dcr) / sum(q$n)
gB <- ggplot(q, aes(bin, pct_dcr)) +
  scale_x_discrete() +
  geom_hline(yintercept = overall, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.45) +
  geom_point(shape = 21, size = 3, stroke = 0.4, fill = DCR_COL[["dcr"]]) +
  scale_y_continuous(limits = c(0, 10.5), breaks = seq(0, 10, 2.5), expand = c(0, 0)) +
  labs(x = "GTEx expression quintile (low to high)", y = "DCR genes (%)") +
  theme_dcr()

# ---- C ----
tau[, bin := factor(sub(" \\(.*\\)", "", expr_bin), levels = levels(q$bin))]
tau[, grp := factor(ifelse(is_dcr %in% c(TRUE, "True"), "DCR", "non-DCR"),
                    levels = c("DCR", "non-DCR"))]
tau <- tau[!is.na(tau) & !is.na(bin)]
chkC <- dcast(tau[, .N, by = .(bin, grp)], bin ~ grp, value.var = "N")
taut$bin <- factor(sub(" \\(.*\\)", "", taut$expr_bin), levels = levels(q$bin))
chkC <- merge(chkC, taut, by = "bin")
stopifnot(all(chkC$DCR == chkC$n_dcr), all(chkC[["non-DCR"]] == chkC$n_non))
brC <- data.frame(x1 = as.numeric(taut$bin) - 0.2, x2 = as.numeric(taut$bin) + 0.2,
                  y = 1.05, label = p_stars(taut$p))
gC <- ggplot(tau, aes(bin, tau, fill = grp, colour = grp)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.4, fatten = 1.6,
               position = position_dodge(0.75)) +
  sig_bracket(brC, tip = 0.015, text_size = DCR_SIZE$annot, gap = 0.012) +
  scale_fill_manual(values = c(DCR = DCR_COL[["dcr"]], "non-DCR" = DCR_COL[["ctrl_light"]]),
                    labels = c("DCR genes", "non-DCR genes"), name = NULL) +
  scale_colour_manual(values = c(DCR = "#B87900", "non-DCR" = "grey40"),
                      labels = c("DCR genes", "non-DCR genes"), name = NULL) +
  scale_y_continuous(limits = c(0, 1.13), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = "GTEx expression quintile (low to high)",
       y = expression(paste("Tissue specificity ", tau))) +
  theme_dcr()

bot <- dcr_grid(list(gB, gC), labels = c("B", "C"), nrow = 1, rel_widths = c(0.8, 1.2),
                align = "h", axis = "b")
fig <- cowplot::plot_grid(dcr_grid(list(gA), labels = "A"), bot, ncol = 1,
                          rel_heights = c(1, 1))
save_dcr(fig, "Supp_splice_site_expression_tissue_specificity", 7.2, 6.2, outdir)
write.csv(chkA, file.path(outdir, "Supp_splice_expression_tissue_checkA.csv"), row.names = FALSE)
