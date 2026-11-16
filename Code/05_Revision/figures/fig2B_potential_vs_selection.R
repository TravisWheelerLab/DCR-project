#!/usr/bin/env Rscript
# fig2B_potential_vs_selection.R -- new Figure 2 panel B (replaces the ">1 ORF by
# length" panel). Three sub-plots i-iii; every number is read from the CSV inputs.
#
# Inputs (in --indir, default "inputs"):
#   definition_reconciliation.csv  B(i)  DCR exons, permissive sub-interval rule
#   xprime_summary.csv             B(i)  N (non-DCR exons with potential), X
#                                  B(ii) D', N'_matched, N'_raw, CMH p
#   gc_audit_summary.csv           B(i)  denominator of non-DCR coding exons
#   selection_omega_v3.csv (preferred if present) else selection_omega.csv   B(iii)
#
# Slot: the panel-B area of GBE/Figures/Figure2.pdf measured with pdftotext -bbox:
#   x 142.7-706.7 pt, y 0-262 pt  ->  7.83 x 3.64 in. Rendered at 7.8 x 3.6 in.
#
# Usage: Rscript Code/05_Revision/figures/fig2B_potential_vs_selection.R [indir] [outdir]

args <- commandArgs(trailingOnly = TRUE)
indir  <- if (length(args) >= 1) args[1] else "inputs"
outdir <- if (length(args) >= 2) args[2] else "figures"
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "dcr_style.R"))

SLOT_W <- 7.8; SLOT_H <- 3.6
# p for B(iii) is not a column of selection_omega.csv v2; it comes from the
# bootstrap difference test (2,000 replicates over intervals) in
# OLGenie_methods_explained.md. Used ONLY if the CSV has no p column.
P_OMEGA_FALLBACK <- 0.001   # reported as p < 0.001
P_OMEGA_FALLBACK_IS_BOUND <- TRUE
N_BOOT <- 2000              # bootstrap replicates of the v3 difference test

wilson <- function(k, n, z = 1.96) {
  p <- k / n; d <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d)
}
pct_num <- function(s) as.numeric(sub("%", "", sub("~", "", s)))

# ---------------- data ----------------
defr <- read_first(file.path(indir, "definition_reconciliation.csv"))
xs   <- read_first(file.path(indir, "xprime_summary.csv"))
gca  <- read_first(file.path(indir, "gc_audit_summary.csv"))
om   <- read_first(c(file.path(indir, "selection_omega_v3.csv"), "selection_omega_v3.csv",
                     file.path(indir, "selection_omega.csv")))

row_perm <- defr[defr$measurement == "human DCR exon (sub-interval >=21nt)", ]
stopifnot(nrow(row_perm) == 1)
dcr_n   <- row_perm$n
dcr_pct <- pct_num(row_perm$value)

# non-DCR coding exons: N_total (exons carrying a dual-coding-capable interval)
# over all non-DCR coding exons (both non-DCR exon classes of the GC audit).
non_den <- sum(gca$n[gca$cls %in% c("nonDCR_exon_DCR_gene", "nonDCR_exon_nonDCR_gene")])
non_k   <- xs$N_total
non_pct <- 100 * non_k / non_den
X_val   <- xs$X
stopifnot(abs(X_val - non_k / (non_k + xs$D_total)) < 1e-6)

d1 <- data.frame(grp = c("DCR\nexons", "non-DCR\nexons"),
                 pct = c(dcr_pct, non_pct), n = c(dcr_n, non_den),
                 fill = c(DCR_COL[["dcr"]], DCR_COL[["ctrl_dark"]]))
ci <- rbind(wilson(dcr_pct / 100 * dcr_n, dcr_n), wilson(non_k, non_den)) * 100
d1$lo <- ci[, 1]; d1$hi <- ci[, 2]
d1$grp <- factor(d1$grp, levels = d1$grp)

d2 <- data.frame(grp = c("DCR", "latent,\nlength-\nmatched", "latent,\nall"),
                 pct = c(xs$D_prime_pct, xs$N_prime_matched_pct, xs$N_prime_raw_pct),
                 fill = c(DCR_COL[["dcr"]], DCR_COL[["ctrl_dark"]], DCR_COL[["ctrl_light"]]))
d2$grp <- factor(d2$grp, levels = d2$grp)
p_cmh <- xs$p_CMH

num <- function(x) suppressWarnings(as.numeric(x))
if ("scope" %in% names(om)) {
  # v3 schema: one row per scope x comparison; D_* = DCR, N_* = latent null.
  # Use the primary-contig main analysis (alt/random/fix/chrUn contigs excluded).
  r <- om[grepl("primary", om$scope) & om$comparison == "main", ]
  stopifnot(nrow(r) == 1)
  om_scope <- r$scope
  d3 <- data.frame(grp = c("DCR", "latent\nnull"),
                   w  = num(c(r$D_w_dSNdSS, r$N_w_dSNdSS)),
                   lo = num(c(r$D_w_dSNdSS_lo, r$N_w_dSNdSS_lo)),
                   hi = num(c(r$D_w_dSNdSS_hi, r$N_w_dSNdSS_hi)),
                   n  = num(c(r$n_D, r$n_N)),
                   col = c(DCR_COL[["dcr"]], DCR_COL[["ctrl_dark"]]))
  # two-sided bootstrap p; when no replicate crosses 0 the CSV holds 1/B, i.e. p < 1/B
  p_om <- num(r$p_w_dSNdSS)
  p_om_bound <- abs(p_om - 1 / N_BOOT) < 1e-12
} else {
  # v2 schema
  om_scope <- "all contigs (v2)"
  dcr_row  <- om[grepl("^DCR", om$set), ]
  null_row <- om[grepl("null", om$set), ]
  stopifnot(nrow(dcr_row) == 1, nrow(null_row) == 1)
  d3 <- data.frame(grp = c("DCR", "latent\nnull"),
                   w  = num(c(dcr_row$w_alt_dSNdSS, null_row$w_alt_dSNdSS)),
                   lo = num(c(dcr_row$dSNdSS_lo, null_row$dSNdSS_lo)),
                   hi = num(c(dcr_row$dSNdSS_hi, null_row$dSNdSS_hi)),
                   n  = num(c(dcr_row$n, null_row$n)),
                   col = c(DCR_COL[["dcr"]], DCR_COL[["ctrl_dark"]]))
  p_om <- P_OMEGA_FALLBACK; p_om_bound <- P_OMEGA_FALLBACK_IS_BOUND
}
d3$grp <- factor(d3$grp, levels = d3$grp)
om_lab <- if (p_om_bound) "***" else p_stars(p_om)

# ---------------- plots ----------------
bar_w <- 0.62
lab_y <- function(d) pmax(d$pct, if (!is.null(d$hi)) d$hi else d$pct) + 2.5

g1 <- ggplot(d1, aes(grp, pct)) +
  geom_col(aes(fill = fill), width = bar_w) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.4) +
  geom_text(aes(y = lab_y(d1), label = sprintf("%.1f%%", pct)), vjust = 0,
            size = DCR_SIZE$annot / PT, family = DCR_FONT) +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "Exons with a dual-coding-\ncapable interval (%)") +
  theme_dcr()

# effect size, not stars: the point is that conservation barely narrows the field
gap_pts <- d2$pct[1] - d2$pct[2]
br2 <- data.frame(x1 = 1, x2 = 2, y = max(d2$pct[1:2]) + 9,
                  label = sprintf("%.1f pts", gap_pts))
g2 <- ggplot(d2, aes(grp, pct)) +
  geom_col(aes(fill = fill), width = bar_w) +
  geom_text(aes(y = pct + 2, label = sprintf("%.1f%%", pct)), vjust = 0,
            size = DCR_SIZE$annot / PT, family = DCR_FONT) +
  sig_bracket(br2, tip = 2.5, text_size = DCR_SIZE$annot) +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "Still dual-coding in mouse (%)") +
  theme_dcr()

ytop3 <- max(c(d3$w, d3$hi), na.rm = TRUE)
br3 <- data.frame(x1 = 1, x2 = 2, y = max(1, ytop3) + 0.12, label = om_lab)
g3 <- ggplot(d3, aes(grp, w)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = DCR_COL[["neutral_line"]],
             linewidth = 0.5) +
  annotate("text", x = 0.5, y = 1, label = "neutral", hjust = 0, vjust = -0.5,
           size = DCR_SIZE$small / PT, family = DCR_FONT, colour = "grey30") +
  geom_errorbar(data = d3[!is.na(d3$lo), ], aes(ymin = lo, ymax = hi), width = 0.15,
                linewidth = 0.45) +
  geom_point(aes(fill = col), shape = 21, size = 3.2, stroke = 0.4, colour = "black") +
  sig_bracket(br3, tip = 0.035) +
  scale_fill_identity() +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  scale_y_continuous(limits = c(0, 1.3), breaks = seq(0, 1.25, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = expression(paste("Alternate-frame ", omega[alt], " (dSN/dSS)"))) +
  theme_dcr()

row <- sub_grid(list(g1, g2, g3), c("i.", "ii.", "iii."), nrow = 1,
                rel_widths = c(1, 1.22, 1), align = "h", axis = "tb")
# reserve a left gutter for the bold "B", as in the original slot
fig <- cowplot::plot_grid(NULL, row, nrow = 1, rel_widths = c(0.045, 1)) +
  cowplot::draw_label("B", x = 0.002, y = 0.995, hjust = 0, vjust = 1, size = DCR_SIZE$letter,
                      fontface = "bold", fontfamily = DCR_FONT)

save_dcr(fig, "Figure2B_dual_coding_potential_vs_selection", SLOT_W, SLOT_H, outdir)

# numbers actually plotted, for the legend
vals <- data.frame(
  item = c("Bi_DCR_pct", "Bi_DCR_n", "Bi_DCR_lo", "Bi_DCR_hi", "Bi_nonDCR_pct", "Bi_nonDCR_k",
           "Bi_nonDCR_n", "X", "Bii_D_prime", "Bii_N_matched", "Bii_N_raw", "Bii_ratio",
           "Bii_ratio_lo", "Bii_ratio_hi", "Bii_OR_MH", "Bii_p_CMH", "Biii_DCR_w", "Biii_DCR_lo",
           "Biii_DCR_hi", "Biii_DCR_n", "Biii_null_w", "Biii_null_lo", "Biii_null_hi",
           "Biii_null_n", "Biii_p", "Biii_source", "Biii_scope", "Bii_gap_pts"),
  value = c(dcr_pct, dcr_n, d1$lo[1], d1$hi[1], non_pct, non_k, non_den, X_val,
            xs$D_prime_pct, xs$N_prime_matched_pct, xs$N_prime_raw_pct, xs$ratio, xs$ratio_lo,
            xs$ratio_hi, xs$OR_MH, p_cmh, d3$w[1], d3$lo[1], d3$hi[1], d3$n[1], d3$w[2],
            d3$lo[2], d3$hi[2], d3$n[2],
            if (p_om_bound) paste0("<", p_om) else p_om, attr(om, "source"), om_scope, gap_pts))
write.csv(vals, file.path(outdir, "Figure2B_values_plotted.csv"), row.names = FALSE)
