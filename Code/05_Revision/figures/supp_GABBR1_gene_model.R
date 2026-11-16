#!/usr/bin/env Rscript
# supp_GABBR1_gene_model.R -- Supplementary figure: GABBR1 DCR 1812 in the style of
# Figure 2E. A: SwissProt isoforms Q9UBS5-1..-5 (DOOVER BED12 models) with zoom on
# the 26-nt DCR; B: the DCR sequence (transcript orientation) and its translation in
# both frames. Inputs from figures/prep_gene_models.py:
#   data/gm_gabbr1_exons.csv, gm_gabbr1_dcr.csv, gm_gabbr1_iso.csv,
#   gm_gabbr1_dcr_seq.csv, gm_gabbr1_frames.csv
# Usage: Rscript Code/05_Revision/figures/supp_GABBR1_gene_model.R [datadir] [outdir]

args <- commandArgs(trailingOnly = TRUE)
datadir <- if (length(args) >= 1) args[1] else "data"
outdir  <- if (length(args) >= 2) args[2] else "figures"
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "dcr_style.R"))
source(file.path(script_dir, "gene_model.R"))

ex  <- fread(file.path(datadir, "gm_gabbr1_exons.csv"))
dcr <- fread(file.path(datadir, "gm_gabbr1_dcr.csv"))
seq <- fread(file.path(datadir, "gm_gabbr1_dcr_seq.csv"))
frs <- fread(file.path(datadir, "gm_gabbr1_frames.csv"))
stopifnot(nrow(dcr) == 1, dcr$end - dcr$start == nrow(seq[track == "nt"]))

row_order <- sort(unique(ex$accession))
dcr_exon_span <- ex[dcr_exon == TRUE, .(s = min(start), e = max(end))]
zoom <- c(dcr_exon_span$s - 60, dcr_exon_span$e + 60)
gA <- gene_model_plot(ex, dcr, zoom, row_order, n_ticks = 4)

# ---- B: sequence + both translations (5'->3' of the transcript; minus strand) ----
nt <- seq[track == "nt"]; aa <- seq[track == "aa"]
L <- nrow(nt)
cw <- 1                                  # one nucleotide = 1 unit
role_y <- c(canonical = -1.15, alternative = -2.05)
aa[, `:=`(y = role_y[role], fill = ifelse(role == "canonical", DCR_COL[["frame_can"]],
                                          DCR_COL[["frame_alt"]]))]
frs[, lab := sprintf("%s  %s", vapply(strsplit(accessions, ","), collapse_iso, ""),
                     frame_lab(offset, frame_role))]
frs[, y := role_y[frame_role]]
xoff <- 9.5                              # room for row labels
gB <- ggplot() +
  annotate("rect", xmin = xoff - 0.5, xmax = xoff + L - 0.5, ymin = -0.45, ymax = 0.45,
           fill = DCR_COL[["dcr_highlight"]]) +
  geom_text(data = nt, aes(x = xoff + i, y = 0, label = char), size = 7 / PT,
            family = "Courier") +
  geom_rect(data = aa, aes(xmin = xoff + i - 0.45, xmax = xoff + i + 2.45, ymin = y - 0.32,
                           ymax = y + 0.32, fill = fill), colour = "white", linewidth = 0.3) +
  geom_text(data = aa, aes(x = xoff + i + 1, y = y, label = char), size = 7 / PT,
            family = DCR_FONT, colour = "white", fontface = "bold") +
  geom_text(data = frs, aes(x = xoff - 1, y = y, label = lab), hjust = 1,
            size = DCR_SIZE$small / PT, family = DCR_FONT) +
  annotate("text", x = xoff - 1, y = 0, label = "DCR 1812 (5'-3')", hjust = 1,
           size = DCR_SIZE$small / PT, family = DCR_FONT) +
  annotate("text", x = xoff + c(0, L - 1), y = 0.75,
           label = fmt_coord(c(nt$gpos[1], nt$gpos[L])), size = DCR_SIZE$small / PT,
           family = DCR_FONT, vjust = 0, hjust = c(0.2, 0.8)) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-9.5, xoff + L + 0.5), ylim = c(-2.6, 1.3), expand = FALSE) +
  theme_void(base_family = DCR_FONT) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(14, 4, 4, 4))

fig <- dcr_grid(list(gA, gB), labels = c("A", "B"), ncol = 1, rel_heights = c(3.2, 1),
                pad_top = 14, pad_left = 4)
save_dcr(fig, "Supp_GABBR1_DCR_gene_model", 7.2, 4.4, outdir)
