#!/usr/bin/env Rscript
# supp_GNAS_gene_model.R -- Supplementary figure replacing "data not shown" (R1-m1):
# GNAS DCR 201, the only DCR > 500 nt conserved as dual-coding in mouse.
#   A  human GNAS SwissProt isoforms (DOOVER BED12) with zoom on DCR 201, Figure 2E style
#   B  mouse ortholog of DCR 201 (liftOver to mm39): the three forward frames, stop codons,
#      and the two full-length ORFs; frame colour = human frame it is homologous to
# Inputs: data/gm_gnas_exons.csv, gm_gnas_dcr.csv, gnas_mouse_frames.csv,
#         gnas_mouse_summary.csv (figures/prep_gene_models.py); inputs/gnas_dcr_gt500.csv
# Usage: Rscript Code/05_Revision/figures/supp_GNAS_gene_model.R [datadir] [indir] [outdir]

args <- commandArgs(trailingOnly = TRUE)
datadir <- if (length(args) >= 1) args[1] else "data"
indir   <- if (length(args) >= 2) args[2] else "inputs"
outdir  <- if (length(args) >= 3) args[3] else "figures"
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "dcr_style.R"))
source(file.path(script_dir, "gene_model.R"))

ex   <- fread(file.path(datadir, "gm_gnas_exons.csv"))
dcr  <- fread(file.path(datadir, "gm_gnas_dcr.csv"))
mfr  <- fread(file.path(datadir, "gnas_mouse_frames.csv"))
msum <- fread(file.path(datadir, "gnas_mouse_summary.csv"))
gt5  <- read_first(file.path(indir, "gnas_dcr_gt500.csv"))
g <- gt5[gt5$gene == "gnas", ]
stopifnot(g$dceID == dcr$id, g$start == dcr$start, g$end == dcr$end,
          g$mouse_full_length_ORFs == sum(!is.na(msum$full_length_orf_aa)))

# display order: NESP55, XLas, ALEX, Gs-alpha (5'->3' order of first exons)
first <- ex[, .(s = min(iso_start)), by = accession][order(s, accession)]
row_order <- first$accession
zoom <- c(dcr$start - 260, dcr$end + 420)
gA <- gene_model_plot(ex, dcr, zoom, row_order, n_ticks = 5)

# ---- B: mouse ortholog frames ----
L <- msum$mouse_len[1]
role_of <- ifelse(msum$kmer_share_canonical > msum$kmer_share_alternative &
                    msum$kmer_share_canonical > 0.1, "canonical",
           ifelse(msum$kmer_share_alternative > 0.1, "alternative", "none"))
msum[, role := role_of]
msum[, fill := ifelse(!is.na(full_length_orf_aa),
                      ifelse(role == "canonical", DCR_COL[["frame_can"]],
                             ifelse(role == "alternative", DCR_COL[["frame_alt"]], "grey60")),
                      DCR_COL[["ctrl_light"]])]
msum[, y := -frame * 0.5]
msum[, lab := ifelse(!is.na(full_length_orf_aa),
                     sprintf("%d aa, no stop", full_length_orf_aa),
                     sprintf("%d stop codons", n_stops))]
mfr <- merge(mfr, msum[, .(frame, y)], by = "frame")
x0 <- 0.14; x1 <- 0.995
mx <- function(p) x0 + p / L * (x1 - x0)
ticks <- c(0, 500, 1000, 1500, L)
gB <- ggplot() +
  annotate("segment", x = mx(0), xend = mx(L), y = 0, yend = 0, linewidth = 0.3) +
  annotate("segment", x = mx(ticks), xend = mx(ticks), y = 0, yend = 0.06, linewidth = 0.3) +
  annotate("text", x = mx(ticks), y = 0.1, vjust = 0, size = DCR_SIZE$small / PT,
           family = DCR_FONT, hjust = c(0, 0.5, 0.5, 0.5, 1),
           label = c(paste0("mm39#", msum$mouse_chrom[1], ":", fmt_coord(msum$mouse_start[1] + 1)),
                     "+500", "+1,000", "+1,500", fmt_coord(msum$mouse_end[1]))) +
  geom_rect(data = msum, aes(xmin = mx(offset), xmax = mx(offset + 3 * n_codons),
                             ymin = y - 0.15, ymax = y + 0.15, fill = fill)) +
  geom_segment(data = mfr, aes(x = mx(pos + 1.5), xend = mx(pos + 1.5), y = y - 0.15,
                               yend = y + 0.15), linewidth = 0.35, colour = "black") +
  geom_text(data = msum, aes(x = x0 - 0.01, y = y, label = ifelse(is.na(full_length_orf_aa),
            paste0("frame ", frame, " (", n_stops, " stops)"), paste("frame", frame))), hjust = 1,
            size = DCR_SIZE$small / PT, family = DCR_FONT) +
  geom_text(data = msum[!is.na(full_length_orf_aa)], aes(x = mx(L / 2), y = y, label = lab),
            size = DCR_SIZE$small / PT, family = DCR_FONT, colour = "white", fontface = "bold") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 1), ylim = c(-1.75, 0.35), expand = FALSE) +
  theme_void(base_family = DCR_FONT) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(14, 4, 4, 4))

fig <- dcr_grid(list(gA, gB), labels = c("A", "B"), ncol = 1, rel_heights = c(3.6, 1),
                pad_top = 14, pad_left = 4)
save_dcr(fig, "Supp_GNAS_DCR_gene_model_mouse_ORFs", 7.2, 5.0, outdir)
write.csv(msum, file.path(outdir, "Supp_GNAS_mouse_frames_plotted.csv"), row.names = FALSE)
