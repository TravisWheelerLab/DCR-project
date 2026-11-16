# dcr_style.R -- house style for the DCR paper (Goubert et al., GBE revision)
#
# One style file for every revision figure. Source it at the top of each figure
# script:   source(file.path(script_dir, "dcr_style.R"))
#
# Rules encoded here (from the authors' house style):
#   * no titles on panels; all description goes in the figure legend text
#   * bold capital panel letters (A, B, C...) top-left; roman sub-parts "i." "ii."
#     in plain text, as in Figure 1
#   * ggplot2 classic theme: white background, L-shaped axes, no grid
#   * Helvetica/Arial-like sans font; legends at the bottom
#   * DCR orange #FFA800 vs greys (gray80 / grey60) for controls
#   * mouse status: DCR #FFA800 / not DCR orange4 / too short grey20 / not mapped #FFF0BD
#   * gene models: canonical frame orange, alternative frame purple, non-DCR exons steel blue
#   * significance as brackets with stars / ns
#
# Type sizes are measured from the published GBE/Figures/Figure2.pdf
# (706.7 x 718.3 pt page; Helvetica widths): panel letters 18 pt bold,
# axis titles 10 pt, tick labels 8 pt, legend / annotation text 8 pt.

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(scales)
})

DCR_FONT <- "Helvetica"

# ---- sizes (pt) ----
DCR_SIZE <- list(letter = 18, sub = 11, axis_title = 10, axis_text = 8,
                 legend = 8, annot = 8, small = 7)
PT <- ggplot2::.pt   # converts pt to ggplot "size" units for geom_text

# ---- palette ----
DCR_COL <- c(
  dcr          = "#FFA800",
  dcr_pale     = "#FFF0BD",   # D* in Fig 2C / "not mapped"
  ctrl_dark    = "grey60",
  ctrl_light   = "gray80",
  ctrl_darker  = "grey40",
  not_dcr      = "orange4",
  too_short    = "grey20",
  not_mapped   = "#FFF0BD",
  frame_can    = "#FFA800",   # canonical reading frame (orange)
  frame_alt    = "#8E44AD",   # alternative reading frame (purple)
  exon_nondcr  = "steelblue",
  dcr_highlight= "#F5EF1A",   # yellow DCR box in the Fig 2E zoom
  marker_fill  = "#F5EF1A",
  marker_line  = "#E8261C",
  neutral_line = "grey30"
)
MOUSE_STATUS_COL <- c("DCR" = "#FFA800", "not DCR" = "orange4",
                      "too short" = "grey20", "not mapped" = "#FFF0BD")
# Fig 2C exon classes (D, D*, ND)
EXON_CLASS_COL  <- c("D" = "#FFA800", "D*" = "#FFF0BD", "ND" = "gray80")
EXON_CLASS_LINE <- c("D" = "#FFA800", "D*" = "#F3DE8A", "ND" = "grey60")

# ---- theme ----
theme_dcr <- function(base_size = DCR_SIZE$axis_text, legend = "bottom") {
  theme_classic(base_size = base_size, base_family = DCR_FONT) %+replace%
    theme(
      axis.title   = element_text(size = DCR_SIZE$axis_title, colour = "black"),
      axis.title.y = element_text(size = DCR_SIZE$axis_title, angle = 90,
                                  margin = margin(r = 3)),
      axis.title.x = element_text(size = DCR_SIZE$axis_title, margin = margin(t = 3)),
      axis.text    = element_text(size = DCR_SIZE$axis_text, colour = "grey30"),
      axis.line    = element_line(colour = "black", linewidth = 0.5),
      axis.ticks   = element_line(colour = "black", linewidth = 0.4),
      legend.position = legend,
      legend.title = element_text(size = DCR_SIZE$legend),
      legend.text  = element_text(size = DCR_SIZE$legend),
      legend.key.size = unit(8, "pt"),
      legend.background = element_blank(),
      plot.title   = element_blank(),
      plot.subtitle= element_blank(),
      strip.background = element_blank(),
      strip.text   = element_text(size = DCR_SIZE$axis_text, colour = "black"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin  = margin(4, 4, 4, 4)
    )
}

# ---- significance ----
p_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

# Bracket layer(s). df: x1, x2 (numeric positions), y (bracket height), label;
# tip = downward tick length in data units.
sig_bracket <- function(df, tip, text_size = DCR_SIZE$axis_title, gap = 0.3 * tip,
                        lw = 0.4) {
  segs <- rbind(
    data.frame(x = df$x1, xend = df$x2, y = df$y, yend = df$y),
    data.frame(x = df$x1, xend = df$x1, y = df$y, yend = df$y - tip),
    data.frame(x = df$x2, xend = df$x2, y = df$y, yend = df$y - tip))
  lab <- data.frame(x = (df$x1 + df$x2) / 2, y = df$y + gap, label = df$label)
  # stars sit closer to the bracket than "ns" (stars have no descender space)
  list(
    geom_segment(data = segs, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, linewidth = lw, colour = "black"),
    geom_text(data = lab, aes(x = x, y = y, label = label), inherit.aes = FALSE,
              vjust = 0, size = text_size / PT, family = DCR_FONT)
  )
}

# ---- panel letters / composition ----
# Bold capital letters, top-left, as in Figures 1-6.
# Each plot gets extra top/left margin so the letter never collides with tick labels.
dcr_grid <- function(plots, labels = LETTERS[seq_along(plots)], pad_top = 14, pad_left = 10, ...) {
  plots <- lapply(plots, function(p) if (inherits(p, "ggplot"))
    p + theme(plot.margin = margin(pad_top, 4, 4, pad_left)) else p)
  cowplot::plot_grid(plotlist = plots, labels = labels, label_size = DCR_SIZE$letter,
                     label_fontfamily = DCR_FONT, label_fontface = "bold",
                     hjust = 0, vjust = 1, label_x = 0, label_y = 1, ...)
}
# Plain roman sub-part tag ("i.", "ii.") placed top-left of a sub-plot.
sub_tag <- function(p, tag, x = 0, y = 1, size = DCR_SIZE$sub) {
  cowplot::ggdraw(p) + cowplot::draw_label(tag, x = x, y = y, hjust = 0, vjust = 1,
                                           size = size, fontfamily = DCR_FONT)
}

# Grid of aligned sub-plots with plain roman sub-part tags ("i.", "ii.") top-left
# (Figure 1 style). Aligns axes across sub-plots.
sub_grid <- function(plots, tags, ...) {
  cowplot::plot_grid(plotlist = plots, labels = tags, label_size = DCR_SIZE$sub,
                     label_fontfamily = DCR_FONT, label_fontface = "plain",
                     hjust = 0, vjust = 1, label_x = 0, label_y = 1, ...)
}

# ---- output ----
# Always vector PDF + 300-dpi PNG. The PDF uses the base pdf() device with
# Helvetica (cairo is unavailable without XQuartz on this machine); Greek letters
# go through plotmath (Symbol font), so keep labels ASCII or plotmath. Fonts are
# embedded afterwards with Ghostscript when it is found.
save_dcr <- function(plot, stem, width, height, outdir = "figures") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  pdf_f <- file.path(outdir, paste0(stem, ".pdf"))
  png_f <- file.path(outdir, paste0(stem, ".png"))
  ggsave(pdf_f, plot = plot, width = width, height = height, units = "in",
         device = grDevices::pdf, useDingbats = FALSE, family = DCR_FONT, bg = "white")
  gs <- Sys.which(c("gs", "/opt/homebrew/bin/gs", "/usr/local/bin/gs"))
  gs <- gs[nzchar(gs)][1]
  if (!is.na(gs)) {
    Sys.setenv(R_GSCMD = gs)
    tmp <- paste0(pdf_f, ".tmp.pdf")
    ok <- tryCatch({ grDevices::embedFonts(pdf_f, outfile = tmp,
                       options = "-dPDFSETTINGS=/prepress -dEmbedAllFonts=true"); TRUE },
                   error = function(e) FALSE)
    if (ok && file.exists(tmp)) file.rename(tmp, pdf_f)
  }
  ggsave(png_f, plot = plot, width = width, height = height, units = "in",
         dpi = 300, device = ragg::agg_png, bg = "white")
  invisible(c(pdf_f, png_f))
}

# Read a CSV, preferring the first existing path in `paths` (e.g. a v3 that may
# appear later); returns the data with attribute "source" = path used.
read_first <- function(paths, ...) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) stop("none of the inputs exist: ", paste(paths, collapse = ", "))
  d <- utils::read.csv(hit, check.names = FALSE, ...)
  attr(d, "source") <- hit
  message("reading ", hit)
  d
}

fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
