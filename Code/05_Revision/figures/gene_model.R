# gene_model.R -- Figure 2E-style gene model (isoform tracks + DCR zoom), drawn on a
# unit canvas with ggplot2. Sourced by the GNAS and GABBR1 figure scripts after
# dcr_style.R. All coordinates come from data/gm_<gene>_exons.csv (DOOVER BED12
# SwissProt isoform models) and data/gm_<gene>_dcr.csv.
#
# Elements reproduced from Figure 2E (MADD): coordinate ruler with
# "hg38#chrom:pos" first label, isoform accessions at left, steel-blue exons on a
# thin intron line with strand chevrons, DCR marker triangle (yellow, red outline),
# black frame around the DCR in the main tracks linked to a zoom box; in the zoom,
# one row per frame group, DCR-containing exon coloured by frame (canonical orange,
# alternative purple), yellow box marking the DCR interval; key at bottom left.
# No title (description goes in the legend).

suppressPackageStartupMessages(library(data.table))

collapse_iso <- function(acc) {
  # "Q9UBS5-1","Q9UBS5-2",.. -> "Q9UBS5 1-4"; mixed bases joined by "; "
  base <- sub("-[0-9]+$", "", acc); num <- as.integer(sub("^.*-", "", acc))
  out <- c()
  for (b in unique(base)) {
    n <- sort(num[base == b]); runs <- split(n, cumsum(c(1, diff(n) != 1)))
    s <- paste(vapply(runs, function(r) if (length(r) > 2) paste0(r[1], "-", r[length(r)])
                      else paste(r, collapse = ","), ""), collapse = ",")
    out <- c(out, paste(b, s))
  }
  paste(out, collapse = "; ")
}

fmt_coord <- function(x) formatC(x, format = "d", big.mark = ",")

frame_lab <- function(off, role) sprintf("%s%d (%s)", ifelse(off == 0, "", "+"), off, role)

# exons: data.frame(isoform, accession, start, end, strand, dcr_exon, frame_offset, frame_role)
# dcr: one-row data.frame(chrom, start, end, id, strand, gene)
# zoom: c(start, end) genomic window for the inset
# row_order: accessions in display order (top to bottom)
gene_model_plot <- function(exons, dcr, zoom, row_order, assembly = "hg38",
                            n_ticks = 5, label_w = 0.11, inset_x = c(0.42, 0.995),
                            min_exon_frac = 0.0012, chevron_step = 0.035,
                            group_label = NULL) {
  ex <- as.data.table(exons)
  ex <- ex[accession %in% row_order]
  g0 <- min(ex$iso_start); g1 <- max(ex$iso_end)
  pad <- 0.004 * (g1 - g0); g0 <- g0 - pad; g1 <- g1 + pad
  x0 <- label_w + 0.01; x1 <- 0.995
  gx <- function(g) x0 + (g - g0) / (g1 - g0) * (x1 - x0)
  nr <- length(row_order)
  ex[, row := match(accession, row_order)]
  ex[, y := -row * 0.42]
  strand <- ex$strand[1]
  eh <- 0.13                                # half exon height (main)
  # --- main tracks ---
  iso <- ex[, .(s = min(iso_start), e = max(iso_end), y = y[1]), by = accession]
  exr <- ex[, .(xmin = gx(start), xmax = pmax(gx(end), gx(start) + min_exon_frac), y)]
  chev <- iso[, {xs <- if (gx(e) - gx(s) > chevron_step)
                   seq(gx(s) + chevron_step / 2, gx(e) - chevron_step / 2, by = chevron_step)
                 else numeric(0)
                 .(x = xs, y = rep(y, length(xs)))}, by = accession]
  # drop chevrons that fall on exons
  on_ex <- function(x, yy) any(exr$y == yy & x >= exr$xmin - 0.004 & x <= exr$xmax + 0.004)
  chev <- chev[!mapply(on_ex, chev$x, chev$y)]
  ytop <- 0.05
  ticks <- round(seq(g0 + pad, g1 - pad, length.out = n_ticks))
  tick_lab <- fmt_coord(ticks); tick_lab[1] <- fmt_coord(ticks[1] + 1)  # 1-based start
  tick_lab[1] <- paste0(assembly, "#", dcr$chrom, ":\n", tick_lab[1])
  ymin_main <- min(ex$y) - 0.2
  # DCR frame in the main tracks
  dx <- c(gx(dcr$start), gx(dcr$end)); dx <- c(min(dx[1], mean(dx) - 0.004), max(dx[2], mean(dx) + 0.004))
  # --- inset ---
  grp <- ex[dcr_exon == TRUE, .(acc = list(sort(unique(accession))), off = frame_offset[1]),
            by = frame_role]
  grp <- grp[order(match(frame_role, c("canonical", "alternative")))]
  grp[, label := vapply(acc, collapse_iso, "")]
  if (!is.null(group_label)) grp[, label := group_label(label, frame_role)]
  ng <- nrow(grp)
  ins_top <- ymin_main - 0.35; rowh <- 0.36
  ins_rows <- ins_top - 0.28 - (seq_len(ng) - 1) * rowh
  ins_bot <- min(ins_rows) - 0.5
  ix0 <- inset_x[1] + 0.13; ix1 <- inset_x[2] - 0.01
  zx <- function(g) ix0 + (pmin(pmax(g, zoom[1]), zoom[2]) - zoom[1]) / diff(zoom) * (ix1 - ix0)
  ins_ex <- rbindlist(lapply(seq_len(ng), function(k) {
    rep_acc <- grp$acc[[k]][1]
    e <- ex[accession == rep_acc & end > zoom[1] & start < zoom[2]]
    e[, `:=`(yy = ins_rows[k], fill = ifelse(dcr_exon,
              ifelse(frame_role == "canonical", DCR_COL[["frame_can"]], DCR_COL[["frame_alt"]]),
              DCR_COL[["exon_nondcr"]]))]
    e
  }))
  ins_line <- rbindlist(lapply(seq_len(ng), function(k) {
    rep_acc <- grp$acc[[k]][1]; i <- iso[accession == rep_acc]
    data.table(x = zx(max(i$s, zoom[1])), xend = zx(min(i$e, zoom[2])), y = ins_rows[k])
  }))
  zt <- round(seq(zoom[1] + 0.12 * diff(zoom), zoom[2] - 0.12 * diff(zoom), length.out = 3))
  ins_yl <- ins_rows[ng] - 0.25
  # key (bottom-left)
  keys <- data.table(lab = c(frame_lab(grp$off, grp$frame_role), "non-DCR"),
                     fill = c(ifelse(grp$frame_role == "canonical", DCR_COL[["frame_can"]],
                                     DCR_COL[["frame_alt"]]), DCR_COL[["exon_nondcr"]]))
  keys[, `:=`(kx = 0.02, ky = ins_top - 0.32 - (seq_len(.N) - 1) * 0.2)]
  tsz <- DCR_SIZE$small / PT
  chev_char <- if (strand == "+") ">" else "<"

  ggplot() +
    # ruler
    annotate("segment", x = x0, xend = x1, y = ytop, yend = ytop, linewidth = 0.3) +
    annotate("segment", x = gx(ticks), xend = gx(ticks), y = ytop, yend = ytop + 0.05,
             linewidth = 0.3) +
    annotate("text", x = gx(ticks), y = ytop + 0.08, label = tick_lab, vjust = 0,
             hjust = c(0, rep(0.5, n_ticks - 2), 1), size = tsz, family = DCR_FONT,
             lineheight = 0.9) +
    # marker
    annotate("point", x = mean(c(gx(dcr$start), gx(dcr$end))), y = ytop + 0.5, shape = 25,
             size = 2.4, stroke = 0.6, fill = DCR_COL[["marker_fill"]],
             colour = DCR_COL[["marker_line"]]) +
    # isoforms
    geom_segment(data = iso, aes(x = gx(s), xend = gx(e), y = y, yend = y),
                 colour = DCR_COL[["exon_nondcr"]], linewidth = 0.35) +
    geom_text(data = chev, aes(x, y, label = chev_char), colour = DCR_COL[["exon_nondcr"]],
              size = 5 / PT, family = DCR_FONT, vjust = 0.42) +
    geom_rect(data = exr, aes(xmin = xmin, xmax = xmax, ymin = y - eh, ymax = y + eh),
              fill = DCR_COL[["exon_nondcr"]], colour = NA) +
    geom_text(data = iso, aes(x = label_w, y = y, label = accession), hjust = 1, size = tsz,
              family = DCR_FONT) +
    # DCR frame + connectors
    annotate("rect", xmin = dx[1], xmax = dx[2], ymin = ymin_main, ymax = -0.2,
             fill = NA, colour = "black", linewidth = 0.35) +
    annotate("segment", x = dx[1], xend = inset_x[1], y = ymin_main, yend = ins_top,
             linewidth = 0.3) +
    annotate("segment", x = dx[2], xend = inset_x[2], y = ymin_main, yend = ins_top,
             linewidth = 0.3) +
    # inset
    annotate("rect", xmin = inset_x[1], xmax = inset_x[2], ymin = ins_bot, ymax = ins_top,
             fill = "white", colour = "black", linewidth = 0.4) +
    annotate("rect", xmin = zx(dcr$start), xmax = zx(dcr$end),
             ymin = min(ins_rows) - 0.17, ymax = max(ins_rows) + 0.17,
             fill = DCR_COL[["dcr_highlight"]], colour = NA) +
    geom_segment(data = ins_line, aes(x = x, xend = xend, y = y, yend = y),
                 colour = DCR_COL[["exon_nondcr"]], linewidth = 0.5) +
    geom_rect(data = ins_ex, aes(xmin = zx(start), xmax = zx(end), ymin = yy - 0.085,
                                 ymax = yy + 0.085, fill = fill), colour = NA) +
    annotate("text", x = inset_x[1] + 0.01, y = ins_rows, label = grp$label, hjust = 0,
             size = tsz, family = DCR_FONT) +
    annotate("segment", x = zx(zt), xend = zx(zt), y = ins_yl, yend = ins_yl + 0.05,
             linewidth = 0.3) +
    annotate("text", x = zx(zt), y = ins_yl - 0.03, label = fmt_coord(zt), vjust = 1,
             size = tsz, family = DCR_FONT) +
    # key
    annotate("text", x = 0.02, y = ins_top - 0.12, label = "exon reading frame", hjust = 0,
             size = DCR_SIZE$annot / PT, family = DCR_FONT) +
    geom_rect(data = keys, aes(xmin = kx, xmax = kx + 0.03, ymin = ky - 0.06, ymax = ky + 0.06,
                               fill = fill)) +
    geom_text(data = keys, aes(x = kx + 0.04, y = ky, label = lab), hjust = 0, size = tsz,
              family = DCR_FONT) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0, 1), ylim = c(ins_bot - 0.05, ytop + 0.62), expand = FALSE) +
    theme_void(base_family = DCR_FONT) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(14, 4, 4, 4))
}
