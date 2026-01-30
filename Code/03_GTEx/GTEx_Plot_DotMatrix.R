#!/usr/bin/env Rscript

#library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
# library(ggpubr)
# library(png)
# library(jpeg)
library(ggh4x)
library(tidyr)
library(stringr)
library(tidyr)
library(pheatmap)
library(reshape2)

# Datasets --> HERE WE USE THE DATA LOADED FROM GTEx_data_plot.R

# sed 's/ : /:/g' output/03_GTEx/output_DCE_All_frames.tsv > output/03_GTEx/output_DCE_All_frames_simple.tsv
main_data <- read.table("output/03_GTEx/output_DCE_All_frames_simple.tsv", header = TRUE) #read.table("data/DCE_GTEx_DE_xnomial.out", header = TRUE)
names(main_data)<-c("pairID", "out_p", "n_tissues", "n_comps", "out_frames", "out_maxdif", "out_pattern", "best_pair")
genes<-read.table("output/03_GTEx/gene_chr_pos_DCRid.tsv")
names(genes)<-c("gname", "chr", "pos", "DCE_ID")
main_data[c("DCE_ID")]<-str_split_fixed(main_data$pairID, pattern =  "_", n = 2)[,1]
main_data<-merge(main_data, genes, by = "DCE_ID")

main_data$key <- rownames(main_data)
minP <- min(main_data[main_data$out_p > 0,]$out_p)
main_data$out_p[main_data$out_p == 0] <- minP
main_data<-separate_wider_delim(main_data, cols = pairID, delim = "_", names = c("id", "end"))
main_data$pairID <- paste(main_data$DCE_ID, main_data$end, sep = "_")
main_data$y <- -log10(main_data$out_p)
main_data$x <- paste(main_data$pairID, main_data$gname, sep = "|")

associated_data <- read.table("output/03_GTEx/AllFrames.pairedCounts_long.tsv", header = TRUE, comment.char = "")
normsFactors<-read.table("output/03_GTEx/GTExMetadataUpdated.txt", h = T, comment.char = "")[,c(1,3,6)]
associated_data<-merge(associated_data[,-11], normsFactors, by = "sample")

### attempt at heatmap
perIsoDCR<-read.table("output/01_Description/isoforms_to_DCR_metadata_canon.tsv", h = T)
updated_frame_assignment<-read.table("TableS1.tsv", h = T)
updated_frame_assignment<-updated_frame_assignment[,c("isoform", "DCR", "frame", "canonical_frame", "frame_status")]
names(updated_frame_assignment)<-c("isoform", "DCR", "frame", "new_canonical_frame", "new_frame_status")
perIsoDCR$Iso_DCR<-paste(perIsoDCR$isoform, perIsoDCR$DCR, sep = "_")
updated_frame_assignment$Iso_DCR<-paste(updated_frame_assignment$isoform, updated_frame_assignment$DCR, sep = "_")
  
perIsoDCR<-merge(perIsoDCR, updated_frame_assignment[,c("Iso_DCR", "new_canonical_frame", "new_frame_status")], by = "Iso_DCR")
perIsoDCR$canonical_frame<-perIsoDCR$new_canonical_frame
perIsoDCR$frame_status<-perIsoDCR$new_frame_status
perIsoDCR<-perIsoDCR[,-c(9,10)]

perIsoDCR$DCR_frame<-paste(perIsoDCR$DCR, perIsoDCR$frame, sep = "_")
# fix the one isoform with no attributed frame
perIsoDCR[perIsoDCR$frame_status == "undetermined",]

# for the figure, we have to assign a canonical frame even if isoform 1 is not involved:
# - canonical goes to frame used by most isoforms
# - if tie, goes with longest isoform (non-truncated form), probably will have to decide manually.
# however we will only do that for the one in the candidates
# unique(perIsoDCR[perIsoDCR$frame_status == "undetermined" & perIsoDCR$DCR %in% main_data$DCE_ID,]$DCR)
# [1] 1196  669 but 669 is NS --> we only assign 1169
perIsoDCR[perIsoDCR$DCR == "1196",]$frame_status<-c("canonical", "alternative", "canonical")

associated_data[associated_data$pairID %in% main_data[main_data$out_pattern == "***",]$pairID,] -> dis.cand
dis.cand$frame<-sub("frame", "", dis.cand$frame)
dis.cand %>%
  group_by(tissue, ID, frame) %>% summarize(total_counts = sum(count), gene = gene) -> ass.dat.sumzed
ass.dat.sumzed<-unique(ass.dat.sumzed)
ass.dat.sumzed$key <- paste(ass.dat.sumzed$ID, ass.dat.sumzed$frame, sep = "_")
perIsoDCR$key <- paste(perIsoDCR$DCR, perIsoDCR$frame, sep = "_")

ass.dat.sumzed$newID<-paste(ass.dat.sumzed$ID, ass.dat.sumzed$gene, sep = "-")

ass.dat.sumzed.status<-unique(merge(ass.dat.sumzed, perIsoDCR[,c("frame_status", "key")], by = "key"))

ass.dat.sumzed.status %>%
  group_by(tissue, newID) %>%
  summarise(
    counts_canonical = sum(total_counts[frame_status == "canonical"], na.rm = TRUE),
    total_counts      = sum(total_counts, na.rm = TRUE),
    .groups = "drop"   # drop all grouping here
  ) %>%
  group_by(newID) %>%                             # regroup by ID for the max
  mutate(max_total_ID = max(total_counts, na.rm = TRUE)) %>%
  ungroup() -> pre.mat.counts

pre.mat.counts$ratio <- pre.mat.counts$counts_canonical/pre.mat.counts$total_counts
pre.mat.counts$scaled_ratio <- 2 * (pre.mat.counts$ratio - 0.5) # center at 0 and [-1;1]
pre.mat.counts$norm.factor <- pre.mat.counts$total_counts/pre.mat.counts$max_total_ID

(1/2)+0.5
(0.25-0.5)*2

### Below is for the heatmap
pre.mat.counts$norm.ratios <- pre.mat.counts$scaled_ratio * pre.mat.counts$norm.factor
mat.counts <- dcast(pre.mat.counts, formula = tissue ~ newID)
rownames(mat.counts)<-mat.counts$tissue

### Other idea, a grid of dots
# we will reuse the heatmap function to get the clustering order
DCRdotMatrix_df<-pre.mat.counts
# export for app

write.table(associated_data, "Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR_Explorer/data/associated_data_new.tsv",
            row.names = F, quote = F)
# gene-level clustering based on frame usage
ratio_mat<-dcast(DCRdotMatrix_df[,c(1,2,7)], newID ~ tissue)
ratio_clust<-hclust(dist(ratio_mat, method = "euclidean"), method = "complete")
saveRDS(ratio_clust, file = "Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR_Explorer/data/ratio_clust.rds")
rows_order<-ratio_clust$order
DCRgenes_order<-unique(DCRdotMatrix_df$newID)[rows_order]
DCRdotMatrix_df$newID<-factor(DCRdotMatrix_df$newID, levels = DCRgenes_order)

name_cor<-main_data[,c("pairID", "DCE_ID", "gname")]
name_cor$newID<-paste(name_cor$DCE_ID, name_cor$gname, sep = "-")
DCRdotMatrix_df<-merge(DCRdotMatrix_df, name_cor, by = "newID")

write.table(DCRdotMatrix_df, "Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR_Explorer/data/DCRdotMatrix_df_new.tsv",
            row.names = F, quote = F)

# circle with border
dotPlot<-ggplot(DCRdotMatrix_df)+
  geom_point(aes(x = factor(tissue), y = newID,
                 fill = scaled_ratio, size = norm.factor, alpha = is.na(scaled_ratio)),
             pch = 21, col = "black", stroke =  0.1)+
  scale_fill_gradientn(colors = rev(c("orange", "#faf4e8", "#e4dae1", "purple")),
                       values = scales::rescale(c(-1, -0.01, 0.01, 1)),
                       breaks = c(-1,0, 1), labels = c("alternative",
                                                       "1:1",
                                                     "canonical"),
                       name = "DCR\nframe proportion")+
  # scale_fill_gradientn(name = "frame ratio", colors = rev(c("orange", "white", "purple")),
  #                       values = scales::rescale(c(-1, 0, 1)))+
  scale_alpha_manual(values = c(1,0))+
  scale_size_continuous(
    name = "relative gene\nexpression\n(across tissues)",
    range = c(1, 4),       # visual range of point sizes
    limits = c(0, 1)      # data range (values outside are clipped)
  ) +
  scale_y_discrete(position = "right") +
  xlab("GTEx tissue")+
  ylab("DCR gene")+
  theme_minimal(base_size = 8)+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        plot.margin = margin(0, 0, 0, 0)
  )+
  guides(alpha = "none")

DCRgene_dendro<-ggdendro::dendro_data(ratio_clust, type = "rectangle")
DgDendroPlot<-ggplot() +
  geom_segment(
    data = DCRgene_dendro$segments,
    aes(x = y, y = x, xend = yend, yend = xend),
    linewidth = 0.3, colour = "grey"
  ) +
  scale_x_reverse(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0.5, length(rows_order) + 0.5)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

gtexMetadata<-read.table("GTExMetadataUpdated.txt", header=TRUE, comment.char = "")
# bar with colors
#levels(factor(DCRdotMatrix_df$tissue))
tissuecolor <- unique(gtexMetadata[c("tissue", "color")])
tissuecolor<-tissuecolor[order(tissuecolor$tissue),]
colorkey<-ggplot(tissuecolor)+
  geom_tile(aes(x = factor(tissue), y = 1, fill = factor(tissue)), show.legend = F)+
  theme_minimal(base_size = 8)+
  scale_fill_manual(values = tissuecolor$color)+
  xlab("GTEx tissue")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.ticks.x = element_line(colour = "grey"),
        panel.grid = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

al_dendro<-cowplot::ggdraw(cowplot::align_plots(DgDendroPlot, dotPlot, align = "h")[[1]])
al_cbar<-cowplot::ggdraw(cowplot::align_plots(colorkey, dotPlot, align = "v")[[1]])

cowplot::plot_grid(al_dendro, dotPlot, NULL, al_cbar, ncol = 2, rel_widths = c(0.1, 0.9), rel_heights = c(0.83, 0.17))

length(unique(DCRdotMatrix_df$newID))

# test interactive

dotPlotLy<-ggplotly(dotPlot)
dotPlotLy<-dotPlotLy %>% layout(yaxis = list(side = "right", automargin = TRUE),
                                legend = list(
                                  x = 2,           # move further right of the plotting area
                                  xanchor = "top",  # anchor legend box to its left edge
                                  y = 1))

for (i in seq_along(dotPlotLy$x$data)) {
  tr <- dotPlotLy$x$data[[i]]
  
  # scatter/markers with per-trace colors
  if (!is.null(tr$marker) && (!is.null(tr$marker$colorscale) || !is.null(tr$marker$color))) {
    dotPlotLy$x$data[[i]]$marker$coloraxis <- "coloraxis"
    dotPlotLy$x$data[[i]]$marker$colorscale <- NULL
    dotPlotLy$x$data[[i]]$marker$colorbar  <- NULL
    dotPlotLy$x$data[[i]]$showscale <- NULL
  }
  
  # heatmap/contour/imshow-like traces
  if (!is.null(tr$colorscale) || !is.null(tr$z)) {
    dotPlotLy$x$data[[i]]$coloraxis <- "coloraxis"
    dotPlotLy$x$data[[i]]$colorscale <- NULL
    dotPlotLy$x$data[[i]]$colorbar   <- NULL
    dotPlotLy$x$data[[i]]$showscale  <- NULL
  }
}

dotPlotLy <- dotPlotLy %>% layout(
  coloraxis = list(
    colorbar = list(
      x = 1.10,        # push to the right of the plotting area
      xanchor = "left",
      y = 0.5,         # center vertically
      len = 0.9,       # optional: length as a fraction of plot height
      thickness = 14   # optional: width in px
    )
  ),
  margin = list(r = 180)  # ensure there’s room to render the colorbar
)



# scale_fill_gradientn(colors = rev(c("orange", "#faf4e8", "#e4dae1", "purple")),
#                      values = scales::rescale(c(-1, -0.01, 0.01, 1)),
#                      breaks = c(-1,0, 1), labels = c("alternative",
#                                                      "1:1",
#                                                      "canonical")
#                      
                     

## 1) Decide the numeric limits and midpoints you used in ggplot
limits <- c(-1, 1)             # e.g., limits= in your ggplot scale
mid1   <- -.01                 # your first midpoint
mid2   <- 0.01                 # your second midpoint
cols   <- c("low"="orange", "mid1"="#faf4e8", "mid2"="#e4dae1", "high"="purple")

## 2) Build a Plotly colorscale with *positions in [0,1]*
pos <- rescale(c(limits[1], mid1, mid2, limits[2]), to = c(0, 1), from = limits)
plotly_colorscale <- list(
  list(pos[1], cols["low"]),
  list(pos[2], cols["mid1"]),
  list(pos[3], cols["mid2"]),
  list(pos[4], cols["high"])
)

## 3) Set breaks & labels to match ggplot
breaks <- c(-1, 0, 1)          # same numeric breaks as ggplot
labels <- c("canonical", "1:1", "alternative")          # same labels as ggplot

## 4) Apply on the shared coloraxis
dotPlotLy <- dotPlotLy %>% layout(
  coloraxis = list(
    cmin = limits[1],
    cmax = limits[2],
    colorscale = plotly_colorscale,
    # set to TRUE if you used scale_*_reverse()
    reversescale = FALSE,
    colorbar = list(
      x = 1.10, xanchor = "left",
      y = 0.9, len = 0.2, thickness = 14,
      title = list(text = "DCR\nframe proportion"),
      tickmode = "array",
      tickvals = breaks,
      ticktext = labels
    )
  ),
  margin = list(r = 180)
)

# --- Helper: add a size legend under the colorbar ----------------------------
add_size_legend_vertical <- function(fig,
                                     sizes_px = NULL,   # e.g., c(8,14,22,30)
                                     labels   = NULL,   # e.g., c("10","50","100","200")
                                     title    = "Point size",
                                     cb_x     = NULL,   # x position (default: next to colorbar)
                                     y_top    = 0.20,   # start below the colorbar top
                                     gap      = 0.06,   # vertical spacing between points
                                     align    = "center" # "center" or "left"
) {
  if (is.null(cb_x)) {
    cb_x <- tryCatch(fig$x$layout$coloraxis$colorbar$x, error = function(e) 1.10)
    if (is.null(cb_x)) cb_x <- 1.10
  }
  
  # derive representative sizes if none provided
  if (is.null(sizes_px)) {
    idx <- which(sapply(fig$x$data, function(tr) !is.null(tr$marker$size)))[1]
    if (length(idx)) {
      s <- as.numeric(fig$x$data[[idx]]$marker$size)
      s <- s[is.finite(s)]
      if (length(s)) {
        q <- stats::quantile(s, c(0.25, 0.5, 0.75), names = FALSE)
        sizes_px <- as.numeric(unique(round(q)))
      }
    }
    if (is.null(sizes_px) || !length(sizes_px)) sizes_px <- c(8,14,22)
  }
  if (is.null(labels)) labels <- as.character(sizes_px)
  
  x0 <- cb_x + 0.04  # move right of colorbar
  y0 <- y_top
  
  anns <- list(
    list(
      x = x0, y = y0 + 0.05, xref = "paper", yref = "paper",
      text = title, showarrow = FALSE, xanchor = align, font = list(size = 12)
    )
  )
  
  for (i in seq_along(sizes_px)) {
    yi <- y0 - (i - 1) * gap
    anns[[length(anns)+1]] <- list(
      x = x0, y = yi, xref = "paper", yref = "paper",
      text = "●", showarrow = FALSE, xanchor = "center",
      font = list(size = sizes_px[i])
    )
    anns[[length(anns)+1]] <- list(
      x = x0 + 0.05, y = yi - 0.01, xref = "paper", yref = "paper",
      text = labels[i], showarrow = FALSE, xanchor = "left",
      font = list(size = 10)
    )
  }
  
  fig %>% plotly::layout(annotations = c(fig$x$layout$annotations %||% list(), anns))
}

dendroPlotLy<-ggplotly(DgDendroPlot)%>% layout(
  xaxis = list(showgrid = FALSE, zeroline = FALSE,
               showticklabels = FALSE, ticks = "", title = "", showline = FALSE),
  yaxis = list(showgrid = FALSE, zeroline = FALSE,
               showticklabels = FALSE, ticks = "", title = "", showline = FALSE),
  showlegend   = FALSE,
  plot_bgcolor = "rgba(0,0,0,0)",
  paper_bgcolor= "rgba(0,0,0,0)"
)



dotPlotLy<-add_size_legend_vertical(dotPlotLy, title = "normalized gene\nexpression",sizes_px = c(0.3, 15, 30),
                               labels = c(0.1, 0.5, 1), y = 0.6, gap = 0.03)
subplot(dendroPlotLy, dotPlotLy,widths = c(0.2, 0.8))
#############################################################################################

# ART
ggplot(DCRdotMatrix_df)+
  geom_point(aes(x = factor(tissue), y = newID, col = scaled_ratio, size = norm.factor))+
  # scale_color_gradientn(colors = rev(c("orange", "bisque1", "thistle2", "purple")),
  #                       values = scales::rescale(c(-1, -0.01, 0,01, 1)))+
  scale_color_gradientn(name = "frame ratio", colors = rev(c("orange", "white", "purple")),
                        values = scales::rescale(c(-1, 0, 1)))+
  scale_size_continuous(
    name = "normalized expression\n(per DCR-gene)",
    range = c(1, 4),       # visual range of point sizes
    limits = c(0, 1)      # data range (values outside are clipped)
  ) +
  xlab("GTEx tissue")+
  ylab("DCR-gene")+
  theme_void()+
  theme(axis.text.x =  element_blank(),legend.position = "none")
