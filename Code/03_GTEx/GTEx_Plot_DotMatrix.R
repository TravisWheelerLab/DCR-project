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
main_data <- read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/03_GTEx/output_DCE_All_frames_simple.tsv", header = TRUE) #read.table("data/DCE_GTEx_DE_xnomial.out", header = TRUE)
names(main_data)<-c("pairID", "out_p", "n_tissues", "n_comps", "out_frames", "out_maxdif", "out_pattern", "best_pair")
genes<-read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/03_GTEx/gene_chr_pos_DCRid.tsv")
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

associated_data <- read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/03_GTEx/AllFrames.pairedCounts_long.tsv", header = TRUE, comment.char = "")
normsFactors<-read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/03_GTEx/GTExMetadataUpdated.txt", h = T, comment.char = "")[,c(1,3,6)]
associated_data<-merge(associated_data[,-11], normsFactors, by = "sample")

### attempt at heatmap
perIsoDCR<-read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/New_TableS1.tsv", h = T)

# for the figure, we have to assign a canonical frame even if isoform 1 is not involved:
# - canonical goes to frame used by most isoforms
# - if tie, goes with longest isoform (non-truncated form), probably will have to decide manually.
# however we will only do that for the one in the candidates
# unique(perIsoDCR[perIsoDCR$frame_status == "undetermined" & perIsoDCR$DCR %in% main_data$DCE_ID,]$DCR)
# [1] 1196
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

write.table(associated_data, "/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR_Explorer/data/associated_data_new.tsv",
            row.names = F, quote = F)
# gene-level clustering based on frame usage
ratio_mat<-dcast(DCRdotMatrix_df[,c(1,2,7)], newID ~ tissue)
ratio_clust<-hclust(dist(ratio_mat, method = "euclidean"), method = "complete")
saveRDS(ratio_clust, file = "/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR_Explorer/data/ratio_clust.rds")
rows_order<-ratio_clust$order
DCRgenes_order<-unique(DCRdotMatrix_df$newID)[rows_order]
DCRdotMatrix_df$newID<-factor(DCRdotMatrix_df$newID, levels = DCRgenes_order)

name_cor<-main_data[,c("pairID", "DCE_ID", "gname")]
name_cor$newID<-paste(name_cor$DCE_ID, name_cor$gname, sep = "-")
DCRdotMatrix_df<-merge(DCRdotMatrix_df, name_cor, by = "newID")

write.table(DCRdotMatrix_df, "/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR_Explorer/data/DCRdotMatrix_df_new.tsv",
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

gtexMetadata<-read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/data/GTExMetadataUpdated.txt", header=TRUE, comment.char = "")
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

# length(unique(DCRdotMatrix_df$newID))
# 65

# ART #########################################################################################
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
