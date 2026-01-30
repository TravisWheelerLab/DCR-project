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

# Datasets
# sed 's/ : /:/g' output/03_GTEx/output_DCE_All_frames.tsv > output/03_GTEx/output_DCE_All_frames_simple.tsv
main_data <- read.table("Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/03_GTEx/output_DCE_All_frames_simple.tsv", header = TRUE) #read.table("data/DCE_GTEx_DE_xnomial.out", header = TRUE)
names(main_data)<-c("pairID", "out_p", "n_tissues", "n_comps", "out_frames", "out_maxdif", "out_pattern")
genes<-read.table("Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/03_GTEx/gene_chr_pos_DCRid.tsv")
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

associated_data <- read.table("Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/03_GTEx/AllFrames.pairedCounts_long.tsv", header = TRUE, comment.char = "")
normsFactors<-read.table("Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/03_GTEx/GTExMetadataUpdated.txt", h = T, comment.char = "")[,c(1,3,6)]
associated_data<-merge(associated_data[,-11], normsFactors, by = "sample")
# head(main_data)
names(DCRarc)[1]<-"DCE_ID"
plot_main_data<-merge(main_data, DCRarc, by = "DCE_ID")

# ggplot(plot_main_data, aes(x = x, y = as.numeric(out_maxdif),
#                       text = paste("pairID:", pairID,
#                                    "\nmin. P-value (Bonf.):", signif(out_p, 3),
#                                    "\njunction end:", end, size = y)))+
#   geom_point(aes(size = -log10(out_p), col = out_pattern == "***", shape = `Downstream:C-Terminal` == TRUE))+
#   #geom_point(aes(size = -log10(out_p), col = chr))+
#   scale_alpha_manual(name = "frame switch pattern", values = c(1, .2))+
#   xlab("DCE_junction|gene") +
#   ylab("Maximum frame expression difference across all \npairwise tissues comparisons (ratio)") +
#   scale_size(name = "-log10(P-value)", range = c(0.2, 2.2)) +
#   #scale_color_manual(name = "Tissue x Junction Interraction", values = c("grey80", "green3")) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90),
#         plot.background = element_blank(),
#         legend.position = "bottom")
# 

plot_main_data$Termini<-factor(ifelse(plot_main_data$`Downstream:C-Terminal` == TRUE & plot_main_data$`Upstream:N-Terminal` == TRUE, "CN", 
                               ifelse(plot_main_data$`Downstream:C-Terminal` == FALSE & plot_main_data$`Upstream:N-Terminal` == FALSE, "Internal",
                                      ifelse(plot_main_data$`Downstream:C-Terminal` == TRUE & plot_main_data$`Upstream:N-Terminal` == FALSE, "C-Terminal", "N-Terminal"))),
                               levels = c("C-Terminal", "N-Terminal", "CN", "Internal"))

ggplot(plot_main_data, aes(x = as.numeric(out_maxdif), y = -log10(out_p),
                          text = paste0("DCE_ID: ", DCE_ID,
                                        "<br>out_p: ", signif(out_p, 3),
                                        "<br>maxdif: ", signif(out_maxdif, 3),
                                        "<br>pattern: ", out_pattern,
                                        "<br>frames: ", out_frames)
))+
  geom_point(aes(shape = out_pattern == "***",
                 col = Termini))+
  #geom_point(aes(size = -log10(out_p), col = chr))+
  scale_alpha_manual(name = "frame switch pattern", values = c(1, .2))+
  xlab("Maximum expression ratio between frames") +
  ylab("Minimum pairwise P-value (Bonferronni)\n-log10(P-value)") +
  scale_shape_manual(name = "Disordinal interaction", values = c(1, 19))+
  #scale_size(name = "-log10(P-value)", range = c(0.2, 2.2)) +
  #scale_color_manual(name = "Tissue x Junction Interraction", values = c("grey80", "green3")) +
  scale_color_manual(name = "Termini", values = c("salmon", "royalblue1", "purple", "grey")) +
  scale_y_continuous(breaks = c(0, 100, 200, 300))+
  coord_cartesian(ylim = c(0, 301))+
  theme_classic() +
  theme(plot.background = element_blank(),
        legend.position = "bottom")

## stats for paper
length(unique(plot_main_data$DCE_ID))
length(unique(plot_main_data[plot_main_data$out_pattern == "***",]$DCE_ID))
table(plot_main_data[plot_main_data$out_pattern == "***",]$Termini)

## per DCR plot function
# pair_id = "1064_left"
# r_threshold = 2
plot_function <- function(pair_id, data = associated_data, table_out = FALSE, r_threshold = 2, TMMfactors = NULL, shortID = FALSE, filter = TRUE, raw = TRUE) {
  done <- data[data$pairID == pair_id, ]
  if (nrow(done) == 0) {
    stop("No data found for the selected pairID.")
  }
  
  if (filter) {
    sRaw <- done %>%
      group_by(tissue) %>%
      summarize(meanRPM = mean(count, na.rm = TRUE), .groups = "drop")
    
    t_to_keep <- unique(sRaw[sRaw$meanRPM > r_threshold, ]$tissue)
    if (length(t_to_keep) > 1) {
      done <- done[done$tissue %in% t_to_keep, ]
      if (!all(c("lib_size", "NormFactors") %in% colnames(done))) {
        stop("Columns 'lib_size' and 'NormFactors' are required for CPM calculation.")
      }
      done$CPM <- (done$count / (done$lib_size * done$NormFactors) * 1e6)
      mean.t.CPM <- done %>%
        group_by(tissue) %>%
        summarize(meanCPM = mean(CPM, na.rm = TRUE), .groups = "drop")
      CPMthreshold <- signif(min(mean.t.CPM$meanCPM), 2)
    } else {
      stop("Not enough tissues pass the filtering threshold.")
    }
  }
  
  myCols <- done$color
  names(myCols) <- done$tissue
  myTiss <- done$tissue
  ct_table <- data.frame(color = myCols, tissue = myTiss)
  ct_table_uniq <- distinct(ct_table)
  labelcol<-ct_table_uniq[order(ct_table_uniq$tissue),]$color
  strips <- strip_themed(background_x = elem_list_rect(fill = labelcol))
  
  if (raw) {
    p <- ggplot(done, aes(x = factor(frame), y = count)) +
      geom_hline(aes(yintercept = r_threshold), linetype = 2) +
      #geom_label(aes(label = tissue, fill = tissue))+
      #scale_fill_manual(values = labelcol)+
      geom_jitter(aes(col = frame, shape = end), width = 0.2) +
      geom_boxplot(aes(col = frame), width = 0.4, fill = "white", alpha = 0.5, outlier.alpha = 0) +
      scale_color_manual(values = c("purple", "orange", "firebrick"))+
      xlab("GTEx tissue/frame") +
      ylab("Raw counts") +
      ggtitle(paste(
        "DCE:", unique(done$ID),
        "| Gene:", unique(done$gene),
        "| End:", unique(done$end),
        "| Raw reads threshold:", r_threshold,
        "| CPM threshold:", CPMthreshold)) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        strip.text.x = element_text(angle = 90)
      )#+
    #facet_wrap2(~tissue, nrow = 1, strip = strips)
  } else {
    p <- ggplot(done, aes(x = factor(frame), y = CPM)) +
      geom_hline(aes(yintercept = CPMthreshold), linetype = 2) +
      geom_jitter(aes(col = frame, shape = end), width = 0.2) +
      geom_boxplot(aes(col = frame), width = 0.4, fill = "white", alpha = 0.5, outlier.alpha = 0) +
      scale_color_manual(values = c("purple", "orange", " "))+
      xlab("GTEx tissue/frame") +
      ylab("Normalized Counts (CPM-TMM)") +
      ggtitle(paste(
        "DCE:", unique(done$ID),
        "| Gene:", unique(done$gene),
        "| End:", unique(done$end),
        "| Raw reads threshold:", r_threshold,
        "| CPM threshold:", CPMthreshold)) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        strip.text.x = element_text(angle = 90)
      )#+
    #facet_wrap2(~tissue, nrow = 1, strip = strips)
  }
  
  if (shortID) {
    p <- p + facet_wrap2(~ShortName, nrow = 1, strip = strips)#facet_grid(~ShortName)
  } else {
    p <- p + facet_wrap2(~tissue, nrow = 1, strip = strips)#facet_grid(~tissue)
  }
  
  return(p)
}

# Create ggplot object with 'text' aesthetic
# Main plot
library(ggrepel)

p1 <- ggplot(main_data, aes(
  x = as.numeric(out_maxdif),
  y = -log10(out_p)
)) +
  geom_point(aes(
    size = -log10(out_p),
    color = out_pattern == "***",
    shape = factor(out_frames),
  )) +
  geom_label_repel(data = main_data[main_data$out_pattern == "***",], aes(label = pairID))+
  scale_size(name = "-log10(P-value)", range = c(0.2, 2.2)) +
  scale_color_manual(
    name = "Tissue x Junction Interaction",
    values = c("grey80", "green3")
  ) +
  xlab("maximum difference between two frames") +
  ylab("-log10(P-value)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90),
    plot.background = element_blank(),
    legend.position = "bottom"
  )

clemsPicks<-c("359_left",
              "824_left",
              "1046_left",
              "1264_left",
              "1270_left",
              "1768_left",
              "1769_left",
              "1784_left",
              "1875_left",
              "7_left",
              "160_left",
              "167_left",
              "247_left",
              "361_left",
              "458_left",
              "472_left",
              "481_left",
              "520_left",
              "523_left",
              "591_left",
              "673_left",
              "762_left",
              "774_left",
              "883_right",
              "907_right",
              "987_left",
              "1003_left",
              "1064_left",
              "1086_left",
              "1215_left",
              "1219_left",
              "1278_left",
              "1287_left",
              "1573_left",
              "1630_left",
              "1701_left",
              "1812_left",
              "2012_left")

# plotly for quick lookup
plot_ly(data = main_data,
        x = ~as.numeric(out_maxdif),
        y = ~-log10(out_p),
        text = ~paste(pairID),
        symbol = ~(out_pattern == "***"),
        color = ~(pairID %in% clemsPicks),
        type = "scatter"
        )

# better
library(plotly)

library(plotly)

# Prepare data: make sure clems_pick is "Yes"/"No" for coloring
main_data$clems_pick <- ifelse(main_data$pairID %in% clemsPicks, "Yes", "No")

# Optionally clean up out_pattern (if needed)
main_data$out_pattern <- factor(main_data$out_pattern)

# Define shape and color mappings
shape_map <- c("***" = "circle", "NS" = "circle-open")
color_map <- c("Yes" = "limegreen", "No" = "grey")

# Start plotly object
p <- plot_ly()

# Add traces manually by shape (out_pattern) and color (clems_pick)
for (pat in unique(main_data$out_pattern)) {
  for (pick in c("Yes", "No")) {
    subset_df <- main_data[main_data$out_pattern == pat & main_data$clems_pick == pick, ]
    
    if (nrow(subset_df) > 0) {
      p <- add_trace(
        p,
        data = subset_df,
        x = ~as.numeric(out_maxdif),
        y = ~-log10(out_p),
        type = 'scatter',
        mode = 'markers',
        marker = list(
          symbol = shape_map[as.character(pat)],
          color = color_map[pick],
          size = 7
        ),
        name = paste("Frame Switch:", pat, "| Clem's Pick:", pick),
        legendgroup = paste0("group_", pat, "_", pick),
        showlegend = TRUE,
        text = ~paste("pairID:", pairID, "<br>Frame Switch:", out_pattern, "<br>Clem's Pick:", clems_pick),
        hoverinfo = 'text'
      )
    }
  }
}

# Layout: axis labels, bottom legend
p <- layout(
  p,
  xaxis = list(title = "Max frame expression difference"),
  yaxis = list(title = "-log10(P-value)"),
  legend = list(
    orientation = "h",
    x = 0.1,
    y = -0.1,
    xanchor = "left"
  )
)

p


#############
library(plotly)

symbols <- c(
  "circle", "square", "diamond", "cross", "x",
  "triangle-up", "triangle-down", "triangle-left", "triangle-right",
  "star", "hexagram", "hourglass", "bowtie"
)

plot_ly(
  x = seq_along(symbols), y = rep(1, length(symbols)),
  type = 'scatter', mode = 'markers+text',
  marker = list(symbol = symbols, size = 20),
  text = symbols, textposition = 'bottom center'
)




plot_function(pair_id = "1812_left", filter = T, raw = F) # gabbr1 brain specific inversion of pattern
plot_function(pair_id = "361_left", filter = T, raw = F) #


###

associated_data %>% group_by(sample) %>% summarize(tissue = unique(tissue),
                                                   lib_size = unique(lib_size),
                                                   ) -> GTEx_samples
write.table(GTEx_samples, "Google Drive/My Drive/Wheeler_Lab/2025/11-12_DCR_Paper/Sup/Table_S1.tsv", quote = F, row.names = F)
associated_data %>% group_by(tissue) %>% summarize(count = length(unique(sample))) -> sample_per_t
sample_per_t[sample_per_t$count < 10,]
