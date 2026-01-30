#!/usr/bin/env Rscript

### Distribution of DCE regions sizes

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: DCR_exon_length_distribution.r <DCEregions_bed_file>")
}
#dce_bed_file <- args[1]
dce_bed_file <- "Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb.bed"
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(cowplot)

# awk '{DCE=sub(/DCEid=/, "", $NF); Frames=sub(/DCE_region;nFrames=/,"",$3); print $1"\t"$4-1"\t"$5"\t"$NF"\t"$7"\t"$3}' DCE_overlapping_regions__FAST_091224.gff > DCE_overlapping_regions__FAST_091224.bed
DCEregions<-read.table(dce_bed_file)
names(DCEregions)<-c("Index", "chr", "start", "end", "DCEindex", "direction", "gene", "frames")
DCEregions$length<- DCEregions$end - DCEregions$start
exons_per_dce <- as.data.frame(tapply(DCEregions$DCEindex, DCEregions$DCEindex, length))
exons_per_dce$DCEindex <- rownames(exons_per_dce)
names(exons_per_dce) <- rev(c("DCEindex", "exons_ct"))
DCEregions <- merge(DCEregions, exons_per_dce, by = "DCEindex")
#DCEregions<-DCEregions[DCEregions$DCEindex %in% unique(humDCE2F_noITR2F$Index),]
DCEregions_Final <- DCEregions %>% group_by(DCEindex) %>%
  summarize(DCEindex = unique(DCEindex),
            exons_ct = unique(exons_ct),
            length = sum(length))

DCEregions_Final[order(DCEregions_Final$length,decreasing = T),]$length[1:3]

library(ggplot2)
library(ggbreak)
print(paste("total exon count:", sum(DCEregions_Final$exons_ct)))
main_hist<-ggplot(DCEregions_Final, aes(x = length))+
  geom_histogram(aes(fill = factor(exons_ct)), binwidth = 15)+
  scale_fill_manual(values = c("#e3a9e3", "#c783d6", "#aa72c4", "#452661", "orange"),
                    name = "exons per DCR")+
  scale_x_break(breaks = c(801, 1850),ticklabels = 1870)+  
  xlab("DCR length")+
  ylab("DCR count")+
  theme_classic()+
  theme(legend.position = "inside",panel.background = element_blank())

zoom1<-ggplot(DCEregions_Final, aes(x = length))+
  geom_histogram(aes(fill = factor(exons_ct)), binwidth = 15)+
  scale_fill_manual(values =  c("#e3a9e3", "#c783d6", "#aa72c4", "#452661", "orange"),
                    name = "exons per DCR")+
  coord_cartesian(ylim = c(0,10), xlim = c(0, 650))+
  xlab("DCR length")+
  ylab("DCR count")+
  theme_classic()+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
        )
# stats for paper:
print(paste("min DCR exon length:", min(DCEregions_Final$length), "bp"))
print(paste("max DCR exon length:", max(DCEregions_Final$length), "bp"))
print(paste("mean DCR exon length:", mean(DCEregions_Final$length), "bp"))
print(paste("sd DCR exon length:", sd(DCEregions_Final$length), "bp"))
print(paste("DCR count:", length(unique(DCEregions_Final$DCEindex))))
table(DCEregions_Final$exons_ct)

# automatic output
F1D<-ggdraw(main_hist) + 
draw_plot(zoom1, x = 0.81,y = 0.70, width = 0.16, height = 0.3,)
ggsave("Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/01_Description/DCR_exon_length_distribution.pdf", F1D, width = 6, height = 4)
