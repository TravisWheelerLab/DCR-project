#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
### Human DCR ranked by PP and mouse status V2

dcePP<-read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation/meanPP_per_DCE.tsv", h = T)
# length(unique(dcePP[dcePP$DCE_ID %in% unique(humDCE2F_noITR$Index), ]$DCE_ID))
names(dcePP)<-c("DCE_ID", "meanPP")
length(unique(dcePP$DCE_ID))
dce2gene<-read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation/DCEtoGene.txt")
#length(unique(dce2gene$DCE_ID))
names(dce2gene)<-c("DCE_ID", "gene")
dcePPwgene<-merge(dcePP, dce2gene, by = "DCE_ID")
#length(unique(dcePPwgene$DCE_ID)) 1992

dcePPwgene$gene <- toupper(dcePPwgene$gene)

DCEframesCount<-read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation/frames_per_DCR.txt")
names(DCEframesCount)<-c("DCE_ID", "frames")
dcePPwgene<-merge(dcePPwgene, DCEframesCount, by = "DCE_ID")
#length(unique(dcePPwgene$DCE_ID)) 1726

nnDCEmean<-mean(dcePP$meanPP)#mean(as.numeric(PP[PP$type == "nDCE_exon_nDCE_gene",]$mean_PhyloP100Way), na.rm = T)
nnDCEmed<-median(dcePP$meanPP)#median(as.numeric(PP[PP$type == "nDCE_exon_nDCE_gene",]$mean_PhyloP100Way), na.rm = T)
nnDCEsd<-sd(dcePP$meanPP)#sd(as.numeric(PP[PP$type == "nDCE_exon_nDCE_gene",]$mean_PhyloP100Way), na.rm = T)
low_1_nnDCEmean<-quantile(x = as.numeric(dcePP$meanPP), prob = c(0.01, 0.99), na.rm = T)[1]#quantile(x = as.numeric(PP[PP$type == "nDCE_exon_nDCE_gene",]$mean_PhyloP100Way), prob = c(0.01, 0.99), na.rm = T)[1]
high_99_nnDCEmean<-quantile(x = as.numeric(dcePP$meanPP), prob = c(0.01, 0.99), na.rm = T)[2]#quantile(x = as.numeric(PP[PP$type == "nDCE_exon_nDCE_gene",]$mean_PhyloP100Way), prob = c(0.01, 0.99), na.rm = T)[2]

infoQT<-as.data.frame(cbind(c("1st percentile", "DCR median", "99th percentile"),
                            c(low_1_nnDCEmean, nnDCEmed, high_99_nnDCEmean),
                            rep(-100, 3)))
names(infoQT)<-c("label", "y", "x")

mice_satus<-read.csv("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation/DCRmiceStatus.csv", sep = ",")

# exon_to_dce <- read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation/mouse_exon_to_dce_id.txt")
# names(exon_to_dce)<-c("ID", "ORF2exon_ratio", "DCE_ID")
# # report if exon is dual coding in mice:
# exon_to_dce[exon_to_dce$ORF2exon_ratio >= 0.8,] %>% group_by(ID) %>%
#   reframe(dual_coding = ifelse(length(ID) > 1, 1, 0),
#           orf_over_80pc = length(ID),
#           DCE_ID = unique(DCE_ID)) %>% group_by(DCE_ID) %>%
#   summarize(exons_dce_level = sum(dual_coding)/length(DCE_ID),
#             mouse_DC_status = ifelse(sum(dual_coding)/length(DCE_ID) == 1, "DCE", 
#                                      ifelse(sum(dual_coding)/length(DCE_ID) == 0, "non-DCE",
#                                             "partial-DCE"))) -> exonID_to_DCEstatus


dcePPwgene_wID <- merge(dcePPwgene, mice_satus, by = "DCE_ID")
print(length(unique(dcePPwgene_wID$DCE_ID)))
dcePPwgene_wID$mice_status <- factor(dcePPwgene_wID$mice_status, levels = c("DCR", "not DCR", "too short", "not mapped"))

#length(unique(dcePPwgene_wID$DCE_ID))

inset<-ggplot(dcePPwgene_wID)+
  geom_bar(aes(x = reorder(as.character(DCE_ID),meanPP),
               y = meanPP,
               fill = (meanPP > high_99_nnDCEmean | meanPP < low_1_nnDCEmean)),
           stat = "identity", show.legend = F)+
  scale_fill_manual(values = c("grey", "#FFA800"))+
  geom_hline(yintercept = nnDCEmed, lty = 2, col = "grey20",)+
  geom_hline(yintercept = low_1_nnDCEmean, lty = 2, col = "grey")+
  geom_hline(yintercept = high_99_nnDCEmean, lty = 2, col = "grey")+
  #geom_text(data = infoQT, aes(y = as.numeric(y)+0.2, x = 50, label = as.character(label),
  #                             hjust = 0), col = c("grey50", "black", "grey50"))+
  xlab("all DCRs\n")+
  ylab("mean PhyloP")+
  coord_cartesian(ylim = c(-1.5,8))+
  geom_rect(aes(xmin = -35, xmax = length(DCE_ID)+35, ymin = -1.2, ymax = 7.9), fill = NA, col = "black")+
  theme_classic()+
  theme(legend.position = "inside",legend.direction = "horizontal",
        legend.position.inside = c(0.5,0.9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_text(size = 10, margin = margin(0,0,0,0, "cm")),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_blank(),
        plot.background = element_blank()
  )


PPoutliers<-dcePPwgene_wID[dcePPwgene_wID$meanPP > high_99_nnDCEmean | dcePPwgene_wID$meanPP < low_1_nnDCEmean,]
plotW<-length(PPoutliers$DCE_ID)
main_p<-ggplot(PPoutliers)+
  geom_bar(aes(x = reorder(as.character(DCE_ID),meanPP),
               y = meanPP,
               fill = mice_status)
           , stat = "identity")+
  scale_fill_manual(values = c("#FFA800", "orange4", "#FFF0BD"), name = "mouse status:")+
  geom_segment(aes(y = nnDCEmed, yend = nnDCEmed, x = 0, xend = plotW+1), lty = 2)+
  geom_segment(aes(y = low_1_nnDCEmean, yend = low_1_nnDCEmean, x = 0, xend = plotW+1), lty = 2, col = "grey")+
  geom_segment(aes(y = high_99_nnDCEmean, yend = high_99_nnDCEmean, x = 0, xend = plotW+1), lty = 2, col = "grey")+
  geom_text(aes(x = reorder(as.character(DCE_ID), meanPP),
                y = ifelse(meanPP > 0, meanPP + 0.1, meanPP - 0.1),  # Adjust position based on sign
                label = gene
  ),
  size = 3, hjust = ifelse(PPoutliers$meanPP > 0, 0, 1), angle = 90, vjust = 0.5) + # Adjust vjust for placement
  geom_text(data = infoQT, aes(y = as.numeric(y), x = plotW + 1, label = as.character(label)),
            hjust = 0, size = 3.4, col = c("grey50", "black", "grey50"))+
  xlab("DCR region (outliers)                   ")+
  ylab("PhyloP score (100-Ways)")+
  coord_cartesian(ylim = c(-2,9), xlim = c(0,plotW+5))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.position.inside = c(0.5,0.07))

figure2C<-cowplot::ggdraw(main_p) +
  cowplot::draw_plot(inset,x = 0.06 ,y = 0.70, width = 0.25, height = 0.3,)# axis.ticks.y = element_blank(),
ggsave("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation/Figure2C_DCR_PhyloP_outliers.pdf", figure2C, width = 8, height = 6)
