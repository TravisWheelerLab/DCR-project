#!/usr/bin/env Rscript
library(ggpubr)
library(dplyr)
library(ggplot2)

### PhyloP

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
       stop("Usage: Rscript PhyloP_boxplots.r <PP_file>")
}
PP_file <- args[1] # file with the phyloP scores
PP_title <- args[2] # title to the plot (use quotes)
PP_way <- args[3] # what is the nth wah (7, 20, 100)

PP <- read.table(PP_file)

###### I think I want to bin sizes, and show phyloP boxplot per size bin
#PP<-read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2024/07-02-DCE_Conservation/PPworkflow_053025/DCEvsNonDCE.meanPhyloP100.out")
names(PP)<-c("exon", "mean_PhyloP100Way", "type", "len")
# remove the DCE lines to parse
PP_noDCR<-PP[PP$type != "DCE_exon",]
PP_DCR<-PP[PP$type == "DCE_exon",]
PP_DCR$DCR_ID<-as.numeric(sapply(strsplit(PP_DCR$exon, split = "_"), function(x) x[[3]]))
PP_DCR_keep<-PP_DCR #[PP_DCR$DCR_ID %in% unique(humDCE2F_noITR$Index),]
PP<-rbind(PP_noDCR, PP_DCR[,-5])
#as.nudo.call('rbind', strsplit(PP$exon, split = "_"))[,3]

PP_nDCE_to_sample<-PP[PP$type == "nDCE_exon_nDCE_gene",]
#set.seed(42)
#nDCE_100k<-PP_nDCE_to_sample[sample(1:length(PP_nDCE_to_sample$exon), 100000),]

#PP<-rbind(PP[PP$type == "DCE_exon",], PP[PP$type == "nDCE_exon_DCE_gene",], nDCE_100k)
# restricted to DCE size range
lowDCE<-min(PP[PP$type == "DCE_exon",]$len)+10
maxDCE<-max(PP[PP$type == "DCE_exon",]$len)
PPsize<-PP[PP$len >= lowDCE & PP$len <= maxDCE,]
PPsize$bin <- ifelse(PPsize$len < 50, "20-49bp",
                     ifelse(PPsize$len < 100, "50-99bp",
                            ifelse(PPsize$len < 250, "100-249bp",
                                   ifelse(PPsize$len < 500, "250-499bp", ">=500bp")
                                   #paste("[500,", max(PPsize), sep = " ")
                                   )
                            )
                     )

                                               

comparisons <- list(c("DCE_exon", "nDCE_exon_DCE_gene"), c("nDCE_exon_DCE_gene", "nDCE_exon_nDCE_gene"),
                    c("DCE_exon", "nDCE_exon_nDCE_gene"))
PPsize$bin <- factor(PPsize$bin, levels = c("20-49bp", "50-99bp", "100-249bp","250-499bp", ">=500bp"))
# PPsize$type <- factor(PPsize$type, levels = c("DCE_exon", "nDCE_exon_DCE_gene", "nDCE_exon_nDCE_gene"),
#                       labels = c("DCE exon (D)", "non-DCE exon \nfrom DCE gene (D*)", "non-DCE exon \nfrom non-DCE genes (ND)"))                                           
PPsize_summary <- PPsize %>%
  group_by(type, bin) %>%
  summarise(count = n())
Figure2B<-ggplot(data = transform(PPsize,
                        bin=factor(bin, levels = c("20-49bp", "50-99bp", "100-249bp",
                                              "250-499bp", ">=500bp")
                                   )
                        ),
       aes(x = as.factor(type), y = as.numeric(mean_PhyloP100Way))
       )+
  geom_jitter(data = PPsize[PPsize$type == "DCE_exon",], 
              aes(x = as.factor(type), y = as.numeric(mean_PhyloP100Way)),
              col = "#FFA800", alpha = 1, width = 0.2)+
  geom_boxplot(aes(fill = as.factor(type)), 
               width = 0.7, outlier.shape = NA, alpha = 0.5, notch = 
                 F) +
  scale_fill_manual(values = c("#FFA800", "#FFF0BD" , "grey80"), name = "exon category") +
  geom_text(data = PPsize_summary, aes(x = type, y = -3, label = paste("n =\n", count)), 
            vjust = 0, size = 3, color = "black")+
  scale_x_discrete(labels=c('D', 'D*', 'ND'))+
  xlab("exon category")+
  ylab(PP_title)+
  facet_grid(.~factor(bin))+
  theme_classic()+
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif",
                     tip.length = 0.01) +  # Add pairwise t-test
  theme(legend.position = "bottom")

ggsave(paste("output/02_Conservation/Figure2B_PhyloP_boxplot-", PP_way, "-ways.pdf", sep = ""), Figure2B, width = 8, height = 4)