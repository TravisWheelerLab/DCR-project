library(tidyr)
library(dplyr)

NMD<-read.csv("output/01_Description/DCR_NMD_analysis.out", header = F, sep = "\t")
names(NMD)<-c("DCR", "outcome", "reason", "frames")
finalDCR<-read.table("output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed")
NMD<-NMD[NMD$DCR %in% finalDCR$V4,]
NMD %>% group_by(DCR) %>%
  summarize(DCEindex = unique(DCR),
            NMDstatus = unique(reason)) -> NMD_per_index


NMD_per_index <- NMD_per_index[NMD_per_index$NMDstatus != 2,]
NMD_per_index$NMDstatus<-factor(NMD_per_index$NMDstatus, 
                                levels = rev(c("last_exon", "coding_long_enough",  "not_C-terminal", "in_NMD_zone")),
                                labels = rev(c("Not enough evidence", "Not enough evidence", "DCR exons not C-Terminal", "NMD" )))

table(NMD_per_index$NMDstatus)
