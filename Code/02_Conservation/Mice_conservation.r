#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(magick)
library(pdftools)
library(cowplot)
library(dplyr) 

option_list = list(
  make_option(c("-d", "--dce"), type="character", default=NULL, help="Path to DCE file", metavar="character"),
  make_option(c("-n", "--nondce"), type="character", default=NULL, help="Path to nonDCE file", metavar="character"),
  make_option(c("-e", "--expdce"), type="character", default=NULL, help="Path to expDCE file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$dce) | is.null(opt$nondce) | is.null(opt$expdce)) {
  print_help(opt_parser)
  stop("All three input files must be supplied (--dce, --nondce, --expdce)", call.=FALSE)
}

# data
DCE <- read.table(opt$dce, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
nonDCE <- read.table(opt$nondce, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
expDCE <- read.table(opt$expdce, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# DCE <- read.table("~/Google Drive/My Drive/Wheeler_Lab/2026/01_22_DCR-Mouse-cons-redux/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# nonDCE <- read.table("~/Google Drive/My Drive/Wheeler_Lab/2026/01_22_DCR-Mouse-cons-redux/fake_DCE_intervals.fa_translated.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# expDCE <- read.table("~/Google Drive/My Drive/Wheeler_Lab/2026/01_22_DCR-Mouse-cons-redux/fake_mm39_exons.fasta_translated.fa_linearized_renamedWithCoord.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# allD<-read.table("~/Google Drive/My Drive/Wheeler_Lab/2026/01_22_DCR-Mouse-cons-redux/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
header <- c("ORF", "chr", "start", "end", "strand", "ExLength", "ORFstart", "ORFend", "ORFlen_nt", "ORF2exon_ratio", "ORF_len_AA", "frame", "AA_seq")
names(DCE)=header
names(nonDCE)=header
names(expDCE)=header
names(allD)=header

# add key to each
DCE$ID<-paste("DCE", DCE$chr, DCE$start, DCE$end, sep=":")
DCE$type<-rep("DCE",length(DCE$ORF))

nonDCE$ID<-paste("nonDCE", nonDCE$chr,nonDCE$start,nonDCE$end, sep=":")
nonDCE$type<-rep("nonDCE",length(nonDCE$ORF))
expDCE$ID<-paste("expDCE", expDCE$chr,expDCE$start,expDCE$end, sep=":")
expDCE$type<-rep("expDCE",length(expDCE$ORF))

#connect DCR id
exon_to_dce <- read.table("~/Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/02_Conservation/mouse_exon_to_dce_id.txt")
names(exon_to_dce)<-c("ID", "ORF2exon_ratio", "DCE_ID")
exon_to_dce$exID<-gsub("DCE:", "", exon_to_dce$ID)
DCE$exID<-paste(DCE$chr, DCE$start, DCE$end, sep = ":")

e2dce<-exon_to_dce[,c("exID", "DCE_ID")]
DCE<-merge(DCE, e2dce, by = "exID")

# filter the DCE to only those with human sub-regions >= 45nt
DCE <- DCE[!DCE$DCE_ID %in% unique(DCE[DCE$ExLength < 21,]$DCE_ID),] 

Exons<-rbind(unique(DCE[,-c(1,17)]), nonDCE, expDCE)

Exons$ORFstart<-Exons$ORFstart-1
Exons$ORFlen_nt<-Exons$ORFend-Exons$ORFstart
Exons$ORF2exon_ratio<-Exons$ORFlen_nt/Exons$ExLength

# now we look for each orf if it is dual coding by checking if ORFstart-3 <= 0 & ORFend-ExLength <= 3
Exons$isDualCoding<-ifelse(Exons$ORFstart-3 <= 0 & Exons$ExLength-Exons$ORFend <= 3, TRUE, FALSE)

# new filter: not 80% of region but all region +/- 3bp on either sides
Exons %>% group_by(ID) %>% filter(isDualCoding == TRUE & type == "DCE") %>%
  summarize(frames = length(unique(frame)),
  ) %>% filter(frames > 1) -> mice_DCR_per_ortholog

# Total mapped DCR exons
print(paste("Total mapped DCR exons:", length(unique(DCE$exID))))
# Total mapped DCR regions to exons
print(paste("Total mapped DCR regions to mouse exons:",length(unique(DCE$DCE_ID))))
# DCR exons passing filter
print(paste("DCRexons passing filter:", length(unique(mice_DCR_per_ortholog$ID))))

############
# function #
############
#threshold=0.8 # filter ORFs >= N% exon length
plotConsDCE<-function(ORFtable = Exons, threshold = 0.8, xmin = 0, xmax = 500, binlen = 1){
thr<-threshold # irrelevant now
# new filter: not 80% of region but all region +/- 6bp on either sides
Exons_f<-Exons[Exons$isDualCoding == TRUE,]
# count how many exons have 2+ orfs per length bin
OpE<-as.data.frame(tapply(Exons_f$ID, Exons_f$ID, length))
OpE_df<-cbind(rownames(OpE), OpE)
names(OpE_df)<-c("ID", "ORFs")
OpE_df_f_s_Len<-unique(Exons_f[,c("ID", "ExLength")])
OpE_df<-merge(OpE_df, OpE_df_f_s_Len, by = "ID")
OpE_df$ORFs_bin<-ifelse(OpE_df$ORFs > 1, 1, 0) # make count binary (1 = 0, 2+ = 1)
OpE_df<-merge(OpE_df, unique(Exons_f[,c("ID", "type")]), by = "ID")
# make bins
#bins<-c(1:100, seq(122,max(DCE$ExLength),round(sqrt(length(unique(Exons_f$ID)))/2)))
breaks<-c() # initialize bin breaks 
#pos<-0
i <- 1
while(i < max(OpE_df$ExLength)){
  if(i <= 100){
    binwidth = binlen
    breaks[i] <- i + binwidth
    #pos = pos + binwidth
    i = i + binwidth
    #print(pos)
  } else {
    binwidth = ceiling(sqrt(i)/2)
    breaks[i] = i + binwidth
    #pos = pos + binwidth
    i = i + binwidth
    #print(pos)
  }
}
breaks<-breaks[!is.na(breaks)]
# a function to get midpoints (https://stackoverflow.com/questions/22312207/how-to-assign-cut-range-midpoints-in-r)
get_midpoint <- function(cut_label) {
  mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ","))))
}
# slice for non DCE
ND_OpE_df<-OpE_df[OpE_df$type == "nonDCE",]
ND_OpE_df$bin_label<-cut(ND_OpE_df$ExLength, breaks)
ND_OpE_df$midpoints<-sapply(cut(ND_OpE_df$ExLength, breaks), get_midpoint)
#ND_OpE_df$bin_number<-as.numeric(cut(ND_OpE_df$ExLength, bins))
ND_OpE_df$bin_number<-as.numeric(cut(ND_OpE_df$ExLength, breaks))
ND_OpE_df %>% group_by(midpoints) %>% summarize(bin_count = length(midpoints), type = "nonDCE") -> ND_bin_counts
ND_ratios<-as.data.frame(tapply(ND_OpE_df$ORFs_bin, ND_OpE_df$midpoints, function(x){sum(x)/length(x)}))
ND_ratios$bins<-rownames(ND_ratios)
names(ND_ratios)<-c("prop", "bins")
ND_ratios$type<-rep("nonDCE", length(ND_ratios$bins))
# slice for DCE
D_OpE_df<-OpE_df[OpE_df$type == "DCE",]
D_OpE_df$bin_label<-cut(D_OpE_df$ExLength, breaks)
D_OpE_df$midpoints<-sapply(cut(D_OpE_df$ExLength, breaks), get_midpoint)
#D_OpE_df$bin_number<-as.numeric(cut(D_OpE_df$ExLength, bins))
D_OpE_df$bin_number<-as.numeric(cut(D_OpE_df$ExLength, breaks))
D_OpE_df %>% group_by(midpoints) %>% summarize(bin_count = length(midpoints), type = "DCE") -> D_bin_counts
D_ratios<-as.data.frame(tapply(D_OpE_df$ORFs_bin, D_OpE_df$midpoints, function(x){sum(x)/length(x)}))
D_ratios$bins<-rownames(D_ratios)
names(D_ratios)<-c("prop", "bins")
D_ratios$type<-rep("DCE", length(D_ratios$bins))
# slice for expDCE
E_OpE_df<-OpE_df[OpE_df$type == "expDCE",]
E_OpE_df$bin_label<-cut(E_OpE_df$ExLength, breaks)
E_OpE_df$midpoints<-sapply(cut(E_OpE_df$ExLength, breaks), get_midpoint)
#E_OpE_df$bin_number<-as.numeric(cut(E_OpE_df$ExLength, bins))
E_OpE_df$bin_number<-as.numeric(cut(E_OpE_df$ExLength, breaks))
E_OpE_df %>% group_by(midpoints) %>% summarize(bin_count = length(midpoints), type = "expDCE") -> E_bin_counts
E_ratios<-as.data.frame(tapply(E_OpE_df$ORFs_bin, E_OpE_df$midpoints, function(x){sum(x)/length(x)}))
E_ratios$bins<-rownames(E_ratios)
names(E_ratios)<-c("prop", "bins")
E_ratios$type<-rep("expDCE", length(E_ratios$bins))

# expected distribution
probab=c()
for(m in 0:5000){
  #n = as.integer(m/3)-2*(m*0.1)
  #probab[m] = 1 - (1-(61/64)^n)^2
  probab[m] = 1 - (1-(61/64)^(m/3))^2
}
probs=as.data.frame(cbind(as.vector(probab), c(1:5000)))
names(probs)<-c("prop", "ExLength")
probs$bins<-cut(probs$ExLength, breaks)
probs$midpoints<-sapply(cut(probs$ExLength, breaks), get_midpoint)
expected<-as.data.frame(tapply(probs$prop, probs$midpoints, mean))
expected$bins<-rownames(expected)
names(expected)<-c("prop", "bins")
expected$type<-rep("expected_basic", length(expected$bins))

# combine
DND_DCE<-rbind(ND_ratios, D_ratios, E_ratios)#, expected)
#DND_DCE$bins<-factor(DND_DCE$bins, levels = unique(DND_DCE$bins))
DND_DCE$type<-factor(DND_DCE$type, levels = c("DCE", "nonDCE", "expDCE"), labels = c("human DCE", "human non DCE [random region]", "human non DCE [shuffled codons]"))

all_bin_counts <- rbind(ND_bin_counts, D_bin_counts, E_bin_counts)
return(list(main = DND_DCE, bin_counts = all_bin_counts))
}

############################# function end #####################################

plotConsDCE(binlen = 3)->DND_DCE
names(DND_DCE$bin_counts) <- c("bins", "bin_count", "type")# DND_DCE$bin_counts %>% group_by(midpoints) %>%
#  2 summarize(bins = unique(bins),
#             bin_count = length(bins)) -> bin_counts
DND_DCE$bin_counts$type<-factor(DND_DCE$bin_counts$type, levels = c("DCE", "nonDCE", "expDCE"), labels = c("human DCE", "human non DCE [random region]", "human non DCE [shuffled codons]"))
DND_DCE$main<-merge(DND_DCE$main, DND_DCE$bin_count, by = c("bins", "type"))
# leggg<-draw_image(magick::image_read_pdf("../../../2024/09-24_DCE-paper/Figures/Figure_2_subfigures/Conservation_mice_legend.pdf", density = 300))


print(paste("Total mapped DCR >21nt exons:", length(unique(DCE[DCE$ExLength >= 21,]$exID))))
# Total mapped DCR regions to exons
print(paste("Total mapped DCR regions >21nt to mouse exons:",length(unique(DCE[DCE$ExLength >= 21,]$DCE_ID))))
# DCR exons filter
print(paste("DCRexons passing filter:", length(unique(mice_DCR_per_ortholog$ID))))

fDCE<-unique(DCE[DCE$ExLength >= 21, c("ID", "DCE_ID")])
fDCE$MiceDualExon<-ifelse(fDCE$ID %in% mice_DCR_per_ortholog$ID, "yes", "no")
fDCE %>% group_by(DCE_ID) %>% summarize(
  MiceDualExon_list = paste0(unique(MiceDualExon), collapse = ","),
  mice_status = ifelse(MiceDualExon_list == "yes", "DCR", 
                       ifelse(MiceDualExon_list == "no", "not_DCR", "partial"))
    ) -> human_to_mice_DCR_status_21plus


allD$exID<-paste(allD$chr, allD$start, allD$end, sep = ":")
allD<-merge(allD, e2dce, by = "exID")
too_short<-sort(unique(allD[!allD$DCE_ID %in% unique(human_to_mice_DCR_status_21plus$DCE_ID),]$DCE_ID))
#final list of DCR
DCR2keep<-read.table("Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb.bed")
DCRlist2keep<-unique(DCR2keep$V1)
not_mapped<-DCRlist2keep[!DCRlist2keep %in% unique(c(too_short, human_to_mice_DCR_status_21plus$DCE_ID))]

rbind(human_to_mice_DCR_status_21plus[,c("DCE_ID", "mice_status")],
data.frame(DCE_ID = too_short, 
           mice_status = rep("too_short", length(too_short))
                             ),
data.frame(DCE_ID = not_mapped,
           mice_status = rep("not_mapped", length(not_mapped)))
) -> DCRmiceStatus

table(DCRmiceStatus$mice_status)
sum(table(DCRmiceStatus$mice_status))

##############
# FIGURE 2 A #
##############

# simplify, get rid of "partial" (turn it into not_DCR)
DCRmiceStatus[DCRmiceStatus$mice_status == "partial",]$mice_status <- "not_DCR"
DCRmiceStatus$group<-rep(1, length(DCRmiceStatus$DCE_ID))
DCRmiceStatus$mice_status <- factor(DCRmiceStatus$mice_status, levels = c("DCR", "not_DCR", "too_short", "not_mapped"),
                                    labels = c("DCR", "not DCR", "too short", "not mapped"))
ggplot(DCRmiceStatus)+
  geom_bar(aes(x = factor(mice_status), fill = factor(mice_status)), show.legend = F)+
  geom_text(aes(x = factor(mice_status), label = factor(mice_status)), 
            stat = "count", 
            vjust = 0.5, 
            angle = 90, 
            hjust = -0.15,
            size = 3)+
  scale_fill_manual(name = "mouse status", values = c("#FFA800","#FFF0BD", "orange4", "grey20"))+
  scale_y_continuous(breaks = c(0, 250, 500, 750, 1000, 1250, 1500))+
  coord_cartesian(ylim = c(0,1200))+
  xlab("mouse status")+
  ylab("human DCR")+
  theme_classic()+
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

write.table(DCRmiceStatus, "~/Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/02_Conservation/mouse_exon_to_dce_id.txt", quote = F, row.names = F)

##############
# FIGURE 2 B #
##############

mouse_cons<-ggplot(DND_DCE$main)+
  geom_point(aes(x = as.numeric(bins), y = prop, col = type))+ #, shape = bin_count < 2))+
  scale_color_manual(values = c("#FFA800", "gray80", "grey60"), guide = "none")+ #, "royalblue"))+
  #scale_shape_manual(values = c(19, 21), name = "single observation")+
  geom_text(
    data = subset(DND_DCE$main, bin_count < 10), 
    aes(x = as.numeric(bins), y = prop, label = bin_count), 
    vjust = -0.6, size = 3, colour = "#FFA800"
  ) +
  #geom_bar(aes(x = as.numeric(bins), y = log10(bin_count)), stat = "identity")+
  coord_cartesian(xlim = c(21,500))+
  scale_x_continuous(breaks = c(21, 100, 200, 300, 400, 500))+
  geom_vline(xintercept = 21, col = "#FFA800", lty = 2)+
  xlab("Exon region length (bp)")+
  ylab("Proportion of regions with >1 ORFs")+
  #coord_cartesian(xlim = c(0,50))+
  theme_classic()+
  theme(legend.position = "inside",
        legend.position.inside = c(.85,0.65), 
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10))

mouse_hist<-ggplot(DND_DCE$main)+
  geom_bar(aes(x = as.numeric(bins), y = bin_count, fill = type), stat = "identity")+
  scale_fill_manual(values = c("#FFA800", "gray80", "grey60"), guide = "none")+ #, "royalblue"))+
  coord_cartesian(xlim = c(21,500))+
  xlab("Exon region length (bp)")+
  ylab("count")+
  #scale_y_log10()+
  #coord_cartesian(xlim = c(0,50))+
  theme_classic()+
  theme(legend.position = "inside",
        legend.position.inside = c(.85,0.65), 
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

full_mouse_hist <- plot_grid(mouse_hist, mouse_cons, ncol = 1, rel_heights = c(0.15, 0.85), align = "v")
#DND_DCE$main[DND_DCE$main$bin_count < 2,]
# figure2A<-plot_grid(ggdraw() + leggg, full_mouse_hist,
#           rel_widths = c(2/6,3/6))
# head(DND_DCE[DND_DCE$type == "Expected shuffled codons" & DND_DCE$,])
#ggsave("Google Drive/My Drive/Wheeler_Lab/2024/07-02-DCE_Conservation/Figures/Figure_2_subfigures/Conservation_mice.pdf", figure2A, width = 6, height = 4)
ggsave("output/02_Conservation/Conservation_mice.pdf", full_mouse_hist , width = 6, height = 4)
