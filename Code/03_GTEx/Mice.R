# mice diff usage
library(reshape2)
library(stringr)
library(ggh4x)
library(scales)
mice_counts<-read.table("Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/data/Mice/Mice.DCRjunctions.counts.tsv", h = T)
mice_meta<-read.table("Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/data/Mice/RNAseq.metadata", h = T)
mice_ct_table<-melt(mice_counts, variable.name = "run", value.name = "count")
mice_ct_table<-merge(mice_ct_table, mice_meta, by = "run")
mice_ct_table<-separate_wider_delim(mice_ct_table, cols = needle, delim = "_", names = c("DCR", "gene", "isoform", "frame_status", "junction"))

miceTMMfactors<-as.data.frame(calcNormFactors(mice_counts[,-1]))
miceTMMfactors$sample<-rownames(miceTMMfactors)
names(miceTMMfactors)<-c("NormFactors", "run")

mice_ct_table<-merge(mice_ct_table, miceTMMfactors, by = "run")
mice_ct_table$pairID<-paste(mice_ct_table$DCR, mice_ct_table$gene, sep = "_")

unique(mice_ct_table$tissue)
mice_cols<-c("#CDB79E", "#EEEE00", "#EED5D2", "#7A67EE", "#DDBB77",
             "#DDBB85", "#FFD39B", "#FFD39B", "#EE82CC","#FFB6C1",
             "lightblue", "#B452CD", "#9ACD32", "#778855", "#8FBC8F",
             "#CDB79E", "#AAD39B", "#A6A6A6")

mice.tissue.colors<-data.frame(tissue = unique(mice_ct_table$tissue), 
                               color = mice_cols)
mice_ct_table<-merge(mice_ct_table, mice.tissue.colors, by = "tissue")
mice_ct_table$frame_status<-factor(mice_ct_table$frame_status, levels = c("can", "alt"), labels = c("canonical", "alternative"))
# function
mice_tissue_plot<-function(pair_id = NULL, data = mice_ct_table){
# data <- mice_ct_table
# pair_id <- "523_MADD"

done <- data[data$pairID == pair_id, ]

# normalize all
done$CPM <- (done$count / (done$lib_size * done$NormFactors) * 1e6)
mean.t.CPM <- done %>%
  group_by(tissue) %>%
  summarize(meanCPM = mean(CPM, na.rm = TRUE), .groups = "drop")
#CPMthreshold <- signif(min(mean.t.CPM$meanCPM), 2)
#done$tissue<-factor(done$tissue, levels = ordered.hcols$tissue)
myCols <- done$color
names(myCols) <- done$tissue
myTiss <- done$tissue
ct_table <- data.frame(color = myCols, tissue = myTiss)
ct_table_uniq <- distinct(ct_table)
labelcol<-ct_table_uniq[order(ct_table_uniq$tissue),]$color
strips <- strip_themed(background_x = elem_list_rect(fill = labelcol))

#done$frame<-sub("frame", "", done$frame)
#done$key<-paste(done$ID, done$frame, sep = "_")
#done<-unique(merge(done, perIsoDCR[,c("frame_status", "key")], by = "key"))
#done$frame_status<-factor(done$frame_status, levels = c("canonical", "alternative")) 
done$junction<-ifelse(done$junction == "left", "upstream", "downstream")
#done$ShortName<-factor(done$ShortName, levels = ordered.hrows$tissue)
p<-ggplot(done, aes(x = factor(frame_status), y = CPM)) +
  #geom_hline(aes(yintercept = CPMthreshold), linetype = 2) +
  geom_jitter(aes(col = frame_status), width = 0.1, size = 1) +
  geom_boxplot(aes(col = frame_status), width = 0.4, fill = "white", alpha = 0.5, outlier.alpha = 0) +
  scale_color_manual(values = c("orange", "purple"))+
  #scale_color_manual(values = c("orange", "purple", "firebrick"))+
  xlab("tissue/frame") +
  ylab("Normalized Counts (CPM-TMM)") +
  ggtitle(paste(
    "DCR ortholog id:", unique(done$DCR),
    "| Gene:", unique(done$gene),
    "| Junction:", unique(done$junction)#,
    #  "| Raw reads threshold:", r_threshold,
    #  "| CPM threshold:", CPMthreshold)) +
  ))+
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    #axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom", legend.direction = "horizontal",
    strip.text.x = element_text(angle = 90)
  )+
  facet_wrap2(~tissue, nrow = 1, strip = strips)
return(p)
}

# simplify human data to match mouse
ad_for_mouse <- associated_data
unique(ad_for_mouse$tissue)
Tmatches <- grep(paste(unique(mice_ct_table$tissue), collapse = "|"), unique(ad_for_mouse$tissue))
unique(ad_for_mouse$tissue)[Tmatches]

madd_human<-plot_function2("523_left")
madd_mice<-mice_tissue_plot("523_MADD")
madd_both<-plot_grid(madd_human, madd_mice, ncol =1)
# works well, brain

human_mknk1<-plot_function2("1287_left")
mice_mknk1<-mice_tissue_plot("1287_MKNK1")
mknk1_both<-plot_grid(human_mknk1, mice_mknk1, ncol =1)
# too specific of some brain tissue in human

human_slc35a2<-plot_function2("638_left")
mice_slc35a2<-mice_tissue_plot("638_SLC35A2")
slc35a2_both<-plot_grid(human_slc35a2, mice_slc35a2, ncol =1)
# too specific of some brain tissue in human

human_mpzl1<-plot_function2("1912_left")
mice_mpzl1<-mice_tissue_plot("1912_MPZL1")
mpzl1_both<-plot_grid(human_mpzl1, mice_mpzl1, ncol =1) 
# note can has 2 isoforms and we can see a bimodal distro

human_ampd2<-plot_function2("907_right")
mice_ampd2<-mice_tissue_plot("907_AMPD2")
ampd2_both<-plot_grid(human_ampd2, mice_ampd2, ncol =1) 
# nope check if both frames are quantified

human_grb7<-plot_function2("1278_left")
mice_grb7<-mice_tissue_plot("1278_GRB7")
# nope check if both frames are quantified

human_atp2b4<-plot_function2("518_left")
mice_atp2b4<-mice_tissue_plot("518_ATP2B4")
atp2b4_both<-plot_grid(human_atp2b4, mice_atp2b4, ncol =1) 
# works, limited but works

human_rgs12<-plot_function2("781_left")
mice_rgs12<-mice_tissue_plot("781_RGS12")
rgs12_both<-plot_grid(human_rgs12, mice_rgs12, ncol =1) 
# limited, but somewhat works

human_pbx1<-plot_function2("1664_right")
mice_pbx1<-mice_tissue_plot("1664_PBX1")
pbx1_both<-plot_grid(human_pbx1, mice_pbx1, ncol =1) 
# problem --> not same gene

human_magi1<-plot_function2("1450_left")
mice_magi1<-mice_tissue_plot("1450_MAGI1")
magi1_both<-plot_grid(human_magi1, mice_magi1, ncol =1) 
# works well (brain)

human_atp2b1<-plot_function2("774_left")
mice_atp2b1<-mice_tissue_plot("774_ATP2B1")
atp2b1_both<-plot_grid(human_atp2b1, mice_atp2b1, ncol =1) 
# works well (brain)

human_mark4<-plot_function2("1630_left")
mice_mark4<-mice_tissue_plot("1630_MARK4")
mark4_both<-plot_grid(human_mark4, mice_mark4, ncol =1) 


# --- 
  
# combining and filtering tissues per species for comparison
hum_tiss <- read.table("Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/data/Mice/human_tissues_correspondance.tsv", h = T)
mou_tiss <- read.table("Google Drive/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/data/Mice/mouse_tissues_correspondance", h = T)

human_rna <- associated_data[associated_data$tissue %in% hum_tiss$tissue & associated_data$pairID %in% DCRdotMatrix_df$pairID,]
human_rna$frame<-sub("frame", "", human_rna$frame)
human_rna$key<-paste(human_rna$ID, human_rna$frame, sep = "_")
human_rna<-unique(merge(human_rna, perIsoDCR[,c("frame_status", "key")], by = "key"))
human_rna<-merge(human_rna, hum_tiss, by = "tissue")
human_rna$tissue <- human_rna$new_tissue
human_rna_simple <- human_rna[,c("sample", "tissue", "ID", "count", "lib_size", "NormFactors", "frame_status", "end")]
names(human_rna_simple) <- c("run", "tissue", "DCR", "count", "lib_size", "NormFactors", "frame_status", "junction")
human_rna_simple$species <- rep("human", length(human_rna_simple$run))
human_rna_colorkey<- unique(human_rna[,c("tissue", "color")])
human_rna_colorkey$color <- c("#aed6ae",
                              "#EEEE00",
                              "#d4bdae",
                              "#d4bdae",
                              "#d2a4de",
                              "#c6cf93",
                              "#dee6b3",
                              "#add3f0",
                              "#bcb5eb",
                              "#FFB6C1",
                              "#bfa695",
                              "#f7c6a6",
                              "#e3a47b",
                              "#A6A6A6",
                              "#EED5D2")

mouse_rna <- mice_ct_table[mice_ct_table$tissue %in% mou_tiss$tissue,]
mouse_rna <- merge(mouse_rna, mou_tiss, by = "tissue")
mouse_rna$tissue <- mouse_rna$new_tissue
mouse_rna_simple <- mouse_rna[,c("run", "tissue", "DCR", "count", "lib_size", "NormFactors", "frame_status", "junction")]
mouse_rna_simple$species <- rep("mouse", length(mouse_rna_simple$run))

humice_rna <- rbind(human_rna_simple, mouse_rna_simple)
humice_rna <- merge(humice_rna, human_rna_colorkey, by = "tissue")
humice_rna$frame_status <- factor(humice_rna$frame_status, levels = c("canonical", "alternative"))

DCR2GENE <- unique(human_rna[,c("ID", "gene")])
names(DCR2GENE) <- c("DCR", "gene")
humice_rna <- merge(humice_rna, DCR2GENE, by = "DCR")
humice_rna$gene <- toupper(humice_rna$gene)


homo_mice_plot<-function(DCR = NULL, data = humice_rna){
  # data <- humice_rna
  # DCR <- "523"
  
  done <- data[data$DCR == DCR, ]
  
  # normalize all
  done$CPM <- (done$count / (done$lib_size * done$NormFactors) * 1e6)

  myCols <- done$color
  names(myCols) <- done$tissue
  myTiss <- done$tissue
  ct_table <- data.frame(color = myCols, tissue = myTiss)
  ct_table_uniq <- distinct(ct_table)
  labelcol<-ct_table_uniq[order(ct_table_uniq$tissue),]$color
  #strips <- strip_themed(background_x = elem_list_rect(fill = labelcol))
  
  done$junction<-ifelse(done$junction == "left", "upstream", "downstream")
  p<-ggplot(done, aes(x = factor(frame_status), y = CPM)) +
    geom_jitter(aes(col = frame_status), width = 0.1, size = 1) +
    geom_boxplot(aes(col = frame_status), width = 0.4, fill = "white", alpha = 0.5, outlier.alpha = 0) +
    scale_color_manual(values = c("orange", "purple"))+
    #scale_color_manual(values = c("orange", "purple", "firebrick"))+
    xlab("tissue/frame") +
    ylab("Normalized Counts (CPM-TMM)") +
    ggtitle(paste(
      "DCR id:", unique(done$DCR),
      "| Gene:", unique(done$gene),
      "| Junction:", unique(done$junction)#,
      #  "| Raw reads threshold:", r_threshold,
      #  "| CPM threshold:", CPMthreshold)) +
    ))+
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      #axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "bottom", legend.direction = "horizontal",
      strip.text.x = element_text(angle = 90)
    )+
    facet_grid(species~tissue, scales = "free_y")
  return(p)
}

library(dplyr)
library(ggplot2)
library(ggh4x)   # for facet strip theming per-panel


homo_mice_plot2 <- function(DCR_id = NULL, data = humice_rna) {
  
  # subset and normalize
  done <- data %>%
    filter(DCR == DCR_id) %>%   # <- critical fix
    mutate(
      CPM = count / (lib_size * NormFactors) * 1e6,
      junction = ifelse(junction == "left", "upstream", "downstream")
    )
  
  # ---------- facet strip color prep ----------
  # unique tissue -> color map (first color wins if duplicated)
  tissue_colors_df <- done %>%
    select(tissue, color) %>%
    filter(!is.na(tissue), !is.na(color)) %>%
    distinct(tissue, .keep_all = TRUE)
  
  # order of tissues in facets (the order they appear in tissue_colors_df)
  # sort tissues alphabetically
  tissue_colors_df <- tissue_colors_df %>% arrange(tissue)
  tissue_levels <- c("Adrenal_Gland", "Brain", "Heart", "Muscle", "Lung",
                     "Spleen", "Stomach", "Intestine_small", "Intestine_large", "Kidney", "Liver",
                     "Ovary", "Uterus", "Testis")#tissue_colors_df$tissue
  
  # apply that order to the plotting data
  done <- done %>%
    mutate(
      tissue = factor(tissue, levels = tissue_levels)
    )
  
  # build strip fill vector in same order
  strip_fill_vec <- tissue_colors_df$color[
    match(tissue_levels, tissue_colors_df$tissue)
  ]
  # optional fallback for missing colors
  strip_fill_vec[is.na(strip_fill_vec)] <- "white"
  
  # ---------- N == 2 labels ----------
  count_df <- done %>%
    group_by(species, tissue, frame_status) %>%
    summarise(
      N = n(),
      y_pos = max(CPM, na.rm = TRUE) + 0.1, #* 1.5,
      .groups = "drop"
    ) %>%
    filter(N == 2) %>%
    mutate(label = paste0("N=", N))
  
  # ---------- plot ----------
  
  p <- ggplot(done, aes(x = factor(frame_status), y = CPM)) +
    geom_jitter(aes(col = frame_status), width = 0.1, size = 1) +
    geom_boxplot(
      aes(col = frame_status),
      width = 0.4, fill = "white", alpha = 0.5, outlier.alpha = 0
    ) +
    geom_text(
      data = count_df,
      aes(x = frame_status, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 3,
      vjust = 0
    ) +
    scale_color_manual(
      name = "DCR frame",
      values = c(
        "alternative" = "purple",
        "canonical"   = "orange"
      )
    ) +
    xlab("Tissue") +
    ylab("Normalized Counts (CPM-TMM)") +
    ggtitle(paste(
      "DCR id:", paste(unique(done$DCR), collapse = ", "),
      "| Gene:", paste(unique(done$gene), collapse = ", "),
      "| Junction:", paste(unique(done$junction), collapse = ", ")
    )) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      legend.direction = "vertical",
      strip.text.x = element_text(angle = 90)
    ) +
    ggh4x::facet_grid2(
      rows = vars(species),
      cols = vars(tissue),
      scales = "free_y",
      #independent = "y",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill  = strip_fill_vec,
          color = "grey20"
        ),
        text_x = ggh4x::elem_list_text(
          colour = rep("black", length(strip_fill_vec)),
          face   = "bold",
          size   = 8
        )
      )
    )
  
  return(p)
}






MADD_MICE <- homo_mice_plot2(DCR = 523) # MADD # pretty good | uterus different
ATP2B4_MICE <- homo_mice_plot2(DCR = 518) # pretty good
MADD_MICE_ALIGN <- cowplot::ggdraw(cowplot::align_plots(MADD_MICE, ATP2B4_MICE, align = "v")[[1]])
DCRbrowser(pairID = "518_left", DCRhighlight = "yellow")

cowplot::plot_grid(MADD_MICE_ALIGN, ATP2B4_MICE, nrow = 2)


homo_mice_plot2(DCR = 1630) # pretty good

homo_mice_plot2(DCR = 774) # atp2b1 # good but not super clean
homo_mice_plot2(DCR = 1450) # MAGI1 # good but brain not so strong in human
homo_mice_plot2(DCR = 1912) # midl
homo_mice_plot2(DCR = 1664) # mild and different in mous

homo_mice_plot2(DCR = 781) # not really
homo_mice_plot2(DCR = 1287) # nope
homo_mice_plot2(DCR = 638) # nope because within brain specific human

homo_mice_plot2(DCR = 907) # alt not quantified in mice
homo_mice_plot2(DCR = 1278) # alt not quantified in mice


mice_ct_table %>% group_by(tissue) %>% summarize(N = length(unique(run)))
length(unique(mice_ct_table$run))
humice_rna %>% filter(species == "mouse") %>% group_by(run) %>% summarize(
                                              tissue = unique(tissue),
                                              lib_size = unique(lib_size)) -> mouse_samples
write.table(mouse_samples, "/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2025/11-12_DCR_Paper/Sup/Table_S2.tsv", row.names = F, quote = F)
