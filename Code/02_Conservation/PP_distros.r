#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(ggside)
library(ggpubr)
library(rstatix)

# input
PP100=read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation/DCEvsNonDCE.meanPhyloP100.out")
headerPP<-c("exon", "mean_PhyloP", "type", "len")
names(PP100)<-headerPP
# restrict to DCR exons length values
PPminL<-min(PP100[PP100$type == "DCE_exon",]$len)
PPmaxL<-max(PP100[PP100$type == "DCE_exon",]$len)
PP100L<-PP100[PP100$len >= PPminL & PP100$len <= PPmaxL,]

# # slice the data
# PP100L$bin<-as.numeric(as.character(cut(PP100L$len, breaks = seq(from = 20, to = PPmaxL, by = 3),
#                                         labels = seq(from = 21.5, to = PPmaxL-1.5, by = 3))))
# #as.numeric(as.character(cut(PP100L$len, breaks = seq(from = 20, to = PPmaxL, by = 3),
# #                            labels = seq(from = 21.5, to = PPmaxL-1.5, by = 3))))
# # resampling to match DCR distribution
# table(PP100L$type)
# PP100L[PP100L$type == "DCE_exon",] %>% group_by(bin) %>%
#   summarize(DCRcount = length(bin)) -> count_per_bin
# # initialize sampling dfs
# s_nDnDG_df <- data.frame(exon = character(),
#                  mean_PhyloP=character(), 
#                  type=character(),
#                  len=numeric(),
#                  bin=numeric()) 
# s_nDDG_df <- data.frame(exon = character(),
#                          mean_PhyloP=character(), 
#                          type=character(),
#                          len=numeric(),
#                          bin=numeric()) 
# for(i in 1:length(count_per_bin$bin)){
# bintodo <- count_per_bin[i,1]
# sample_size <- as.numeric(count_per_bin[i,2])
# 
# # sample nDCE_exon_nDCE_gene
# sliced_nDnDG_df<-PP100L %>% filter(type == "nDCE_exon_nDCE_gene",
#                                    bin == as.numeric(bintodo))
# nDnDG_L <- length(sliced_nDnDG_df$exon)
# sSize <- min(nDnDG_L, sample_size)
# sampled_nDnDG<-sliced_nDnDG_df[sample(1:nDnDG_L, as.numeric(sSize), replace = F),]
# s_nDnDG_df<-rbind(s_nDnDG_df, sampled_nDnDG)
# 
# # sample nDCE_exon_DCE_gene
# sliced_nDDG_df<-PP100L %>% filter(type == "nDCE_exon_DCE_gene",
#                                    bin == as.numeric(bintodo))
# nDDG_L <- length(sliced_nDDG_df$exon)
# sSize <- min(nDDG_L, sample_size)
# sampled_nDDG<-sliced_nDDG_df[sample(1:nDDG_L, as.numeric(sSize), replace = F),]
# s_nDDG_df<-rbind(s_nDDG_df, sampled_nDDG)
# }
# 
# PP100Ls <- rbind(PP100L[PP100L$type == "DCE_exon",], s_nDnDG_df, s_nDDG_df)
# #table(is.na(PP100Ls$mean_PhyloP))
# 
# 
# PP100Ls <- PP100Ls %>%
#   mutate(
#     type = factor(type, levels = c("DCE_exon","nDCE_exon_DCE_gene","nDCE_exon_nDCE_gene")),
#     mean_PhyloP = as.numeric(as.character(mean_PhyloP))
#   ) %>%
#   filter(!is.na(mean_PhyloP),
#          !is.na(bin))
# 
# PPmain<-ggplot(PP100Ls[order(PP100Ls$type, decreasing = T),])+
#   geom_point(aes(x = len, y = mean_PhyloP, col = type), 
#              alpha = 0.3)+
#   scale_color_manual(values = c("#FFA800", "#FFF0BD" , "grey80"))+
#   xlab("Exon length")+
#   ylab("PhyloP score (100-Ways)")+
#   coord_cartesian(ylim = c(-2.5,10))+
#   scale_x_log10()+
#   theme_classic()+
#   theme(legend.position = "inside",
#         legend.position.inside = c(0.75, 0.85))
# 
# # Add p-values for pairwise comparisons
# # Define the comparisons to be made
# my_comparisons <- list(
#   c("DCE_exon", "nDCE_exon_DCE_gene"),
#   c("DCE_exon", "nDCE_exon_nDCE_gene"),
#   c("nDCE_exon_DCE_gene", "nDCE_exon_nDCE_gene")
# )
# 
# PP100Ls %>% t_test(mean_PhyloP ~ type) -> stat.test
# stat.test <- stat.test %>% add_xy_position(x = "type")
# 
# PPboxplot<-ggplot(PP100Ls, aes(x = type, y = mean_PhyloP))+
#   geom_boxplot(aes(col = type, fill = type), alpha = 0.3, notch = T, show.legend = F)+
#   scale_color_manual(values = c("#FFA800", "#FFF0BD" , "grey80"))+
#   scale_fill_manual(values = c("#FFA800", "#FFF0BD" , "grey80"))+
#   stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, y.position = c(8.5, 9, 9.5))+
#   coord_cartesian(ylim = c(-2.5,10))+
#   theme_classic()+
#   theme(axis.title.y = element_blank(),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.x = element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank())
# 
# PPdensity<-ggplot(PP100Ls[order(PP100Ls$type, decreasing = T),])+
#   geom_density(aes(x = len, col = type, fill = type), 
#              alpha = 0.5, show.legend = F)+
#   scale_color_manual(values = c("#FFA800", "#FFF0BD" , "grey80"))+
#   scale_fill_manual(values = c("#FFA800", "#FFF0BD" , "grey80"))+
#   scale_x_log10()+
#   theme_classic()+
#   theme(axis.title.x = element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank())
# 
# cowplot::plot_grid(PPdensity, NULL, PPmain, PPboxplot,
#                    ncol = 2, axis = "b", align = "h", rel_heights = c(0.2, 0.8),
#                    rel_widths = c(0.8, 0.2))


### function version

plot_len_conservation <- function(
    dataset,
    ways,
    len_col   = "len",
    score_col = "mean_PhyloP",
    type_col  = "type",
    exon_col  = "exon",
    # categories and display
    type_levels = c("DCE_exon","nDCE_exon_DCE_gene","nDCE_exon_nDCE_gene"),
    palette      = c("#FFA800", "#FFF0BD", "grey80"),
    # binning
    breaks_from  = 20,
    bin_width    = 3,
    max_len      = NULL,   # if NULL, computed from data
    # stats + reproducibility
    seed         = 42
){
  stopifnot(all(type_levels %in% unique(dataset[[type_col]])))
  if (is.null(max_len)) {
    max_len <- ceiling(max(dataset[[len_col]], na.rm = TRUE) / bin_width) * bin_width
  }
  
  suppressPackageStartupMessages({
    require(dplyr); require(ggplot2); require(ggpubr)
    require(rstatix); require(cowplot)
  })
  
  set.seed(seed)
  
  # ---------- Bin lengths ----------
  breaks  <- seq(from = breaks_from, to = max_len, by = bin_width)
  labels  <- seq(from = breaks_from + bin_width/2, to = max_len - bin_width/2, by = bin_width)
  df <- dataset %>%
    mutate(
      bin = as.numeric(as.character(cut(.data[[len_col]], breaks = breaks, labels = labels)))
    )
  
  # ---------- Count per bin in the DCE_exon stratum ----------
  count_per_bin <- df %>%
    filter(.data[[type_col]] == type_levels[1]) %>%
    group_by(bin) %>%
    summarize(DCRcount = dplyr::n(), .groups = "drop")
  
  # helpers
  sample_like <- function(.df, .type, .bin, .n){
    slice_pool <- .df %>% filter(.data[[type_col]] == .type, bin == .bin)
    pool_n <- nrow(slice_pool)
    if (pool_n == 0 || .n == 0) return(slice_pool[0,])
    sSize <- min(pool_n, .n)
    slice_pool[sample.int(pool_n, sSize, replace = FALSE), , drop = FALSE]
  }
  
  # ---------- Length-matched resampling for the two comparison groups ----------
  sampled_nDnDG <- bind_rows(lapply(seq_len(nrow(count_per_bin)), function(i){
    bintodo     <- count_per_bin$bin[i]
    sample_size <- count_per_bin$DCRcount[i]
    sample_like(df, type_levels[3], bintodo, sample_size)
  }))
  
  sampled_nDDG <- bind_rows(lapply(seq_len(nrow(count_per_bin)), function(i){
    bintodo     <- count_per_bin$bin[i]
    sample_size <- count_per_bin$DCRcount[i]
    sample_like(df, type_levels[2], bintodo, sample_size)
  }))
  
  # ---------- Assemble analysis set ----------
  df_analysis <- bind_rows(
    df %>% filter(.data[[type_col]] == type_levels[1]),
    sampled_nDnDG,
    sampled_nDDG
  ) %>%
    mutate(
      !!type_col  := factor(.data[[type_col]], levels = type_levels),
      !!score_col := as.numeric(as.character(.data[[score_col]]))
    ) %>%
    filter(!is.na(.data[[score_col]]), !is.na(bin))
  
  # ---------- Stats ----------
  stat.test <- df_analysis %>%
    rstatix::t_test(as.formula(paste(score_col, "~", type_col))) %>%
    add_xy_position(x = type_col)
  
  # ---------- Plots ----------
  minCoord<-min(df_analysis$mean_PhyloP) - 0.5
  maxCoord<-max(df_analysis$mean_PhyloP) + 0.5
  
  PPmain <- ggplot(df_analysis[order(df_analysis[[type_col]], decreasing = TRUE), ]) +
    geom_point(aes(x = .data[[len_col]], y = .data[[score_col]], col = .data[[type_col]]),
               alpha = 0.3) +
    scale_color_manual(values = palette) +
    xlab("Exon length") +
    ylab(paste("PhyloP score (", ways,"-Ways)", sep = "")) +
    coord_cartesian(ylim = c(minCoord, maxCoord)) +
    scale_x_log10() +
    theme_classic() +
    theme(legend.position = "inside", legend.position.inside = c(0.75, 0.85))
  
  PPboxplot <- ggplot(df_analysis, aes(x = .data[[type_col]], y = .data[[score_col]])) +
    geom_boxplot(aes(col = .data[[type_col]], fill = .data[[type_col]]),
                 alpha = 0.3, notch = TRUE, show.legend = FALSE, outlier.alpha = 0) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    ggpubr::stat_pvalue_manual(stat.test, label = "p.adj.signif",
                               tip.length = 0.01, y.position = c(1.1 * maxCoord,
                                                                 1.2 * maxCoord ,
                                                                 maxCoord)) +
    coord_cartesian(ylim = c(minCoord, maxCoord)) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  PPdensity <- ggplot(df_analysis[order(df_analysis[[type_col]], decreasing = TRUE), ]) +
    geom_density(aes(x = .data[[len_col]], col = .data[[type_col]], fill = .data[[type_col]]),
                 alpha = 0, show.legend = FALSE) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    scale_x_log10() +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  combined <- cowplot::plot_grid(
    PPdensity, NULL, PPmain, PPboxplot,
    ncol = 2, axis = "b", align = "hv",
    rel_heights = c(0.3, 0.7),
    rel_widths  = c(0.8, 0.2)
  )
  
  list(
    data_processed = df_analysis,
    stats          = stat.test,
    plots          = list(density = PPdensity, main = PPmain, box = PPboxplot, combined = combined)
  )
}


Plot100W<-plot_len_conservation(PP100L, 100)












### see if we keep that 
# 
# PP100Ls %>% group_by(bin, type) %>% summarize(
#   meanPP = mean(na.omit(as.numeric(mean_PhyloP))),
#   binSize = length(len),
# ) -> PP100Lstat
# 
# 
# # order factor levels
# #typelevels<-
# PP100Lstat$type<-factor(PP100Lstat$type, levels = rev(c("nDCE_exon_nDCE_gene", "nDCE_exon_DCE_gene", "DCE_exon")),
#                         labels = c("DCR exon", "non-DCR exon, DCR gene", "non-DCR exon, non-DCR gene"))
# # plot
# sidebp<-ggplot(PP100Lstat[order(PP100Lstat$type, decreasing = T),], aes(x = type, y = meanPP, col = type))+
#   geom_boxplot(aes(fill = type), notch = T, col = "black")+
#   #scale_color_manual(values = c("#FFA800", "#FFF0BD" , "grey80"))+
#   scale_fill_manual(values = c("#FFA800", "#FFF0BD" , "grey80"))+
#   theme_classic()
# 
# ggplot(PP100Lstat[order(PP100Lstat$type, decreasing = T),])+
#   geom_density(aes(x = as.numeric(as.character(bin)), y = meanPP, col = type, group = type))
# 
# plot100PP<-ggplot(PP100Lstat[order(PP100Lstat$type, decreasing = T),], 
#                   aes(x = as.numeric(as.character(bin)), y = meanPP, col = type))+
#   geom_point()+
#   #geom_xsidedensity(aes(y = after_stat(density), fill = type), alpha = 0.3) +
#   #geom_ysidedensity(aes(x = after_stat(density), fill = type), alpha = 0.3) +
#   scale_color_manual(values = c("#FFA800", "#FFF0BD" , "grey80"))+
#   scale_fill_manual(values = c("#FFA800", "#FFF0BD" , "grey80"))+
#   scale_x_log10()+
#   xlab("exon length bin (midpoint, bin size = 3bp)")+
#   ylab("average exon-level PhyloP score")+
#   #ggtitle("PhyloP 100-ways (Vertebrates)")+
#   theme_classic()#+
#   #theme(legend.position = "inside",legend.position.inside = c(0.12,0.18))
# 
# cowplot::plot_grid(plot100PP, sidebp, nrow = 1, axis = "b", align = "h")
# 
# 
# cowplot::plot_grid(plot30PP, plot100PP, nrow = 2)
# ggsave("output/02_Conservation/Figure2B_PhyloP100-30_distros.pdf", width = 16, height = 9)
