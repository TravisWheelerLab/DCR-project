#!/usr/bin/env Rscript

# libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(pheatmap)
library(viridis)
library(grid)

### data #######################
junctions<-read.table("output/03_GTEx/demultiplexed.tallyman.counts.unique.wHeader.tsv", header=TRUE)
gtexMetadata<-read.table("data/GTExMetadataUpdated.txt", header=TRUE, comment.char = "")
#names(gtexMetadata)=c("sample", "tissue", "lib_size", "color", "ShortName", "NormFactors")

### color palette ##############
ann_colors<-list(
  tissue=c(Adipose_Subcutaneous="#FF6600",
           Adipose_Visceral="#FFAA00",
           Adrenal_Gland="#33DD33",
           Artery_Aorta="#FF5555",
           Artery_Coronary="#FFAA99",
           Artery_Tibial="#FF0101",
           Bladder="#AA0000",
           Brain_Amygdala="#EEEE00",
           Brain_Anterior_cingulate_cortex="#EEEE00",
           Brain_Caudate="#EEEE00",
           Brain_Cerebellar_Hemisphere="#EEEE00",
           Brain_Cerebellum="#EEEE00",
           Brain_Frontal_Cortex="#EEEE00",
           Brain_Hippocampus="#EEEE00",
           Brain_Hypothalamus="#EEEE00",
           Brain_Nucleus_accumbens="#EEEE00",
           Brain_Putamen="#EEEE00",
           Brain_Spinal_cord="#EEEE00",
           Brain_Substantia_nigra="#EEEE00",
           Breast_Mammary_Tissue="#00CDCD",
           Cells_Cultured_fibroblasts="#AAEEFF",
           Cells_EBV_transformed_lymphocytes="#EE82EE",
           Cervix_Ectocervix="#EED5D2",
           Cervix_Endocervix="#CCAADD",
           Colon_Sigmoid="#EEBB77",
           Esophagus_Gastroesophageal_Junction="#8B7355",
           Esophagus_Mucosa="#552200",
           Esophagus_Muscularis="#BB9988",
           Fallopian_Tube="#FFCCCC",
           Heart_Atrial_Appendage="#9917FF",
           Heart_Left_Ventricle="#660799",
           Kidney_Cortex="#21FFDD",
           Minor_Salivary_Gland="#99BB88",
           Muscle_Skeletal="#AAAAFF",
           Nerve_Tibial="#FFD700",
           Pancreas="#995522",
           Pituitary="#AAFF99",
           Prostate="#D9D9D9",
           Skin_Not_Sun_Exposed="#3A5FCD",
           Skin_Sun_Exposed="#1E90FF",
           Small_Intestine_Terminal_Ileum="#555622",
           Spleen="#778855",
           Stomach="#FFD39B",
           Thyroid="#008B45",
           Uterus="#FD66FF",
           Vagina="#FF5599",
           Brain_cortex="#EEEE00",
           Colon_transverse="#CC9955",
           Liver="#AABB66",
           Lung="#99FF01",
           Ovary="#FEAAFF",
           Testis="#A6A6A6",
           Whole_blood="#FF00BB")
)

### heatmap function ###############
DCE_heatmaps<-function(input_mat=AllnoProbes, GTEx_metadata=GTEx_metadata, min_reads = 2, min_samples = 10, clustering = T, per_sample = F){
  AllnoProbes<-input_mat
  #print(head(AllnoProbes[1:10,1:10]))
  # Subset your count data
  count_data <- unique(AllnoProbes[, -1])  # We remove duplicate needles due to demultiplexing
  rownames(count_data)<-AllnoProbes[,1]
  print(head(count_data[1:10,1:10]))
  count_data <- as.data.frame(lapply(count_data, as.numeric))
  to_keep<-apply(count_data, 1, function(x){length(x[x>min_reads])}) > min_samples # filter needles based on expresion thresholds
  #print(to_keep)
  count_data<-count_data[to_keep,]
  kept_samples<-data.frame(needle = AllnoProbes[,1], kept = to_keep)
  write.table(kept_samples, "output/03_GTEx/filtered_counts.tsv", quote = FALSE, sep = "\t", row.names = TRUE)
  samples <- colnames(count_data)
  
  # Subset metadata for the samples in count_data
  metadata <- GTEx_metadata %>%
    filter(sample %in% samples) %>%
    arrange(match(sample, samples))
  
  # Ensure the sample order is the same in both tables
  count_data <- count_data[, metadata$sample]
  
  # Normalize the counts by the library sizes
  Nlib_sizes <- metadata$lib_size * metadata$NormFactors
  #done$CPM <- (done$count / (done$lib_size * done$NormFactors) * 1e6)
  normalized_counts <- sweep(count_data, 2, Nlib_sizes, FUN = "/") * 1e6
  
  # grouped by tissue
  if(per_sample == F){
    # Transpose to have samples as rows (so we can join metadata)
    normalized_counts_df <- as.data.frame(t(normalized_counts))
    
    # Add tissue information from the metadata to the normalized counts
    normalized_counts_df$sample <- rownames(normalized_counts_df)
    normalized_counts_df <- merge(normalized_counts_df, metadata, by = "sample")
    #tail(names(normalized_counts_df))
    # write normalized counts to file to extract statistics for paper
    write.table(normalized_counts_df, "output/03_GTEx/normalized_counts.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
    
    # Group by tissue and calculate mean normalized counts per tissue
    mean_normalized_by_tissue <- normalized_counts_df %>%
      select(-sample, -lib_size, -color, -ShortName) %>%  # Remove non-numeric columns
      group_by(tissue) %>%
      summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = 'drop')  # Ensure to drop the grouping
    
    # Convert 'tissue' column to row names (using base R)
    mean_normalized_by_tissue<-data.frame(mean_normalized_by_tissue)
    row.names(mean_normalized_by_tissue) <- mean_normalized_by_tissue$tissue
    mean_normalized_by_tissue<-t(mean_normalized_by_tissue[,-1])
    
    # Ensure it's a numeric matrix before transposing
    mean_normalized_matrix <- as.matrix(mean_normalized_by_tissue)
    
    # Compute the Pearson correlation matrix between tissues
    correlation_matrix <- cor(mean_normalized_matrix, method = "pearson")
    
    # get the colors
    t_colors<-metadata %>% group_by(tissue) %>%
      summarize(tissue = unique(tissue),
                color = unique(as.factor(color)),
                ) %>% data.frame()
    row.names(t_colors)<-t_colors$tissue
    t_colors$color<-NULL
    
    # get the shortnames
    t_short<-metadata %>% group_by(tissue) %>%
      summarize(tissue = unique(tissue),
                short = unique(ShortName)) %>% data.frame()
    row.names(t_short)<-t_short$tissue
    
    if(clustering == T){ 
      # Plot the heatmap using the viridis palette
      pheatmap(correlation_matrix, 
               clustering_distance_rows = "euclidean", 
               clustering_distance_cols = "euclidean", 
               color = viridis(50),  # Use viridis color scale
               main = paste("Tissue Pearson Correlation Heatmap | min.", min_reads, "reads in at least", min_samples, "samples" ),
               annotation_row = t_colors, annotation_colors = ann_colors,
               annotation_col = t_colors,
               labels_col = t_short$short
      )
    } else {
      # unclustered
      pheatmap(correlation_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
               #clustering_distance_rows = "euclidean", 
               #clustering_distance_cols = "euclidean", 
               color = viridis(50),  # Use viridis color scale
               main = paste("Tissue Pearson Correlation Heatmap | min.", min_reads, "reads in at least", min_samples, "samples" ),
               annotation_row = t_colors, annotation_colors = ann_colors,
               annotation_col = t_colors, cellheight = 10, cellwidth = 10
      )
    }
    
    # per sample
  } else {
    ### if we want for each for each replicate
    replicates_cor<-cor(as.matrix(normalized_counts), method = "pearson")
    rep_col=metadata[,c("sample", "tissue")]
    #row.names(rep_col)=rep_col$sample
    rep_col=data.frame(rep_col[,-1], row.names = rep_col[,1])
    names(rep_col)<-"tissue"
    # rep_col=rep_col[,-c(1)]
    # rep_colors=metadata[,c("tissue", "color")]
    #row.names(rep_colors)=rep_col$sample
    if(clustering == T){
      pheatmap(replicates_cor, 
               clustering_distance_rows = "euclidean", 
               clustering_distance_cols = "euclidean", 
               color = viridis(50),  # Use viridis color scale
               main = paste("Samples Pearson Correlation Heatmap | min.", min_reads, "reads in at least", min_samples, "samples" ),
               annotation_row = rep_col,
               annotation_col = rep_col,
               annotation_colors = ann_colors,show_rownames = FALSE,show_colnames = FALSE
      )
    } else {
      # without clustering
      pheatmap(replicates_cor, cluster_rows = FALSE, cluster_cols = FALSE,
               #clustering_distance_rows = "euclidean", 
               #clustering_distance_cols = "euclidean", 
               color = viridis(50),  # Use viridis color scale
               main = paste("Samples Pearson Correlation Heatmap | min.", min_reads, "reads in at least", min_samples, "samples" ),
               annotation_row = rep_col,
               annotation_col = rep_col,
               annotation_colors = ann_colors,
               show_rownames = FALSE,show_colnames = FALSE
      )
    }
  }
}

### run ############################
# 2 reads in 2 sample to keep needle
pdf(width = 20, height = 15, file= "output/03_GTEx/Heatmap.pdf")
# note: the layout and annotation are optimized for the parameters below.
#       If clustering = F, annotation will be with full names on the bott
#       -om row
heatG<-DCE_heatmaps(input_mat = junctions, GTEx_metadata = gtexMetadata,
             min_reads = 2, min_samples = 2,clustering = T)
grid.draw(heatG$gtable)
dev.off()
