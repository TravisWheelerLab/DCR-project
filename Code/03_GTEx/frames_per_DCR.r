#!/usr/bin/env Rscript
# calculate frames per DCE for needles
df<-read.table("output/03_GTEx/homo_sapiens.IREs-removed.duplicates-removed.tsv")
fr_per_pair<-tapply(df$V4, df$V1, function(x){length(unique(x))})
df_F<-df[df$V1 %in% names(fr_per_pair[fr_per_pair > 1]),] # keep pairs with 2+ frames
write.table(df_F, "output/03_GTEx/homo_sapiens.IREs-removed.duplicates-removed.2plus-frames-per-pair.tsv", quote=F, row.names=F, col.names=F, sep="\t")