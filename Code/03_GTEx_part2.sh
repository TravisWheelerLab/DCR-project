######################
# 03 GTEx | Figure 3 #
######################

mkdir -p output/03_GTEx/

#### Data prep #######

# demultiplex the count and reattribute the junction names
join -11 -22 <(sort -k1,1 data/final_output.tsv) <(awk '{sub(/>/, "", $1); print $1"\t"$3}' output/03_GTEx/junctions_name_correspondance_for_demultiplexing.tsv | sort -k2,2) | awk '{printf "%s", $NF; for (i = 2; i < NF; i++) printf " %s", $i; printf "\n";}' > output/03_GTEx/demultiplexed.tallyman.counts.tsv
#NOTE: will need a header!

#### Figure S5 (GTEx heatmap) ####

# keep only unique count for heatmap
awk '{key = ""; for (i = 2; i <= NF; i++) key = key $i OFS;  if (!seen[key]++) print}' output/03_GTEx/demultiplexed.tallyman.counts.tsv > output/03_GTEx/demultiplexed.tallyman.counts.unique.tsv
# add header
cat <(head -n 1 data/final_output.tsv) output/03_GTEx/demultiplexed.tallyman.counts.unique.tsv > output/03_GTEx/demultiplexed.tallyman.counts.unique.wHeader.tsv

# Filter for final set of DCR
cat <(head -n 1 data/demultiplexed.tallyman.counts.unique.wHeader.tsv) <(grep -f <(cut -f 4 output/01_Description/filtered_DCR_overlapping_regions_with_dceID_formatted.bed | sort | uniq | awk '{print "^"$0"_"}') data/demultiplexed.tallyman.counts.unique.wHeader.tsv) > output/03_GTEx/demultiplexed.tallyman.counts.unique.wHeader.tsv

# plot
Rscript scripts/03_GTEx/heatmaps.r
# count genes in heatmap
awk '{print $2}' output/03_GTEx/filtered_counts.tsv | sed 's/_/\t/g' | awk '{print $2}' | sort | uniq | wc -l
#    1263

#### Differential Usage Analysis ####
Rscript scripts/03_GTEx/DiffUsageTest.R
sed 's/ /\t/g' output/03_GTEx/output_DCE_All_frames.tsv | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8$9$10}' > output/03_GTEx/output_DCE_All_frames_simple.tsv


# get gene approximate position for plot 
sed 's/\|/\t/g;s/\//\t/g' output/01_Description/UniprotKB_isoforms_mapping.bed | awk '{print $5"\t"$1"\t"$2}' | sort -u -k1,1 > output/03_GTEx/gene_to_pos.txt
join -11 -22 <(sort -k1,1 output/03_GTEx/gene_to_pos.txt) <(sort -k2,2 output/02_Conservation/DCEtoGene.txt) | sed 's/ /\t/g' > output/03_GTEx/gene_chr_pos_DCRid.tsv

# Figure 4 # 

# Figure 4 is generated with the script scripts/3_GTEx/GTEx_Plot_DotMatrix.R executed manually

# find genes with Diff expr + NMD 
grep -wf <(grep -wf <(grep "\*\*\*" output/03_GTEx/output_DCE_All_frames_simple.tsv | cut -f 1 | sed 's/_/\t/g' | cut -f 1 | sort | uniq) data/Table_S1.tsv | grep -w NMD | awk '{print $2}' | sort | uniq) output/01_Description/filtered_homo_sapiens_DCR.tsv | cut -f 2 | sort | uniq