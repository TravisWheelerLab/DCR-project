######################
# 03 GTEx | Figure 3 #
######################

#### Data prep #######

# demultiplex the count and reattribute the junction names
join -11 -22 <(sort -k1,1 data/final_output.tsv) <(awk '{sub(/>/, "", $1); print $1"\t"$3}' output/03_GTEx/junctions_name_correspondance_for_demultiplexing.tsv | sort -k2,2) | awk '{printf "%s", $NF; for (i = 2; i < NF; i++) printf " %s", $i; printf "\n";}' > output/03_GTEx/demultiplexed.tallyman.counts.tsv
#NOTE: will need a header!

#### Figure S5 (GTEx heatmap) ####

# keep only unique count for heatmap
awk '{key = ""; for (i = 2; i <= NF; i++) key = key $i OFS;  if (!seen[key]++) print}' output/03_GTEx/demultiplexed.tallyman.counts.tsv > output/03_GTEx/demultiplexed.tallyman.counts.unique.tsv
# add header
cat <(head -n 1 data/final_output.tsv) output/03_GTEx/demultiplexed.tallyman.counts.unique.tsv > output/03_GTEx/demultiplexed.tallyman.counts.unique.wHeader.tsv
# Addition 10-07-25: remove DCR with ambiguous isoforms mappings (wrong start or end due to small first or last exon)
cat <(head -n 1 output/03_GTEx/demultiplexed.tallyman.counts.unique.wHeader.tsv) <(grep -f <(cut -f 5 output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb.bed | sort | uniq | awk '{print "^"$0"_"}') output/03_GTEx/demultiplexed.tallyman.counts.unique.wHeader.tsv) > output/03_GTEx/demultiplexed.tallyman.counts.unique.wHeader.tsv_fixed
mv output/03_GTEx/demultiplexed.tallyman.counts.unique.wHeader.tsv_fixed output/03_GTEx/demultiplexed.tallyman.counts.unique.wHeader.tsv
# plot
Rscript scripts/03_GTEx/heatmaps.r

#### Differential Usage Analysis ####
Rscript scripts/03_GTEx/DiffUsageTest.R
sed 's/ /\t/g' output/03_GTEx/output_DCE_All_frames.tsv | cut -f 1-7 > output/03_GTEx/output_DCE_All_frames_simple.tsv

# get gene approximate position for plot (not used)
sed 's/\|/\t/g;s/\//\t/g' output/01_Description/UniprotKB_isoforms_mapping.bed | awk '{print $5"\t"$1"\t"$2}' | sort -u -k1,1 > output/03_GTEx/gene_to_pos.txt
join -11 -22 <(sort -k1,1 output/03_GTEx/gene_to_pos.txt) <(sort -k2,2 output/02_Conservation/DCEtoGene.txt) | sed 's/ /\t/g' > output/03_GTEx/gene_chr_pos_DCRid.tsv

### Others
sed 's/ : /:/g' output/03_GTEx/output_DCE_All_frames.tsv > output/03_GTEx/output_DCE_All_frames_fixedColumns.tsv