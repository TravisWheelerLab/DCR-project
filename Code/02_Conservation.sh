#!/usr/bin/env bash

#####################################
# 02 Conservation of DCR | Figure 2 #
#####################################

# create output and log dirs of do not exist
mkdir -p output/02_Conservation
mkdir -p logs/02_Conservation

############## SCRIPT STARTS HERE ##############

#### FIGURE 2 A ####

## 0. convert the bed into a UCSC compatible bed file
awk '{print $1"\t"$2"\t"$3"\tDCRid="$4";nFrames="$7"\t.\t"$5}' output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed > output/02_Conservation/filtered_DCR_overlapping_regions_to_LO.bed

## 1. LifOver the human DCR bed regions to the mouse genome
DCR2LIFT="output/02_Conservation/filtered_DCR_overlapping_regions_to_LO.bed"
CHAINFILE="data/hg38ToMm39.over.chain"
OUT1="output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed"
OUT2="output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.unlifted.bed"
INTIN=$(wc -l < ${DCR2LIFT} | awk '{print $1}')
DCRIN=$(cut -f 4 ${DCR2LIFT} | sed 's/;/\t/g' | cut -f 1 | sort | uniq | wc -l | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] The input file contains ${INTIN} intervals for ${DCRIN} DCRs"
# LO
~/Documents/liftOver ${DCR2LIFT} ${CHAINFILE} ${OUT1} ${OUT2}
# assess how many regions were successfully lifted over
INOUT=$(wc -l output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed | awk '{print $1}')
DCROUT=$(cut -f 4 output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed | sed 's/;/\t/g' | cut -f 1 | sort | uniq | wc -l | awk '{print $1}')
INTMISS=$(grep -v '#' output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.unlifted.bed | wc -l | awk '{print $1}')
DCRinMISS=$(grep -v '#' output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.unlifted.bed | cut -f 4 | sed 's/;/\t/g' | cut -f 1 | sort | uniq | wc -l | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] After LO there are ${INOUT} intervals for ${DCROUT} DCRs."
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] ${INTMISS} intervals, part of ${DCRinMISS} DCRs were not lifted over "

# [info 2026-02-17 15:57:08] The input file contains 1741 intervals for 1407 DCRs
# [info 2026-02-17 15:57:24] After LO there are 1603 intervals for 1295 DCRs.
# [info 2026-02-17 15:57:24] 138 intervals, part of 134 DCRs were not lifted over

## 2. Convert to fasta and translate
OUT3="output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa"
OUT4="output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa_translated.fa"
OUT5="output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv"

# grab fasta
bedtools getfasta -fi data/mm39.fa -bed ${OUT1} -s > ${OUT3}
# translate (top strand only)
esl-translate -l 1 --watson ${OUT3} > ${OUT4}
# convert to table so I can do stats | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' # from https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${OUT4} | \
awk 'getline seq {print $0"\t"seq}' | \
sed 's/>//g;s/source=//g;s/coords=//g;s/\.\./\t/g;s/length=//g;s/frame=//g;s/desc=//g;s/ /\t/g;s/(-)/\t-/g;s/(+)/\t+/g' | \
awk '{split($2,coord,"[(,),:,-]"); print $1"\t"coord[1]"\t"coord[2]"\t"coord[3]"\t"$3"\t"(coord[3]-coord[2])"\t"$4"\t"$5"\t"($5-$4)"\t"(($5-$4)/(coord[3]-coord[2]))"\t"$6"\t"$7"\t"$8}' | \
sort -k2,2 -k3,3n -k4,4n > ${OUT5}


## 3. Generate random intervals on mouse exons
# get the latest mm39 ENSEMBL gene annotation 
wget https://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz -O data/Mus_musculus.GRCm39.115.gtf.gz
# get list of human/mice orthologs for mm39 by ENSEMBL
# http://www.ensembl.org/biomart/martview/08b0c6e2ab12f5364483585fdeb9b711?VIRTUALSCHEMANAME=default&ATTRIBUTES=mmusculus_gene_ensembl.default.homologs.ensembl_gene_id|mmusculus_gene_ensembl.default.homologs.ensembl_gene_id_version|mmusculus_gene_ensembl.default.homologs.hsapiens_homolog_associated_gene_name|mmusculus_gene_ensembl.default.homologs.hsapiens_homolog_ensembl_gene|mmusculus_gene_ensembl.default.homologs.hsapiens_homolog_chromosome&FILTERS=&VISIBLEPANEL=resultspanel
# file downloaded and placed in data/human-mice-orthologs.tsv 
# filter out non standard human chromosomes and human DCR genes and make a list of mice orthologs
grep -vwf <(cut -f 6 output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed | awk '{print toupper($0)}' | sort | uniq) <(awk '$5 !~ /GL|KI|MT/ && NF == 5' data/human-mice-orthologs.tsv) | awk '{print $1}' > output/02_Conservation/mm39.genes.orthologs.of.non-DCR.genes.txt
# make a bed file of mice mm39 CDS for these control exons
gzcat data/Mus_musculus.GRCm39.115.gtf.gz | awk '$3 == "CDS"' | awk '{gsub("\"",""); gsub(";",""); print "chr"$1"\t"$4-1"\t"$5"\tgene="$10";exon="$18"\t"$6"\t"$7}' | grep -f output/02_Conservation/mm39.genes.orthologs.of.non-DCR.genes.txt > output/02_Conservation/mm39.exons.nonOrthoToHumanDCRintervals.bed

# generate random intervals on mouse exons (1-to-1 orthologous non-DCE exons, and with the same region size distro that DCE regions)
OUT6="output/02_Conservation/fake_DCE_intervals.bed"
python3.11 scripts/02_Conservation/random_CE_regions.py -n 100000 -e  output/02_Conservation/mm39.exons.nonOrthoToHumanDCRintervals.bed -d output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed -f data/mm39.fa


## 4. Convert random intervals to fasta and translate
OUT8="output/02_Conservation/fake_DCE_intervals.fa_translated.tsv"
# grab fasta
OUT7="output/02_Conservation/fake_DCE_intervals.fa"
bedtools getfasta -fi data/mm39.fa -bed <(awk 'NR > 1' ${OUT6} | sort -k1,1 -k2,2n) -s > ${OUT7}
# translate (top strand only)
esl-translate -l 1 --watson ${OUT7} > ${OUT7}.translated.fa
# convert to table so I can do stats | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' # from https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${OUT7}.translated.fa | \
awk 'getline seq {print $0"\t"seq}' | \
sed 's/>//g;s/source=//g;s/coords=//g;s/\.\./\t/g;s/length=//g;s/frame=//g;s/desc=//g;s/ /\t/g;s/(-)/\t-/g;s/(+)/\t+/g' | \
awk '{split($2,coord,"[(,),:,-]"); print $1"\t"coord[1]"\t"coord[2]"\t"coord[3]"\t"$3"\t"(coord[3]-coord[2])"\t"$4"\t"$5"\t"($5-$4)"\t"(($5-$4)/(coord[3]-coord[2]))"\t"$6"\t"$7"\t"$8}' | \
sort -k2,2 -k3,3n -k4,4n > output/02_Conservation/fake_DCE_intervals.fa_translated.tsv

## 5. Shuffle mouse exons for random model
# make a bed with mouse phased mouse CDS > 21 nt (i.e. CDS trimmed to start and end at codons edges, so they can be shuffled)
python3.11 scripts/02_Conservation/extract_codon_phased_CDS.py
# keep < 500 bp, shuffle codons and convert to fasta. We do 503 to tolerate stop codons that will be removed.
awk '$3-$2 <= 503' output/02_Conservation/mm39_CDS_codon_phased_21bpMin_stdChrom.bed > output/02_Conservation/mm39_CDS_codon_phased_21bpMin_stdChrom_500less.bed
python3.11 scripts/02_Conservation/generate_fake_mouse_exons.py
awk '/>/{split($1,header,"_"); sub("exon", "", header[3]); print ">"header[2]"::"header[3]"_shuffled"; next} !/>/{print $0}' output/02_Conservation/fake_mm39_exons.fa > output/02_Conservation/fake_mm39_exons_formatted.fa
# translate and format
OUT9="output/02_Conservation/fake_mm39_exons_formatted.fa"
OUT10="output/02_Conservation/fake_mm39_exons.fasta_translated.fa_linearized_renamedWithCoord.tsv"
# translate (top strand only)
esl-translate -l 1 --watson ${OUT9} > ${OUT9}.translated.fa
# get real length of fake exon
awk 'BEGIN {OFS = "\n"}; /^>/ {print(substr(sequence_id, 2)" "sequence_length); sequence_length = 0; sequence_id = $0}; /^[^>]/ {sequence_length += length($0)}; END {print(substr(sequence_id, 2)" "sequence_length)}' ${OUT9} > ${OUT9}.len
# accolate length to translated orf name
join -12 -21 <(awk '/>/ {sub(">",""); sub("source=",""); print $1"\t"$2}' ${OUT9}.translated.fa | sort -k2,2) <(grep 'shuffled' ${OUT9}.len | sort -k1,1) | sed 's/ /\t/g' > output/02_Conservation/fake_mm39_exons.len2orf.correspondance.txt
# redo the tsv
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${OUT9}.translated.fa | awk 'getline seq {print $0"\t"seq}' > ${OUT9}.translated.fa_linearized
join -11 -22 <(awk '{sub(">", ""); print $0}' ${OUT9}.translated.fa_linearized | sort -k1,1) <(sort -k2,2 output/02_Conservation/fake_mm39_exons.len2orf.correspondance.txt) > output/02_Conservation/fake_mm39_exons.fasta_translated.fa_linearized_with_sizes
# final:
awk '{sub("source=",""); sub("coords=", ""); sub("\\.\\.", "\t"); sub("length=",""); sub("frame=", ""); print $1"\t"$2"\tNA\tNA\tNA\t"$NF"\t"$3"\t"$4"\t"$4-$3"\t"($4-$3)/$NF"\t"$5"\t"$6"\t"$8}' output/02_Conservation/fake_mm39_exons.fasta_translated.fa_linearized_with_sizes | sort -k2,2 > output/02_Conservation/fake_mm39_exons.fasta_translated.fa_linearized_renamedWithCoord.tsv


#### FIGURE 2 A and B ####  

## 6. Create the mouse conservation plot ################

# create correspondance between fake mouse exons and DCRs
bedtools intersect -wao -a <(cut -f 2- output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv) -b output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed | awk '$14 != "-1" {print "DCE:"$1":"$2":"$3"\t"$9"\t"$16}' | sed 's/DCRid=//g;s/;nFrames=/\t/g' | cut -f 1-3 > output/02_Conservation/mouse_exon_to_dce_id.txt
# plot
Rscript scripts/02_Conservation/Mice_conservation.r -d output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv -n output/02_Conservation/fake_DCE_intervals.fa_translated.tsv -e output/02_Conservation/fake_mm39_exons.fasta_translated.fa_linearized_renamedWithCoord.tsv

# [1] "Total mapped DCR exons: 1599"
# [1] "Total mapped DCR to mouse exons: 1295"
# [1] "Total mapped DCR exons >21nt: 1486"
# [1] "Total mapped DCR to mouse exons (>21nt): 1232"
# [1] "DCRexons passing filter: 1216"

# --> DCR too shorts = 1295 - 1232 = 63 --> OK
# --> DCR not mapped = 1407 - 1295 = 112
# --> not DCR in mouse = not_DCR + partial = 176 + 74 = 250 ==> 0.203 --> 20%
# --> DCR in mouse = 982
# --> total = 1407 --> OK!

       # DCR    not_DCR not_mapped    partial  too_short
       # 982        176        112         74         63

#### FIGURE 2 C ####

# format filtered DCR bed properly
awk '{print $1"\t"$2"\t"$3"\t"$4"\t\.\t"$5}' output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed > output/01_Description/filtered_DCR_overlapping_regions_with_dceID_formatted.bed

echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Running PhyloP workflow..."
GENECODE_CDS="data/hg38.GENECODE.V48.bed"
DCRregions="output/01_Description/filtered_DCR_overlapping_regions_with_dceID_formatted.bed"
output_dir="output/02_Conservation"
ENSG2ENST="data/correspondance_ENSG_ENST_v114.hg38.14"
PPOUT1="output/02_Conservation/DCR_exons_DCR_genes_UNIQ.bed"
PPOUT2="output/02_Conservation/nonDCR_exon_DCR_gene_UNIQ.bed"
PPOUT3="output/02_Conservation/nonDCR_exon_nonDCR_gene_UNIQ.bed"

# file prep for PhyloP

DCRin=$(cut -f 4 ${DCRregions} | sed 's/DCRid=//g;s/;nFrames/\t/g' | cut -f 1 | sort | uniq | wc -l | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] DCR in the input file: ${DCRin}"
## 1. intersect DCE regions with the CDS annotation file #
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Intersecting DCE regions with the CDS annotation file"
bedtools intersect -wao -b ${DCRregions} -a ${GENECODE_CDS} | awk '$8 != "-1" && $6 == $12' | cut -f 1-4,10,13 |\
sed 's/DCRid=//g;s/;nFrames=/\t/g' > ${output_dir}/DCE_CDS_intersect.bed
# count the number of DCR in the intersected file
DCR_intersect=$(cut -f 5 ${output_dir}/DCE_CDS_intersect.bed | sort | uniq | wc -l | awk '{print $1}')
difference=$((DCRin - DCR_intersect))
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Number of DCR in the intersected file: ${DCR_intersect} (difference with filtered = ${difference})"
# drop exons that < 20bp
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Dropping exons that are smaller than 20bp"
awk '{if ($3 - $2 >= 20) print $0}' ${output_dir}/DCE_CDS_intersect.bed > ${output_dir}/DCE_CDS_intersect_filtered.bed
# count the number of DCR in the filtered file
DCR_size_filtered=$(cut -f 5 ${output_dir}/DCE_CDS_intersect_filtered.bed | sort | uniq | wc -l | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Number of DCR in the size filtered file: ${DCR_size_filtered} (difference with filtered = $((DCRin - DCR_size_filtered)))"

# [info 2026-02-17 19:33:34] DCR in the input file: 1407
# [info 2026-02-17 19:33:34] Intersecting DCE regions with the CDS annotation file
# [info 2026-02-17 19:33:37] Number of DCR in the intersected file: 1387 (difference with filtered = 20)
# [info 2026-02-17 19:33:37] Dropping exons that are smaller than 20bp
# [info 2026-02-17 19:33:37] Number of DCR in the size filtered file: 1386 (difference with filtered = 21)

## 2. Add gene annotation #
# link each exon to GENE ID using the transcript ID
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Adding gene annotation..."
awk '{split($4, a, "\."); print a[1]"\t"$0}' ${output_dir}/DCE_CDS_intersect_filtered.bed > ${output_dir}/DCE_CDS_intersect_filtered_with_transcriptID.bed
# join Gene ID
join -12 -21 <(sort -k2,2 ${ENSG2ENST}) <(sort -k1,1 ${output_dir}/DCE_CDS_intersect_filtered_with_transcriptID.bed) | sed 's/ /\t/g' | awk '{print $2"\tDCR"}' | sort | uniq > ${output_dir}/DCR_geneID.txt
DCRwGname=$(join -12 -21 <(sort -k2,2 ${ENSG2ENST}) <(sort -k1,1 ${output_dir}/DCE_CDS_intersect_filtered_with_transcriptID.bed) | sed 's/ /\t/g' | cut -f 7 | sort | uniq | wc -l | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Number of DCR with gene annotation: ${DCRwGname} (difference with filtered = $((DCRin - DCRwGname)))"
# also add gene annotation to the whole exon files
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] adding gene annotation to the whole exon files"
bedtools intersect -wao -b ${DCRregions} -a ${GENECODE_CDS} | awk '$8 == "-1"' | cut -f 1-4 > ${output_dir}/nonDCE_CDS_intersect.bed
join -12 -21 <(sort -k2,2 ${ENSG2ENST}) <(awk '{split($4, a, "\."); print a[1]"\t"$0}' ${output_dir}/nonDCE_CDS_intersect.bed | sort -k1,1) | sed 's/ /\t/g' > ${output_dir}/nonDCE_CDS_intersect_WgeneID.txt

# [info 2026-02-17 19:34:33] Adding gene annotation...
# [info 2026-02-17 19:34:37] Number of DCR with gene annotation: 1386 (difference with filtered = 21)
# [info 2026-02-17 19:34:37] adding gene annotation to the whole exon files

## 3. Extract nonDCE regions
# extract nonDCR_exon_DCR_gene
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Extracting nonDCR_exon_DCR_gene"
join -11 -22 <(sort -k1,1 ${output_dir}/DCR_geneID.txt) <(sort -k2,2 ${output_dir}/nonDCE_CDS_intersect_WgeneID.txt) | sed 's/ /\t/g' | cut -f 4- > ${output_dir}/nonDCR_exon_DCR_gene.bed
# extract nonDCR_exon_nonDCR_gene
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Extracting nonDCR_exon_nonDCR_gene"
join -a1 -12 -21 <(sort -k2,2 ${output_dir}/nonDCE_CDS_intersect_WgeneID.txt)  <(sort -k1,1 ${output_dir}/DCR_geneID.txt) | sed 's/ /\t/g' | grep -v 'DCR' | cut -f 3- > ${output_dir}/nonDCR_exon_nonDCR_gene.bed

## Step 4: remove duplicate exons #
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Removing duplicate exons from files..."
# remove duplicate exons in the DCR file
awk '{print $1"_"$2"_"$3"\t"$0}' ${output_dir}/DCE_CDS_intersect_filtered.bed | sort -k2,2 -k3,3n | sort -u -k1,1 | cut -f 2- > ${output_dir}/DCR_exons_DCR_genes_UNIQ.bed
beforePPcount=$(cut -f 5 output/02_Conservation/DCR_exons_DCR_genes_UNIQ.bed | sort | uniq | wc -l | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Number of DCR exons before PhyloP: ${beforePPcount}"
# remove duplicate exons in the nonDCR files
awk '{print $1"_"$2"_"$3"\t"$0}' ${output_dir}/nonDCR_exon_DCR_gene.bed | sort -k2,2 -k3,3n | sort -u -k1,1 | cut -f 2- > ${output_dir}/nonDCR_exon_DCR_gene_UNIQ.bed
awk '{print $1"_"$2"_"$3"\t"$0}' ${output_dir}/nonDCR_exon_nonDCR_gene.bed | sort -k2,2 -k3,3n | sort -u -k1,1 | cut -f 2- > ${output_dir}/nonDCR_exon_nonDCR_gene_UNIQ.bed

# [info 2026-02-17 19:36:00] Removing duplicate exons from files...
# [info 2026-02-17 19:36:00] Number of DCR exons before PhyloP: 1386

# 5. Run PhyloP
PPOUT4_1="output/02_Conservation/DCEvsNonDCE.meanPhyloP100.out"
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Running PhyloP per exon..."
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] PP 100-ways..."
python3.11 scripts/02_Conservation/getPhyloPmeanScore.py data/hg38.phyloP100way.bw > output/02_Conservation/DCEvsNonDCE.meanPhyloP100.out
#   .y.         group1             group2                 n1    n2 statistic    df        p p.adj p.adj.signif y.position groups        xmin  xmax
#   <chr>       <chr>              <chr>               <int> <int>     <dbl> <dbl>    <dbl> <dbl> <chr>             <dbl> <named list> <dbl> <dbl>
# 1 mean_PhyloP DCE_exon           nDCE_exon_DCE_gene   3315  3285     -2.58 6593. 0.01     0.02  *                  8.40 <chr [2]>        1     2
# 2 mean_PhyloP DCE_exon           nDCE_exon_nDCE_gene  3315  3312     -3.57 6614. 0.000358 0.001 **                 8.59 <chr [2]>        1     3
# 3 mean_PhyloP nDCE_exon_DCE_gene nDCE_exon_nDCE_gene  3285  3312     -1.01 6594. 0.314    0.314 ns                 8.78 <chr [2]>        2     3

# 6. Plot PP graphs
PPOUT5="output/02_Conservation/Figure2B_PhyloP_boxplot-30-ways.pdf"
PPOUT6="output/02_Conservation/Figure2B_PhyloP_boxplot-100-ways.pdf"
PPOUT7="output/02_Conservation/Figure2B_PhyloP100-30_distros.pdf"
if [ "$FLAG" = true ] && [ -s "$PPOUT5" ] && [ -s "$PPOUT6" ] && [ -s "$PPOUT7" ]; then
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping PhyloP boxplot generation as output files already exist and are not empty."
else
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Generating PhyloP plots..."
  Rscript scripts/02_Conservation/PhyloP_boxplots.r output/02_Conservation/DCEvsNonDCE.meanPhyloP100.out "mean PhyloP 100-ways conservation score" 100
  Rscript scripts/02_Conservation/PhyloP_boxplots.r output/02_Conservation/DCEvsNonDCE.meanPhyloP30.out "mean PhyloP 30-ways conservation score" 30
  Rscript scripts/02_Conservation/PP_distros.r
fi

# 7. Get DCE PP ranking and mouse status [Figure 2 C]
# calculate score per DCR
PPOUT8="output/02_Conservation/meanPP_per_DCE.tsv"
if [ "$FLAG" = true ] && [ -s "$PPOUT8" ]; then
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping mean PP per DCR calculation as output file already exists and is not empty."
else
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Calculate mean PP per DCR @100-Ways."
  python3.11 scripts/02_Conservation/getPhyloPmeanScore_perDCE.py
fi

# 8. Plot DCR PP outliers
PPOUT9="output/02_Conservation/Figure2C_DCR_PhyloP_outliers.pdf"
if [ "$FLAG" = true ] && [ -s "$PPOUT9" ]; then
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping DCR PhyloP outliers plot generation as output file already exists and is not empty."
else
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Generating DCR PhyloP outliers plot"
  # get correspondance DCR to gene
  awk 'NR > 1 {print $1"\t"$2}' output/01_Description/homo_sapiens.tsv | sort | uniq > output/02_Conservation/DCEtoGene.txt
  # get frame per DCR
  awk 'NR > 1' output/01_Description/homo_sapiens.tsv | cut -f 1,4 | sort | uniq | cut -f 1 | sort | uniq -c | awk '{print $2"\t"$1}' > output/02_Conservation/frames_per_DCR.txt
  # get mouse exon to DCR correspondance
  bedtools intersect -wao -a <(cut -f 2- output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv) -b output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed | awk '$14 != "-1" {print "DCE:"$1":"$2":"$3"\t"$9"\t"$16}' | sed 's/DCRid=//g;s/;nFrames=/\t/g' | cut -f 1-3 > output/02_Conservation/mouse_exon_to_dce_id.txt 
  # plot
  Rscript scripts/02_Conservation/PPoutliers.r
fi

echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Conservation analysis completed"