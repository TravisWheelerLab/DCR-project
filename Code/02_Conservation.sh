#!/usr/bin/env bash

#####################################
# 02 Conservation of DCR | Figure 2 #
#####################################

# Parse arguments
usage() {
    echo "Usage: $0 [-s]"
    echo "  -s    Skip command if output file exists and is not empty (default = false)"
    echo "  -i    Input dir (required)"
    exit 1
}

FLAG=false
INPUT_FILE=""

while getopts ":si:I:" opt; do
  case $opt in
    s)
      FLAG=true
      ;;
    \?)
      usage
      ;;
  esac
done

# create output and log dirs of do not exist
mkdir -p output/02_Conservation
mkdir -p logs/02_Conservation

############## SCRIPT STARTS HERE ##############

#### FIGURE 2 A ####

## 1. LifOver the human DCR bed regions to the mouse genome
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Lifting over DCR regions from human hg38 to mouse mm39 genome"
# DCRBED="output/01_Description/filtered_DCR_overlapping_regions.bed"
# DCR2LIFT="output/02_Conservation/merged_DCR_overlapping_regions.bed"
# for now, I keep all DCR regions, even if they overlap with other DCRs
DCR2LIFT="output/01_Description/filtered_DCR_overlapping_regions.bed"
CHAINFILE="data/hg38ToMm39.over.chain"
OUT1="output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed"
OUT2="output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.unlifted.bed"
INTIN=$(wc -l < ${DCR2LIFT} | awk '{print $1}')
DCRIN=$(cut -f 4 ${DCR2LIFT} | sed 's/;/\t/g' | cut -f 1 | sort | uniq | wc -l | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] The input file contains ${INTIN} intervals for ${DCRIN} DCRs"
# Merge the input file to avoid duplicates
#bedtools merge -s -c 4,5,6 -o distinct -i <(sort -k1,1 -k2,2n ${DCRBED}) > ${DCR2LIFT}
if [ "$FLAG" = true ] && [ -s "$OUT1" ]; then
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping liftOver step as output file already exists and is not empty."
else
# LO
~/Documents/liftOver ${DCR2LIFT} ${CHAINFILE} ${OUT1} ${OUT2}
# assess how many regions were successfully lifted over
  INOUT=$(wc -l output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed | awk '{print $1}')
  DCROUT=$(cut -f 4 output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed | sed 's/;/\t/g' | cut -f 1 | sort | uniq | wc -l | awk '{print $1}')
  INTMISS=$(grep -v '#' output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.unlifted.bed | wc -l | awk '{print $1}')
  DCRinMISS=$(grep -v '#' output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.unlifted.bed | cut -f 4 | sed 's/;/\t/g' | cut -f 1 | sort | uniq | wc -l | awk '{print $1}')
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] After LO there are ${INOUT} intervals for ${DCROUT} DCRs."
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] ${INTMISS} intervals, part of ${DCRinMISS} DCRs were not lifted over "
fi

## 2. Convert to fasta and translate
OUT3="output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa"
OUT4="output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa_translated.fa"
OUT5="output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv"
if [ "$FLAG" = true ] && [ -s "$OUT5" ]; then
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping fasta conversion and translation step as output file already exists and is not empty."
else
  # grab fasta
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Converting lifted over DCR regions to fasta"
  bedtools getfasta -fi data/mm39.fa -bed ${OUT1} -s > ${OUT3}
  # translate (top strand only)
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Translating lifted over DCR regions to amino acids"
  esl-translate -l 1 --watson ${OUT3} > ${OUT4}
  # convert to table so I can do stats | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' # from https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file
  awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${OUT4} | \
  awk 'getline seq {print $0"\t"seq}' | \
  sed 's/>//g;s/source=//g;s/coords=//g;s/\.\./\t/g;s/length=//g;s/frame=//g;s/desc=//g;s/ /\t/g;s/(-)/\t-/g;s/(+)/\t+/g' | \
  awk '{split($2,coord,"[(,),:,-]"); print $1"\t"coord[1]"\t"coord[2]"\t"coord[3]"\t"$3"\t"(coord[3]-coord[2])"\t"$4"\t"$5"\t"($5-$4)"\t"(($5-$4)/(coord[3]-coord[2]))"\t"$6"\t"$7"\t"$8}' | \
  sort -k2,2 -k3,3n -k4,4n > ${OUT5}
fi

## 3. Generate random intervals on mouse exons
# I will not repeat part with EGGIO to id orthologous exons. See README in data/ for details.
# generate random intervals on mouse exons (1-to-1 orthologous non-DCE exons, and with the same region size distro that DCE regions)
OUT6="output/02_Conservation/fake_DCE_intervals.bed"
if [ "$FLAG" = true ] && [ -s "$OUT6" ]; then
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping random intervals generation step as output file already exists and is not empty."
else
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Generating random intervals on mouse exons (1-to-1 orthologous non-DCE)"
  python3.11 scripts/02_Conservation/random_CE_regions.py -n 100000 -e  data/mm39.exons.nonOrthoToHumanDCRintervals.bed -d output/01_Description/filtered_DCR_overlapping_regions.bed
fi

## 4. Convert random intervals to fasta and translate
OUT8="output/02_Conservation/fake_DCE_intervals.fa_translated.tsv"
if [ "$FLAG" = true ] && [ -s "$OUT8" ]; then
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping fasta conversion and translation step for random intervals as output file already exists and is not empty."
else
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Converting random intervals to fasta and translating"
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
fi

## 5. Shuffle mouse exons for random model
OUT9="output/02_Conservation/fake_mm39_exons.fa"
if [ "$FLAG" = true ] && [ -s "$OUT9" ]; then   
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping fake mouse exons generation step as output file already exists and is not empty."
else
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Generating fake mouse exons" 
  python3.11 scripts/02_Conservation/generate_fake_mouse_exons.py
fi
# translate and format
OUT10="output/02_Conservation/fake_mm39_exons.fasta_translated.fa_linearized_renamedWithCoord.tsv"
if [ "$FLAG" = true ] && [ -s "$OUT10" ]; then
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping translation step for shuffled mouse exons as output file already exists and is not empty."
else
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Translting shuffled mouse exons "
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
fi  

## 6. Create the mouse conservation plot
Fig2A="output/02_Conservation/Conservation_mice.pdf"
if [ "$FLAG" = true ] && [ -s "$Fig2A" ]; then
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping mouse conservation plot generation step as output file already exists and is not empty."
else
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Generating mouse conservation plot"
  # create correspondance between fake mouse exons and DCRs
  bedtools intersect -wao -a <(cut -f 2- output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv) -b output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed | awk '$14 != "-1" {print "DCE:"$1":"$2":"$3"\t"$9"\t"$16}' | sed 's/DCEid=//g;s/;nFrames=/\t/g' | cut -f 1-3 > output/02_Conservation/mouse_exon_to_dce_id.txt
  # plot
  Rscript scripts/02_Conservation/Mice_conservation.r -d output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv -n output/02_Conservation/fake_DCE_intervals.fa_translated.tsv -e output/02_Conservation/fake_mm39_exons.fasta_translated.fa_linearized_renamedWithCoord.tsv
fi

#### FIGURE 2 B ####
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Running PhyloP workflow..."
GENECODE_CDS="data/hg38.GENECODE.V48.bed"
DCRregions="output/01_Description/filtered_DCR_overlapping_regions.bed"
output_dir="output/02_Conservation"
ENSG2ENST="data/correspondance_ENSG_ENST_v114.hg38.14"
PPOUT1="output/02_Conservation/DCR_exons_DCR_genes_UNIQ.bed"
PPOUT2="output/02_Conservation/nonDCR_exon_DCR_gene_UNIQ.bed"
PPOUT3="output/02_Conservation/nonDCR_exon_nonDCR_gene_UNIQ.bed"

# file prep for PhyloP
if [ "$FLAG" = true ] && [ -s "$PPOUT1" ] && [ -s "$PPOUT2" ] && [ -s "$PPOUT3" ]; then
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping PhyloP track preparation as output files already exist and are not empty."
  # still print stats
  DCRin=$(cut -f 4 ${DCRregions} | sed 's/DCEid=//g;s/;nFrames/\t/g' | cut -f 1 | sort | uniq | wc -l | awk '{print $1}')
  DCR_intersect=$(cut -f 5 ${output_dir}/DCE_CDS_intersect.bed | sort | uniq | wc -l | awk '{print $1}')
  difference=$((DCRin - DCR_intersect))
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Number of DCR in the intersected file: ${DCR_intersect} (difference with filtered = ${difference})"
  DCR_size_filtered=$(cut -f 5 ${output_dir}/DCE_CDS_intersect_filtered.bed | sort | uniq | wc -l | awk '{print $1}')
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Number of DCR in the size filtered file: ${DCR_size_filtered} (difference with filtered = $((DCRin - DCR_size_filtered)))" 
  DCRwGname=$(join -12 -21 <(sort -k2,2 ${ENSG2ENST}) <(sort -k1,1 ${output_dir}/DCE_CDS_intersect_filtered_with_transcriptID.bed) | sed 's/ /\t/g' | cut -f 7 | sort | uniq | wc -l | awk '{print $1}')
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Number of DCR with gene annotation: ${DCRwGname} (difference with filtered = $((DCRin - DCRwGname)))"
  beforePPcount=$(cut -f 5 output/02_Conservation/DCR_exons_DCR_genes_UNIQ.bed | sort | uniq | wc -l | awk '{print $1}')
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Number of DCR exons before PhyloP: ${beforePPcount}"
else
  DCRin=$(cut -f 4 ${DCRregions} | sed 's/DCEid=//g;s/;nFrames/\t/g' | cut -f 1 | sort | uniq | wc -l | awk '{print $1}')
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] DCR in the input file: ${DCRin}"
  ## 1. intersect DCE regions with the CDS annotation file #
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Intersecting DCE regions with the CDS annotation file"
  bedtools intersect -wao -b ${DCRregions} -a ${GENECODE_CDS} | awk '$8 != "-1" && $6 == $12' | cut -f 1-4,10,13 |\
  sed 's/DCEid=//g;s/;nFrames=/\t/g' > ${output_dir}/DCE_CDS_intersect.bed
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
  #echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Number of DCR in the size filtered file: ${DCR_size_filtered} (difference with filtered = $((DCR_in - DCR_size_filtered)))" >> ${log_file}

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
fi

# Addition 09-26-25: remove DCR with ambiguous isoforms mappings (wrong start or end due to small first or last exon)
grep -vf <(awk '{print "\\t"$0"\\t"}' DCR_to_remove_due_to_isoforms_problematic_end_regions.txt) output/02_Conservation/DCR_exons_DCR_genes_UNIQ.bed > output/02_Conservation/DCR_exons_DCR_genes_UNIQ.bed_fixed
mv output/02_Conservation/DCR_exons_DCR_genes_UNIQ.bed_fixed output/02_Conservation/DCR_exons_DCR_genes_UNIQ.bed

# 5. Run PhyloP
PPOUT4_1="output/02_Conservation/DCEvsNonDCE.meanPhyloP100.out"
#PPOUT4_2="output/02_Conservation/DCEvsNonDCE.meanPhyloP30.out"
# PPOUT5="output/02_Conservation/meanPP_per_DCE.tsv"
if [ "$FLAG" = true ] && [ -s "$PPOUT4_1" ] && [ -s "$PPOUT4_2" ]; then
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping PhyloP analysis as output files already exist and are not empty."
else
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Running PhyloP per exon..."
  echo "[info $(date '+%Y-%m-%d %H:%M:%S')] PP 100-ways..."
  python3.11 scripts/02_Conservation/getPhyloPmeanScore.py data/hg38.phyloP100way.bw > output/02_Conservation/DCEvsNonDCE.meanPhyloP100.out
  # echo "[info $(date '+%Y-%m-%d %H:%M:%S')] PP 30-ways..."
  # python3.11 scripts/02_Conservation/getPhyloPmeanScore.py data/hg38.phyloP30way.bw > output/02_Conservation/DCEvsNonDCE.meanPhyloP30.out
fi

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
  bedtools intersect -wao -a <(cut -f 2- output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv) -b output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed | awk '$14 != "-1" {print "DCE:"$1":"$2":"$3"\t"$9"\t"$16}' | sed 's/DCEid=//g;s/;nFrames=/\t/g' | cut -f 1-3 > output/02_Conservation/mouse_exon_to_dce_id.txt 
  # plot
  Rscript scripts/02_Conservation/PPoutliers.r
fi

echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Conservation analysis completed"