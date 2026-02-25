#!/usr/bin/env bash

#####################################
# 01 Descriptions of DCR | Figure 3 #
#####################################

# Parse arguments
usage() {
    echo "Usage: $0 [-s] -i <path/to/homo_sapiens.csv>"
    echo "  -s    Skip command if output file exists and is not empty (default = false)"
    echo "  -i    Input file (required)"
    echo "  -I    Sequence identity to retain mapping (default = 80)"
    exit 1
}

FLAG=false
INPUT_FILE=""

while getopts ":si:I:" opt; do
  case $opt in
    s)
      FLAG=true
      ;;
    i)
      INPUT_FILE="$OPTARG"
      ;;
    I)
      IDTY="$OPTARG"
      ;;
    \?)
      usage
      ;;
  esac
done

# set defaults
if [ -z "$IDTY" ]; then
    IDTY=80
    echo "[warning $(date '+%Y-%m-%d %H:%M:%S')] No identity (-I) value specified, defaulting to 80% identity"
fi

if [ -z "$INPUT_FILE" ]; then
    echo "[warning $(date '+%Y-%m-%d %H:%M:%S')] No input file specified"
    usage
fi

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') Error: Input file '$INPUT_FILE' does not exist."
    exit 1
fi

# create output and log dirs of do not exist
mkdir -p output/01_Description
mkdir -p logs/01_Description

############## SCRIPT STARTS HERE ##############

## First we generate a DCR track (.bed) with all candidate regions from MIRAGE 2 ##
# we need a .genome file from the reference genome 
samtools faidx data/hg38.fasta
awk '{print $1"\t"$2}' data/hg38.fasta.fai > data/hg38.genome

# Create DCR regions bed ########################
INDIR=data/Extracted-DCEs_06202025
BEDOUT="output/01_Description/raw_DCR_overlapping_regions.bed"
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Extracting exon phases from mapping data..."
cat data/Species-MSAs/homo_sapiens/mappings/*.out |\
 awk '/Sequence ID:/ {print $0; next} /Aminos/{getline coord; n=split(coord, Scoord, ","); print $0"\t"Scoord[1]"\t"Scoord[n]}' |\
  awk '{if($0~/\[revcomp\]/) {print $0"\t"$NF} else {print $0"\t"$(NF-1)}}' | sed 's/\.\./-/g' |\
   awk '/Sequence ID/ {sequenceID = $3; next} /\*/{aminoPos = $NF; split($4, coord, "[:-]"); print sequenceID".."coord[2]".."coord[3]"\t"aminoPos%3}' > output/01_Description/raw_DCR_exons_with_phase.tsv
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Generating bed track for DCR. This can take a while..."
touch logs/01_Description/DCR_regions_bed.log
ls ${INDIR}/homo_sapiens/* |\
parallel scripts/01_Description/extract_DCE_regions_v2_fast.sh 2>> logs/01_Description/DCR_regions_bed.log |\
grep 'DCRid=' 2>> logs/01_Description/DCR_regions_bed.log > "$BEDOUT"

# check count of DCR at this point after some fail (due to chromosome assignment error)
# cut -f 4 output/01_Description/raw_DCR_overlapping_regions.bed | sed 's/DCRid=//g;s/;/\t/g' | cut -f 1 | sort | uniq | wc -l
#    2021

# Self-overlap to id nested DCRs
bedtools intersect -wao -a output/01_Description/raw_DCR_overlapping_regions.bed -b output/01_Description/raw_DCR_overlapping_regions.bed | awk '$4 != $10' > output/01_Description/raw_DCR_overlapping_regions.self_overlap.out
# add total DCR length per DCR
## calculate Total DCR length
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $3-$2}' output/01_Description/raw_DCR_overlapping_regions.bed | bedtools groupby -i - -g 4 -c 5 -o sum > output/01_Description/raw_DCR_total_length.txt
## group OL file intervals per DCR pairs (i.e. if two DCR overlap at two regions, group and sum the total overlap to compare with total DCR)
awk '{print $4"\t"$10"\t"$4"_"$10"\t"$13}' output/01_Description/raw_DCR_overlapping_regions.self_overlap.out | bedtools groupby -i - -g 3 -c 4 -o sum > output/01_Description/rawDCRpair_totalOL.txt
## add total length per DCR
join -11 -21 <(sed 's/_/\t/g' output/01_Description/rawDCRpair_totalOL.txt | sort -k1,1) <(sort -k1,1 output/01_Description/raw_DCR_total_length.txt) | awk '{if($4 > $3){print $0"\tKEEP"} else if($4 == $3){print $0"\tREMOVE"} else {print $0"\tCHECK"}}' > output/01_Description/raw_DCR_post_overlap_adjudication.txt
## make list
grep 'REMOVE' output/01_Description/raw_DCR_post_overlap_adjudication.txt | sed 's/;/\t/g;s/DCRid=//g' | cut -f 1 | sort | uniq | awk '{print $0"\tsub-DCR"}' > output/01_Description/TO_REMOVE.sub-DCR.txt

# we have 2 genes / 5 DCR with "CHECK" --> these will go away with later filters

######### back to OG filters #################

# 0. list Intron-retention events DCR to remove
grep -vwf <(awk 'NR > 1 {gsub(",",""); print $1}' data/Extracted-DCEs_06202025/homo_sapiens.IREs-removed.csv | sort | uniq) <(awk 'NR > 1 {gsub(",",""); print $1}' data/Extracted-DCEs_06202025/homo_sapiens.csv | sort | uniq) | sort | uniq | awk '{print $1"\tIRE"}' > output/01_Description/TO_REMOVE.IRE.txt

# 1. convert to TSV and add generic connection information
#INPUT_FILE=data/Extracted-DCEs_06202025/homo_sapiens.IREs-removed.csv
sed 's/, /\t/g' ${INPUT_FILE} |\
 awk '{if(length($5) == 32 || $5 ~ /Insufficient/) {print $0"\tConnected"} else {print $0"\t"$5}}' |\
 awk -F "\t" '{if(length($6) == 32 || $6 ~ /Insufficient/) {print $0"\tConnected"} else {print $0"\t"$6}}' |\
 awk -F "\t" '{gsub(/Insufficient Nucls \(/, "Insufficient_Nucls_\("); print $0}' > output/01_Description/homo_sapiens.tsv

# 2. list DCR with hits with < 80% identity to remove
awk -F "\t" -v div=80 'NR > 1{sub(/%/, "", $7); if($7+0 < div) {print $1"\t"$4}}' output/01_Description/homo_sapiens.tsv | cut -f 1 | sort -n | uniq | awk '{print $1"\tIDTY_BELOW_80"}' > output/01_Description/TO_REMOVE.IDTY80.txt

# 3. list DCR that do not have 2 frames in the table (mapping issues / can overlap with identity / other extraction issues)
awk 'NR > 1' data/Extracted-DCEs_06202025/homo_sapiens.csv | awk '{print $1"\t"$4}' | sed 's/,//g' | sort | uniq | cut -f 1 | sort | uniq -c | awk '$1 < 2 {print $2"\tONLY_ONE_FRAME"}' > output/01_Description/TO_REMOVE.only1frame.txt

# 4. filter to retain DCR at non-ambiguous genes only
cat data/Species-MSAs/homo_sapiens/mappings/*.out |\
awk '/Sequence ID:/ {print $0; next} /Aminos/{getline coord; n=split(coord, Scoord, ","); print $0"\t"Scoord[1]}' |\
sed 's/\.\./-/g' | awk '/Sequence ID/ {sequenceID = $3; next} /\*/{aminoPos = $NF; split($4, coord, "[:-]"); print coord[1]"\t"coord[2]"\t"coord[3]"\t"sequenceID"\t"aminoPos}' |\
awk '/revcomp/ {strand = "-"; sub(/\[revcomp\]/,""); print $1"\t"$3-1"\t"$2"\t"$4"\t1000\t"strand"\t"$2-$NF; next} !/revcomp/ {chrom=$1; strand = "+"; print $1"\t"$2-1"\t"$3"\t"$4"\t1000\t"strand"\t"$NF-$2}' > output/01_Description/UniprotKB_isoforms_mapping.bed


# check if each gene has consistent chromosome mapping
awk '{split($4, gene, "[\\|/]"); print $1"\t"gene[2]}' output/01_Description/UniprotKB_isoforms_mapping.bed | sort | uniq | cut -f 2 | sort | uniq -d > output/01_Description/isoTrack_ambiguous_gene_mapping.txt
misses=$(wc -l output/01_Description/isoTrack_ambiguous_gene_mapping.txt | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] ${misses} genes have ambiguous chromosome mapping"
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Removing associated DCRS..."
#grep -vwf output/01_Description/isoTrack_ambiguous_gene_mapping.txt output/01_Description/homo_sapiens_filter.${IDTY}pcnt.tsv > output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv
# store retained DCR
#RETAINED=$(awk 'NR > 1 {print $1}' output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv | sort | uniq)
grep -wf output/01_Description/isoTrack_ambiguous_gene_mapping.txt output/01_Description/homo_sapiens.tsv | cut -f 1 | sort -n | uniq | awk '{print $1"\tAMBIGUOUS_GENE_MAPPING"}' > output/01_Description/TO_REMOVE.AMB_GENE_MAP.txt

# convert the bed track for the NMD script
awk '{split($4, iso, "[\\|/]"); print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"iso[2]}' output/01_Description/UniprotKB_isoforms_mapping.bed > output/01_Description/UniprotKB_isoforms_mapping_forNMD.bed
# remove ambiguous genes
grep -vwf output/01_Description/isoTrack_ambiguous_gene_mapping.txt output/01_Description/UniprotKB_isoforms_mapping_forNMD.bed > output/01_Description/UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed

python3.11 scripts/01_Description/convert_bed_to_ExtendedBed.py -i output/01_Description/UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed |\
 awk '{$9 = "255,0,0"; print $0}' | sed 's/ /\t/g' > output/01_Description/Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguous.extBed.bed

# Addition 09-26-25: id DCR with ambiguous isoforms mappings (wrong start or end due to small first or last exon)
# get the problematic isoforms mapping
bedtools getfasta -s -name -fi data/hg38.fa -bed <(cut -f 1-6 output/01_Description/Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguous.extBed.bed | awk '{if($6 == "+"){print $1"\t"$2"\t"$2+3"\t"$4"|start\t\.\t"$6; print $1"\t"$3"\t"$3+3"\t"$4"|end\t\.\t"$6} else {print $1"\t"$2-3"\t"$2"\t"$4"|end\t\.\t"$6; print $1"\t"$3-3"\t"$3"\t"$4"|start\t\.\t"$6}}') | awk 'getline seq {print $0"\t"toupper(seq)}' | awk '/start/{if($2 == "ATG"){print $0"\tPASS"} else {print $0"\tFAIL"}; next} /end/{if($2 == "TAA" || $2 == "TAG" || $2 == "TGA"){print $0"\tPASS"} else {print $0"\tFAIL"}}' | grep 'FAIL' | awk '{split($1, name, "::"); print name[2]"\t"name[1]"\t"$2"\t"$3}' | sed 's/(/\t/g;s/)//g' | awk '{split($1,chrpos,":"); print chrpos[1]"\t"chrpos[2]"\t"$3"\t"$2"\t"$4"\t"$5}' | awk '{split($2,pos,"-"); print $1"\t"pos[1]"\t"pos[2]"\t"$3"|"$5"|"$6"\t\.\t"$4}' > output/01_Description/problematic_end_regions.bed
# assign DCR to isoform
./scripts/01_Description/DCR_isoform_group_to_frame.sh data/Extracted-DCEs_06202025/homo_sapiens > output/01_Description/isoforms_to_DCR_metadata.tsv

# find the DCRs to remove
join -11 -21 <(cut -f 4 output/01_Description/problematic_end_regions.bed | sort | uniq | sed 's/|end/\t|end/g;s/|start/\t|start/g;s/>//g' | sort -k1,1) <(sort -k1,1 output/01_Description/isoforms_to_DCR_metadata.tsv) | awk '{print $NF"\tMAPPING_ENDS_ERROR"}' | sort | uniq > output/01_Description/TO_REMOVE.MAPPING_errors.txt

# # remove them
# cat <(head -n1 output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv) <(grep -vf <(awk '{print "^"$0"\\b"}' DCR_to_remove_due_to_isoforms_problematic_end_regions.txt) <(awk 'NR > 1' output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv)) > output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes_fixed.tsv
# mv output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes_fixed.tsv output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv

########### END OF ID DCR TO FILTER OUT #################

# we now apply the filters #

# first to the isoform track: we only remove those with problematic end region as we will remove the DCRs involving them
grep -vFf <(cut -f 4 output/01_Description/problematic_end_regions.bed | sort | uniq | sed 's/|end/\t|end/g;s/|start/\t|start/g;s/>//g' | cut -f 1 | awk 'length($1)>15') output/01_Description/Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguous.extBed.bed > output/01_Description/FINAL_Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguousGenes_noAmbiguousEnds.extBed.bed

# next to the DCR track
# first compile all DCR to remove
cat output/01_Description/TO_REMOVE.* | sort -k1,1 -k2,2 | sort -u -k1,1 | cut -f 2 | sort | uniq -c
 #  17 AMBIGUOUS_GENE_MAPPING
 #  38 IDTY_BELOW_80
 # 237 IRE
 # 229 MAPPING_ENDS_ERROR
 #  24 ONLY_ONE_FRAME
 #  62 sub-DCR
cat output/01_Description/TO_REMOVE.* | sort -k1,1 -k2,2 | sort -u -k1,1 | cut -f 1 > output/01_Description/TO_REMOVE.all.txt

# make filtered track
grep -vf <(awk '{print "DCRid="$1";"}' output/01_Description/TO_REMOVE.all.txt) output/01_Description/raw_DCR_overlapping_regions.bed | sed 's/DCRid=//g;s/;nFrames=/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7}' > output/01_Description/DCR_overlapping_regions_filtered.bed
# we have 1432

# get geneID and DCRid from the filtered DCRs
sed 's/|/\t/g;s/\//\t/g' output/01_Description/isoforms_to_DCR_metadata.tsv | cut -f 2,7 | sort | uniq | awk 'NR > 1' > output/01_Description/geneID_DCRid.txt
join -14 -22 <(sort -k4,4 output/01_Description/DCR_overlapping_regions_filtered.bed) <(sort -k2,2 output/01_Description/geneID_DCRid.txt) | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$7"\t"$8"\t"$5}' > output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed
# check which 2 are gone --> two ervs that were not in the raw table but still appeared in the .out so made their way to the raw track.

# filter and format the table for analyses
header=$(echo -e "Index\tGene\tGroups\tFrames\tLeftWindow\tRightWindow\tTransAliPctID\tLW\tRW")
cat <(echo "${header}") <(grep -wf <(cut -f 4 output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed | awk '{print "^"$0}') <(awk 'NR > 1' output/01_Description/homo_sapiens.tsv | sort -k1,1n)) > output/01_Description/filtered_homo_sapiens_DCR.csv

# trace filtering steps
grep -vwf <(awk 'NR > 1' output/01_Description/filtered_homo_sapiens_DCR.csv | cut -f 1 | sort | uniq) <(awk 'NR > 1' output/01_Description/homo_sapiens.tsv | cut -f 1 | sort | uniq) > output/01_Description/DCR_filtered_out_from_row_table.txt
cat output/01_Description/TO_REMOVE.AMB_GENE_MAP.txt output/01_Description/TO_REMOVE.IDTY80.txt output/01_Description/TO_REMOVE.IRE.txt output/01_Description/TO_REMOVE.MAPPING_errors.txt output/01_Description/TO_REMOVE.only1frame.txt output/01_Description/TO_REMOVE.sub-DCR.txt | sort -k1,1 -k2,2 | sort -u -k1,1 > output/01_Description/TO_REMOVE.all_with_reason.txt

join -a1 -11 -21 <(sort -k1,1 output/01_Description/DCR_filtered_out_from_row_table.txt) <(sort -k1,1 output/01_Description/TO_REMOVE.all_with_reason.txt) | awk '{print $2}' | sort | uniq -c
# From a total of 1800 DCR in the raw table +237 IRE = 2037 DCR candidates
#   237 IRE
#   17 AMBIGUOUS_GENE_MAPPING
#   38 IDTY_BELOW_80
#  229 MAPPING_ENDS_ERROR
#   24 ONLY_ONE_FRAME <---- this category highlight unaccounted cases in the filtering, which are then manually inspected and described below:
#   62 sub-DCR

# SLATED FOR REMOVAL:
# -------------------
# 71: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 309: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 537: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 588: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 727: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 1007: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 1310: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 1315: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 1321: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 1521: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 1606: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 1749: only one frame in the .out --> NOT DCR AFTER MANUAL REVIEW
# 206:  sub of 205 that is IRE --> remove
# 231:  continuation of 232 that is IRE --> remove
# 1238: continuation of 1239 that is IRE --> remove
# 1761: sub of 1760 but was not counted due to single frame junction issue --> remove
# 1778: sub of 1777 but was not counted due to single frame junction issue --> remove
# 1923: sub of 1922 but was not counted due to single frame junction issue --> remove
# 845:  double alignment of isoform 1 in two frames --> remove
# 1352: double alignment of isoform 1 in two frames --> remove

# ADJACENT DCR WITH ABSOLUTE FRAME CHANGE BETWEEN WITHIN ISOFORM
# --------------------------------------------------------------
# these two cases can be folded into the upstream DCR as isoforms remained grouped together in spite of using a second absolute frame downstream 
# 1088: continuation of 1087, upstream junctions are identical between frames --> merge into 1087
# 2034: continuation of 2035, upstream junctions are identical between frames --> merge into 2035
# 2001: continuation of 2000, fold in 2000
# these two cases cannot be folded into the same DCRs as there is a crossover between which isoform have distinct peptides between the upstream and downstream DCR
# 1171: continuation of 1170, upstream junctions are identical between frames --> will be show in "special case paragraph"


# MANUAL EDITION OF output/01_Description/filtered_homo_sapiens_DCR.csv

# 1087    agtrap  1   3   ggctgctgacaacctgGGGCTGCATTGTATTC    CTCCTGGTCCACACTGgtttccttgggtcttc    100.0%  Connected   Connected
# 1087    agtrap  2   3   ggctgctgacaacctgGGGCTGCATTGTATTC    TACCGGGAGCGCGGGGgtttccttgggtcttc    100.0%  Connected   Connected
# 1087    agtrap  3   2   gcacagggcaacagctGGGCTGCATTGTATTC    C-Terminal  100.0%  Connected   C-Terminal
# 1087    agtrap  4   2   gcacagggcaacagctGGGCTGCATTGTATTC    TACCGGGAGCGCGGGGgtttccttgggtcttc    100.0%  Connected   Connected
# 1087    agtrap  5   2   tgctgtgaacctgaagGGGCTGCATTGTATTC    TACCGGGAGCGCGGGGgtttccttgggtcttc    100.0%  Connected   Connected
# > to become
# 1087    agtrap  1   3   ggctgctgacaacctgGGGCTGCATTGTATTC    C-Terminal  100.0%  Connected   C-Terminal
# 1087    agtrap  2   3   ggctgctgacaacctgGGGCTGCATTGTATTC    C-Terminal  100.0%  Connected   C-Terminal
# 1087    agtrap  3   2   gcacagggcaacagctGGGCTGCATTGTATTC    C-Terminal  100.0%  Connected   C-Terminal
# 1087    agtrap  4   2   gcacagggcaacagctGGGCTGCATTGTATTC    TACCGGGAGCGCGGGGgtttccttgggtcttc    100.0%  Connected   Connected
# 1087    agtrap  5   2   tgctgtgaacctgaagGGGCTGCATTGTATTC    TACCGGGAGCGCGGGGgtttccttgggtcttc    100.0%  Connected   Connected

# awk '$4 == "DCRid=1088;nFrames=2"' output/01_Description/raw_DCR_overlapping_regions.bed | sed 's/DCRid=1088;nFrames=2/DCRid=1087;nFrames=2/g'
# chr1  11750076    11750101    DCRid=1087;nFrames=2    .   +

echo -e "chr1\t11750076\t11750101\t1087\tagtrap\t." >> output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed

# MANUAL EDITION OF output/01_Description/filtered_homo_sapiens_DCR.csv

# 2035, icam4, 1, 1, aggatcaccgcctacaAACCGCCCCACAGCGT, ATTACACTGATGCTCGcttggagccccgcgcc, 100.0%
# 2035, icam4, 2, 1, aggatcaccgcctacaAACCGCCCCACAGCGT, C-Terminal, 100.0%
# 2035, icam4, 3, 3, aggatcaccgcctacaGTGTTCCCGGTGGGCT, ATTACACTGATGCTCGcttggagccccgcgcc, 98.6%
# > to become
# 2035, icam4, 1, 1, aggatcaccgcctacaAACCGCCCCACAGCGT, C-Terminal, 100.0%
# 2035, icam4, 2, 1, aggatcaccgcctacaAACCGCCCCACAGCGT, C-Terminal, 100.0%
# 2035, icam4, 3, 3, aggatcaccgcctacaGTGTTCCCGGTGGGCT, C-Terminal, 98.6%

# awk '$4 == "DCRid=2035;nFrames=2" || $4 == "DCRid=2034;nFrames=2"' output/01_Description/raw_DCR_overlapping_regions.bed | sed 's/DCRid=2034;nFrames=2/DCRid=2035;nFrames=2/g'
# chr19   10287985    10288101    DCRid=2035;nFrames=2    .   +
# chr19   10287612    10287838    DCRid=2035;nFrames=2    .   +
echo -e "chr19\t10287612\t10287838\t2035\ticam4\t." >> output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed

# awk '$4 == "DCRid=2001;nFrames=2" || $4 == "DCRid=2000;nFrames=2"' output/01_Description/raw_DCR_overlapping_regions.bed | sed 's/DCRid=2001;nFrames=2/DCRid=2000;nFrames=2/g'
# chr6    25913319    25913451    DCRid=2000;nFrames=2    .   -
# chr6    25914579    25914670    DCRid=2000;nFrames=2    .   -
echo -e "chr6\t25914579\t25914670\t2000\tslc17a2\t." >> output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed


# here we perform a blank run of the architecture classification to check for other false positives
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' output/01_Description/filtered_homo_sapiens_DCR.csv > output/01_Description/filtered_homo_sapiens_DCR.tsv
python3.11 scripts/01_Description/categorize_DCR_with_up_down_exons_edges_k16.py --input output/01_Description/filtered_homo_sapiens_DCR.tsv --out-detailed output/01_Description/NEW_DCR_architectures_exon_edges.tsv --out-upset output/01_Description/NEW_DCR_architectures_exon_edges_upset.tsv

# we identify cases that have indistinguishable upstream exon-to-DCR-exon-edge junction between frame: this indicate that it is in continuity with another event

# ADDITIONAL MANUAL REMOVAL / EDITS:
# ----------------------------------

# 1958: continuation of 1957 removed due to incorrect mapping.
# 44,45,43 for il23r: complex pattern with frame switching at 3 DCRs --> all pulled out and added to special cases
# 370, 371 for larp1b: complex pattern with frame switching at 2 DCRs --> both pulled out and added to special cases
# 540, 541 for fut8: complex pattern with frame switching at 2 DCRs --> both pulled out and added to special cases
# 670, 671 for nrg2: complex pattern with frame switching at 2 DCRs --> both pulled out and added to special cases
# 972, 973 for ubac2: complex pattern with frame switching at 2 DCRs --> both pulled out and added to special cases
# 2025, 2026 for mf170: complex pattern with frame switching at 2 DCRs --> both pulled out and added to special cases


# 304, 305: complex downstream case: isoform 4 switch group between 304 and 305 due to AS --> both pulled out and added to special cases
# 519: isoform 2 is AS into the same junction that the isoforms in the other frame, creating a group switch that did not trigger a DCR call but is equivalent to a complex case. --> pulled and added to special cases
# 624, 625 group switch for isoform 5 between the DCRs --> both pulled out and added to special cases
# 1180, 1181 isoform 4 does a group switch going into a 3rd frame in 1181 --> both pulled out and added to special cases
# 1281, 1282 isoform 2 does a group switch --> both pulled out and added to special cases

# 576: in this case, I think the DCR could have been extended to be a all C-term: isoform 5 (alt frame) use the same AS than canonical frame but make a premature end in short DCR region while the whole DCR can still be quantified from the 5' junctions 
# --> chr5:33,989,487-33,989,502 == chr5    33989486    33989502    576 (bed record to add)
# 576 amacr   1/4 2   cattgatgcaaatatgGTGGAAGGAACAGCAT    GAGCTGCTGATCAAAGgacttggactaaagtc    100.0%  Connected   Connected
# 576 amacr   2   1   tatttggctttgtcagGTGGAAGGAACAGCAT    GAGCTGCTGATCAAAGatctttgataccatac    100.0%  Connected   Connected
# 576 amacr   3   1   tatttggctttgtcagGTGGAAGGAACAGCAT    GAGCTGCTGATCAAAGgacttggactaaagtc    100.0%  Connected   Connected
# > to become
# 576 amacr   1 2   cattgatgcaaatatgGTGGAAGGAACAGCAT    C-Terminal    100.0%  Connected   C-Terminal
# 576 amacr   4 2   cattgatgcaaatatgGTGGAAGGAACAGCAT    AGAGACTAGCTTTTACcttattactttcaatg    100.0%  Connected   Connected
# 576 amacr   2   1   tatttggctttgtcagGTGGAAGGAACAGCAT    GAGCTGCTGATCAAAGatctttgataccatac    100.0%  Connected   Connected
# 576 amacr   3   1   tatttggctttgtcagGTGGAAGGAACAGCAT    C-Terminal    100.0%  Connected   C-Terminal

# 1033 pax3 to extend C-term iso 5 and 8 --> chr2:222,201,224-222,201,240 == chr2   222201223   222201240   1033    pax3 (bed record to add)
# 1033    pax3    1   2   tggcctctcacctcagGTGCCTTTCATTATCT    C-Terminal  100.0%  Connected   C-Terminal
# 1033    pax3    2   2   tggcctctcacctcagGTGCCTTTCATTATCT    TCAAGCCAGATATCGCttggtttcaaatcctt    100.0%  Connected   Connected
# 1033    pax3    3   1   gggcagtatggacaaaGTGCCTTTCATTATCT    C-Terminal  100.0%  Connected   C-Terminal
# 1033    pax3    4   1   gggcagtatggacaaaGTGCCTTTCATTATCT    TCAAGCCAGATATCGCttggtttcaaatcctt    100.0%  Connected   Connected
# > to become
# 1033    pax3    1,2   2   tggcctctcacctcagGTGCCTTTCATTATCT    C-Terminal  100.0%  Connected   C-Terminal
# 1033    pax3    3,4   1   gggcagtatggacaaaGTGCCTTTCATTATCT    C-Terminal  100.0%  Connected   C-Terminal

# create list of additional to remove

echo -e "71\n309\n537\n588\n727\n1007\n1310\n1315\n1321\n1521\n1606\n1749\n206\n231\n1238\n1761\n1778\n1923\n845\n1352\n1958" | awk '{print $1"\tMANUAL_SELECTION"}' > output/01_Description/TO_REMOVE.manual_selection.txt
echo -e "1171\n1170\n44\n45\n43\n370\n371\n540\n541\n670\n671\n972\n973\n2025\n2026\n304\n305\n519\n624\n625\n1180\n1181\n1281\n1282" | awk '{print $1"\tCOMPLEX_CASES"}' > output/01_Description/TO_REMOVE.complex_cases.txt

### re-filter

# clean

mv output/01_Description/TO_REMOVE.all_with_reason.txt output/01_Description/_TO_REMOVE.all_with_reason.txt
mv output/01_Description/TO_REMOVE.all_with_reason_allsets.txt output/01_Description/_TO_REMOVE.all_with_reason_allsets.txt
mv output/01_Description/TO_REMOVE.all.txt output/01_Description/_TO_REMOVE.all.txt
mv output/01_Description/filtered_homo_sapiens_DCR.tsv output/01_Description/_filtered_homo_sapiens_DCR.tsv
mv output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed output/01_Description/_filtered_DCR_overlapping_regions_with_gene_and_dceID.bed

cat output/01_Description/TO_REMOVE.* | sort -k1,1 -k2,2  > output/01_Description/TO_REMOVE.all_with_reason_allsets.txt
cat output/01_Description/TO_REMOVE.* | sort -k1,1 -k2,2 | sort -u -k1,1 > output/01_Description/TO_REMOVE.all_with_reason.txt
cat output/01_Description/TO_REMOVE.* | sort -k1,1 -k2,2 | sort -u -k1,1 | cut -f 1 > output/01_Description/TO_REMOVE.all.txt

# make filtered track
grep -vf <(awk '{print "DCRid="$1";"}' output/01_Description/TO_REMOVE.all.txt) output/01_Description/raw_DCR_overlapping_regions.bed | sed 's/DCRid=//g;s/;nFrames=/\t/g' > output/01_Description/DCR_overlapping_regions_filtered.bed
# we have 1409

# get geneID and DCRid from the filtered DCRs
join -14 -22 <(sort -k4,4 output/01_Description/DCR_overlapping_regions_filtered.bed) <(sort -k2,2 output/01_Description/geneID_DCRid.txt) | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$7"\t"$8"\t"$5}' > output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed
# we have 1407 (1296 genes) --> this is normal, removed the 2 ervs that are not DCR and not in the table

# filter and format the table for analyses
header=$(echo -e "Index\tGene\tGroups\tFrames\tLeftWindow\tRightWindow\tTransAliPctID\tLW\tRW")
cat <(echo "${header}") <(grep -wf <(cut -f 4 output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed | awk '{print "^"$0}') <(awk 'NR > 1' output/01_Description/homo_sapiens.tsv | sort -k1,1n)) > output/01_Description/filtered_homo_sapiens_DCR.csv
# 1407 --> OK

# MANUAL EDITIONS OF output/01_Description/filtered_homo_sapiens_DCR.csv

# 1087    agtrap  1   3   ggctgctgacaacctgGGGCTGCATTGTATTC    CTCCTGGTCCACACTGgtttccttgggtcttc    100.0%  Connected   Connected
# 1087    agtrap  2   3   ggctgctgacaacctgGGGCTGCATTGTATTC    TACCGGGAGCGCGGGGgtttccttgggtcttc    100.0%  Connected   Connected
# 1087    agtrap  3   2   gcacagggcaacagctGGGCTGCATTGTATTC    C-Terminal  100.0%  Connected   C-Terminal
# 1087    agtrap  4   2   gcacagggcaacagctGGGCTGCATTGTATTC    TACCGGGAGCGCGGGGgtttccttgggtcttc    100.0%  Connected   Connected
# 1087    agtrap  5   2   tgctgtgaacctgaagGGGCTGCATTGTATTC    TACCGGGAGCGCGGGGgtttccttgggtcttc    100.0%  Connected   Connected
# > to become
# 1087    agtrap  1/2   3   ggctgctgacaacctgGGGCTGCATTGTATTC    C-Terminal  100.0%  Connected   C-Terminal
# 1087    agtrap  3   2   gcacagggcaacagctGGGCTGCATTGTATTC    C-Terminal  100.0%  Connected   C-Terminal
# 1087    agtrap  4   2   gcacagggcaacagctGGGCTGCATTGTATTC    TACCGGGAGCGCGGGGgtttccttgggtcttc    100.0%  Connected   Connected
# 1087    agtrap  5   2   tgctgtgaacctgaagGGGCTGCATTGTATTC    TACCGGGAGCGCGGGGgtttccttgggtcttc    100.0%  Connected   Connected
# ------------------------------------------------------------------------------------------------------------------------------
# 2035, icam4, 1, 1, aggatcaccgcctacaAACCGCCCCACAGCGT, ATTACACTGATGCTCGcttggagccccgcgcc, 100.0%
# 2035, icam4, 2, 1, aggatcaccgcctacaAACCGCCCCACAGCGT, C-Terminal, 100.0%
# 2035, icam4, 3, 3, aggatcaccgcctacaGTGTTCCCGGTGGGCT, ATTACACTGATGCTCGcttggagccccgcgcc, 98.6%
# > to become
# 2035, icam4, 1/2 1, aggatcaccgcctacaAACCGCCCCACAGCGT, C-Terminal, 100.0%
# 2035, icam4, 3, 3, aggatcaccgcctacaGTGTTCCCGGTGGGCT, C-Terminal, 98.6%
# ------------------------------------------------------------------------------------------------------------------------------
# 576 amacr   1/4 2   cattgatgcaaatatgGTGGAAGGAACAGCAT    GAGCTGCTGATCAAAGgacttggactaaagtc    100.0%  Connected   Connected
# 576 amacr   2   1   tatttggctttgtcagGTGGAAGGAACAGCAT    GAGCTGCTGATCAAAGatctttgataccatac    100.0%  Connected   Connected
# 576 amacr   3   1   tatttggctttgtcagGTGGAAGGAACAGCAT    GAGCTGCTGATCAAAGgacttggactaaagtc    100.0%  Connected   Connected
# > to become
# 576 amacr   1 2   cattgatgcaaatatgGTGGAAGGAACAGCAT    C-Terminal    100.0%  Connected   C-Terminal
# 576 amacr   4 2   cattgatgcaaatatgGTGGAAGGAACAGCAT    AGAGACTAGCTTTTACcttattactttcaatg    100.0%  Connected   Connected
# 576 amacr   2   1   tatttggctttgtcagGTGGAAGGAACAGCAT    GAGCTGCTGATCAAAGatctttgataccatac    100.0%  Connected   Connected
# 576 amacr   3   1   tatttggctttgtcagGTGGAAGGAACAGCAT    C-Terminal    100.0%  Connected   C-Terminal
# ------------------------------------------------------------------------------------------------------------------------------
# 1033    pax3    1   2   tggcctctcacctcagGTGCCTTTCATTATCT    C-Terminal  100.0%  Connected   C-Terminal
# 1033    pax3    2   2   tggcctctcacctcagGTGCCTTTCATTATCT    TCAAGCCAGATATCGCttggtttcaaatcctt    100.0%  Connected   Connected
# 1033    pax3    3   1   gggcagtatggacaaaGTGCCTTTCATTATCT    C-Terminal  100.0%  Connected   C-Terminal
# 1033    pax3    4   1   gggcagtatggacaaaGTGCCTTTCATTATCT    TCAAGCCAGATATCGCttggtttcaaatcctt    100.0%  Connected   Connected
# > to become
# 1033    pax3    1/2   2   tggcctctcacctcagGTGCCTTTCATTATCT    C-Terminal  100.0%  Connected   C-Terminal
# 1033    pax3    3/4   1   gggcagtatggacaaaGTGCCTTTCATTATCT    C-Terminal  100.0%  Connected   C-Terminal
# ------------------------------------------------------------------------------------------------------------------------------
# 2000    slc17a2 1   3   tagatatcgcccccagATATGCAAGTTTCCTC    C-Terminal  100.0%  Connected   C-Terminal
# 2000    slc17a2 2   1   ctcttttcatctcttgATATGCAAGTTTCCTC    ATTCCTCATCAGTCAGgattttgagtctggtt    100.0%  Connected   Connected
# 2000    slc17a2 3   3   tagatatcgcccccagATATGCAAGTTTCCTC    ATTCCTCATCAGTCAGgattttgagtctggtt    100.0%  Connected   Connected
# > to become
# 2000    slc17a2 1/3   3   tagatatcgcccccagATATGCAAGTTTCCTC    C-Terminal  100.0%  Connected   C-Terminal
# 2000    slc17a2 2   1   ctcttttcatctcttgATATGCAAGTTTCCTC    C-Terminal  100.0%  Connected   C-Terminal


# EDITS TO ADD NEW DCR SUB-REGIONS TO THE FILTERED DCR TRACK
echo -e "chr1\t11750076\t11750101\t1087\t+\tagtrap\t2" >> output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed
echo -e "chr19\t10287612\t10287838\t2035\t+\ticam4\t2" >> output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed
echo -e "chr5\t33989486\t33989502\t576\t-\tamacr\t2" >> output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed
echo -e "chr2\t222201223\t222201240\t1033\t-\tpax3\t2" >> output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed
echo -e "chr6\t25914579\t25914670\t2000\t-\tslc17a2\t2" >> output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed

# classify with python script ###########################

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' output/01_Description/filtered_homo_sapiens_DCR.csv > output/01_Description/filtered_homo_sapiens_DCR.tsv
python3.11 scripts/01_Description/categorize_DCR_with_up_down_exons_edges_k16.py --input output/01_Description/filtered_homo_sapiens_DCR.tsv --out-detailed output/01_Description/NEW_DCR_architectures_exon_edges.tsv --out-upset output/01_Description/NEW_DCR_architectures_exon_edges_upset.tsv

# At this stage, we confirm that the 1407 retained DCRs (1296 genes) are independent, i.e. Dual coding at DCR can be directly explained by incoming (or N-terminal) junctions and outgoing junctions will not cause further dual coding between the isoforms still connected (next AS event in same or other isoform may cause another independent DCR). 
# We retain an additional 24 DCR (12 genes) as "special cases" of dependent DCRs with "group switching": i.e. the groups of isoform that share a frame switch between DCRs
# Total reported is thus:
cat output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed output/01_Description/filtered_DCR_complex_cases_with_gene_and_dceID.bed | cut -f 6 | sort | uniq | wc -l
#    1309 genes

# store special cases
cat <(echo "${header}")  <(grep -wf<(cut -f 1 output/01_Description/TO_REMOVE.complex_cases.txt) output/01_Description/homo_sapiens.tsv) > output/01_Description/homo_sapiens_DCR_complex_cases.tsv
grep -f <(awk '{print "DCRid="$1";"}' output/01_Description/TO_REMOVE.complex_cases.txt) output/01_Description/raw_DCR_overlapping_regions.bed | sed 's/DCRid=//g;s/;nFrames=/\t/g' > output/01_Description/DCR_overlapping_complex_cases.bed
join -14 -22 <(sort -k4,4 output/01_Description/DCR_overlapping_complex_cases.bed) <(sort -k2,2 output/01_Description/geneID_DCRid.txt) | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$7"\t"$8"\t"$5}' > output/01_Description/filtered_DCR_complex_cases_with_gene_and_dceID.bed

 #  17 AMBIGUOUS_GENE_MAPPING
 #  24 COMPLEX_CASES
 #  38 IDTY_BELOW_80
 # 237 IRE
 #  21 MANUAL_SELECTION
 # 229 MAPPING_ENDS_ERROR
 #   3 ONLY_ONE_FRAME --> manual selection
 #  61 sub-DCR

################### END OF FILTERING ###################

# NMD susceptibility ###################################

python3.11 scripts/01_Description/isoform_NMD.py
# add uncertainty to per-isoform DCR assessment
awk '{if ($5 == "no_NMD" && $4 <= 50) print $0 "\tnot_enough_evidence"; else if ($5 == "no_NMD" && $4 == "last_exon") print $0 "\tno_NMD"; else print $0 "\tNMD"}' output/01_Description/NMDtable.tsv > output/01_Description/NMDtable_updated.tsv
# Figure 3.C is then produced using the script in Rscript scripts/01_Description/NMD.R
Rscript scripts/01_Description/NMD.R
# NMD                      DCR exons not C-Terminal      Not enough evidence 
# 466                      264                           677 
               
# DCR exons count and distribution #####################

# NOTE/TO DO: for now I let the error print in stdout: the # of exon and distribution will need to be adjusted manually in the R script for aestetic purposes.
Rscript scripts/01_Description/DCR_exon_length_distribution.r output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed
# [1] "total exon count: 1741"
# [1] "min DCR exon length: 20 bp"
# [1] "max DCR exon length: 1878 bp"
# [1] "mean DCR exon length: 95.4228855721393 bp"
# [1] "sd DCR exon length: 92.0820934642561 bp"
# [1] "DCR count: 1407"
#    1    2    3    4    6
# 1131  227   42    6    1

# Splice site analysis #################################

cut -f 4 output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed | sort | uniq > output/01_Description/filtered.DCR.id.txt

# DCR
./scripts/01_Description/extract_DCE_splice_sites.sh data/Extracted-DCEs_06202025/homo_sapiens/ > output/01_Description/SPLICE_SITES_filtered_DCE_splice_sites.bed

# back to normal
bedtools getfasta -s -bed <(grep 'junction=upstream;site=donor' output/01_Description/SPLICE_SITES_filtered_DCE_splice_sites.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t\.\t"$5}'|\
 awk '{print $1"_"$2"_"$3"\t"$0}' | sort -u -k1,1 | cut -f 2-) -fi data/hg38.fa -name > output/01_Description/SPLICE_SITES_UPSTREAM_DONOR_filtered_DCE_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=upstream;site=acceptor' output/01_Description/SPLICE_SITES_filtered_DCE_splice_sites.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t\.\t"$5}'|\
 awk '{print $1"_"$2"_"$3"\t"$0}' | sort -u -k1,1 | cut -f 2-) -fi data/hg38.fa -name > output/01_Description/SPLICE_SITES_UPSTREAM_ACCEPTOR_filtered_DCE_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=downstream;site=donor' output/01_Description/SPLICE_SITES_filtered_DCE_splice_sites.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t\.\t"$5}'|\
 awk '{print $1"_"$2"_"$3"\t"$0}' | sort -u -k1,1 | cut -f 2-) -fi data/hg38.fa -name > output/01_Description/SPLICE_SITES_DOWNSTREAM_DONOR_filtered_DCE_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=downstream;site=acceptor' output/01_Description/SPLICE_SITES_filtered_DCE_splice_sites.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t\.\t"$5}'|\
 awk '{print $1"_"$2"_"$3"\t"$0}' | sort -u -k1,1 | cut -f 2-) -fi data/hg38.fa -name > output/01_Description/SPLICE_SITES_DOWNSTREAM_ACCEPTOR_filtered_DCE_splice_sites.fa
# print logos
/opt/homebrew/bin/weblogo -f output/01_Description/SPLICE_SITES_UPSTREAM_DONOR_filtered_DCE_splice_sites.fa -o output/01_Description/LOGO_DCR_upstream_donor.pdf -A DNA -U probability -F PDF -t "upstream donor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/SPLICE_SITES_UPSTREAM_ACCEPTOR_filtered_DCE_splice_sites.fa -o output/01_Description/LOGO_DCR_upstream_acceptor.pdf -A DNA -U probability -F PDF -t "upstream acceptor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/SPLICE_SITES_DOWNSTREAM_DONOR_filtered_DCE_splice_sites.fa -o output/01_Description/LOGO_DCR_downstream_donor.pdf -A DNA -U probability -F PDF -t "downstream donor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/SPLICE_SITES_DOWNSTREAM_ACCEPTOR_filtered_DCE_splice_sites.fa -o output/01_Description/LOGO_DCR_downstream_acceptor.pdf -A DNA -U probability -F PDF -t "downstream acceptor site" -s large -c classic

# controls
python3.11 scripts/01_Description/extract_control_exons_splice_site.py
bedtools getfasta -s -bed <(grep 'junction=upstream;site=donor' output/01_Description/SPLICE_SITES_control_splice_sites.bed |\
 sort | uniq) -fi data/hg38.fa -name > output/01_Description/SPLICE_SITES_control_UPSTREAM_DONOR_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=upstream;site=acceptor' output/01_Description/SPLICE_SITES_control_splice_sites.bed |\
 sort | uniq) -fi data/hg38.fa -name > output/01_Description/SPLICE_SITES_control_UPSTREAM_ACCEPTOR_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=downstream;site=donor' output/01_Description/SPLICE_SITES_control_splice_sites.bed |\
 sort | uniq) -fi data/hg38.fa -name > output/01_Description/SPLICE_SITES_control_DOWNSTREAM_DONOR_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=downstream;site=acceptor' output/01_Description/SPLICE_SITES_control_splice_sites.bed |\
    sort | uniq) -fi data/hg38.fa -name > output/01_Description/SPLICE_SITES_control_DOWNSTREAM_ACCEPTOR_splice_sites.fa
/opt/homebrew/bin/weblogo -f output/01_Description/SPLICE_SITES_control_UPSTREAM_DONOR_splice_sites.fa -o output/01_Description/LOGO_control_upstream_donor.pdf -A DNA -U probability -F PDF -t "upstream donor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/SPLICE_SITES_control_UPSTREAM_ACCEPTOR_splice_sites.fa -o output/01_Description/LOGO_control_upstream_acceptor.pdf -A DNA -U probability -F PDF -t "upstream acceptor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/SPLICE_SITES_control_DOWNSTREAM_DONOR_splice_sites.fa -o output/01_Description/LOGO_control_downstream_donor.pdf -A DNA -U probability -F PDF -t "downstream donor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/SPLICE_SITES_control_DOWNSTREAM_ACCEPTOR_splice_sites.fa -o output/01_Description/LOGO_control_downstream_acceptor.pdf -A DNA -U probability -F PDF -t "downstream acceptor site" -s large -c classic

### count multi DCR genes

awk 'NR > 1' output/01_Description/filtered_homo_sapiens_DCR.tsv | cut -f 1,2 | sort | uniq | cut -f 2 | sort | uniq -c | awk '$1 > 1 {print $2}' > output/01_Description/filtered_multi_DCR_genes.txt
join -11 -22 <(awk 'NR > 1' output/01_Description/NEW_DCR_architectures_exon_edges_upset.tsv| sort -k1,1) <(sort -k2,2 output/01_Description/geneID_DCRid.txt) | grep 'C-Terminal' | awk '{print $NF}' | sort | uniq | grep -f output/01_Description/filtered_multi_DCR_genes.txt | wc -l
# 92 have C-Terminal
wc -loutput/01_Description/filtered_multi_DCR_genes.txt
# 97 genes have 2+ DCR --> 94.8%


### update sup tables

# table S1

join -12 -21 <(awk 'NR > 1 {print $0}' data/Table_S1.tsv | sort -k2,2) <(cut -f 4 output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed | sort | uniq) | awk '{print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$11}' > output/New_TableS1.tsv
