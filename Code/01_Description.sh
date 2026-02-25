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

# Create DCR regions bed ########################
INDIR=Data/Extracted-DCEs_06202025
BEDOUT="output/01_Description/raw_DCR_overlapping_regions.bed"
#awk 'NR > 1 {print $1}' output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv | sort | uniq > output/01_Description/finalDCR
# create a list of DCR to parallelize


if [ "$FLAG" = true ] && [ -s "$BEDOUT" ]; then
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping DCR regions bed (already exists and not empty)"
else
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Extracting exon phases from mapping data..."
    cat data/Species-MSAs/homo_sapiens/mappings/*.out |\
     awk '/Sequence ID:/ {print $0; next} /Aminos/{getline coord; n=split(coord, Scoord, ","); print $0"\t"Scoord[1]"\t"Scoord[n]}' |\
      awk '{if($0~/\[revcomp\]/) {print $0"\t"$NF} else {print $0"\t"$(NF-1)}}' | sed 's/\.\./-/g' |\
       awk '/Sequence ID/ {sequenceID = $3; next} /\*/{aminoPos = $NF; split($4, coord, "[:-]"); print sequenceID".."coord[2]".."coord[3]"\t"aminoPos%3}' > output/01_Description/DCE_exons_with_phase.tsv
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Generating bed track for DCR. This can take a while..."
    touch logs/01_Description/DCR_regions_bed.log
    ls ${INDIR}/homo_sapiens/* |\
    parallel scripts/01_Description/extract_DCE_regions_v2_fast.sh 2>> logs/01_Description/DCR_regions_bed.log |\
    grep 'DCEid=' 2>> logs/01_Description/DCR_regions_bed.log > "$BEDOUT"
fi
DCRbedCount=$(sed 's/;/\t/g' ${OUT6} | cut -f 4 | sort | uniq | wc -l | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] The DCR bed track has ${DCRbedCount} DCRs (/${finalDCRcnt} retained DCR used in input))"
# # make filterd track
# grep -f <(awk '{print "DCEid="$0";"}' output/01_Description/DCRs_with_2or3_frames.txt) output/01_Description/DCR_overlapping_regions.bed > output/01_Description/DCR_overlapping_regions_filtered.bed
# # count
# DCRbedCountFiltered=$(sed 's/;/\t/g'  output/01_Description/DCR_overlapping_regions_filtered.bed | cut -f 4 | sort | uniq | wc -l | awk '{print $1}')
# echo "[info $(date '+%Y-%m-%d %H:%M:%S')] The filterd DCR bed track has ${DCRbedCountFiltered} DCRs (/${KEPT} in input))"
# save missing DCRs in the between the filtered csv and the filtered track.
grep -wvf <(sed 's/;/\t/g' ${OUT6} | cut -f 4 | sort | uniq) <(awk '{print "DCEid="$0}' <(echo ${finalDCR})) > output/01_Description/DCRs_missing_from_filterd_bed.txt
missCount=$(sort output/01_Description/DCRs_missing_from_filterd_bed.txt| uniq | wc -l | awk '{print $1}')
if [ -s "output/01_Description/DCRs_missing_from_filterd_bed.txt" ]; then
    echo "[warning $(date '+%Y-%m-%d %H:%M:%S')] ${missCount} DCRs are missing between the filtered tsv and the bed track (see output/01_Description/DCRs_missing_from_filterd_bed.txt)."
fi
# add gene name and frame count
OUT7="output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb.bed"
if [ "$FLAG" = true ] && [ -s "$OUT7" ]; then
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb.bed (already exists and not empty)"
else
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Adding gene names and frame counts to DCR bed track"
    # get geneID and DCRid from the filtered DCRs
    awk 'NR > 1 {print $1"\t"$2}' ${OUT2_2} | sort | uniq > output/01_Description/filtered_DCRid_geneID.tsv
    join -14 -21 <(sed 's/DCEid=//g;s/;nFrames=/\t/g' ${OUT6} | sort -k4,4) <(sort -k1,1 output/01_Description/filtered_DCRid_geneID.tsv) | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$7"\t"$8"\t"$5}' > output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb.bed
fi

# Addition 09-26-25: remove DCR with ambiguous isoforms mappings (wrong start or end due to small first or last exon)
grep -vf <(awk '{print "^"$0"\\b"}' DCR_to_remove_due_to_isoforms_problematic_end_regions.txt) <(awk '{print $4"\t"$0}' output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb.bed) > output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb_fixed.bed
mv output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb_fixed.bed output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb.bed
grep -vf <(awk '{print "DCEid="$0";"}' DCR_to_remove_due_to_isoforms_problematic_end_regions.txt) output/01_Description/filtered_DCR_overlapping_regions.bed > output/01_Description/filtered_DCR_overlapping_regions.bed_fixed
mv output/01_Description/filtered_DCR_overlapping_regions.bed_fixed output/01_Description/filtered_DCR_overlapping_regions.bed




















## Count DCR in ##
INDIR=$(dirname "$INPUT_FILE")
DCRIN=$(awk 'NR > 1{print $1}' ${INPUT_FILE}| sort | uniq | wc -l | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] input file has ${DCRIN} DCRs"

## Process the input file ######################

# 1. convert to TSV and add generic connection information
sed 's/, /\t/g' ${INPUT_FILE} |\
 awk '{if(length($5) == 32 || $5 ~ /Insufficient/) {print $0"\tConnected"} else {print $0"\t"$5}}' |\
 awk -F "\t" '{if(length($6) == 32 || $6 ~ /Insufficient/) {print $0"\tConnected"} else {print $0"\t"$6}}' |\
 awk -F "\t" '{gsub(/Insufficient Nucls \(/, "Insufficient_Nucls_\("); print $0}' > output/01_Description/homo_sapiens.tsv

# 2. retain only hits with 80% identity
OUT1="output/01_Description/DCRs_with_2or3_frames.txt"
if [ "$FLAG" = true ] && [ -s "$OUT1" ]; then
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping creation of $OUT1 (already exists and not empty)"
else
    awk -F "\t" -v div=${IDTY} '{sub(/%/, "", $7); if($7+0 >= div) {print $1"\t"$4}}' output/01_Description/homo_sapiens.tsv | sort | uniq | cut -f 1 | uniq -d > "$OUT1"
fi

# 3. filter for mapping identity
OUT2="output/01_Description/homo_sapiens_filter.${IDTY}pcnt.tsv"
if [ "$FLAG" = true ] && [ -s "$OUT2" ]; then
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping creation of $OUT2 (already exists and not empty)"
else
    cat data/input_header.txt <(join -11 -21 output/01_Description/DCRs_with_2or3_frames.txt <(sort -k1,1 output/01_Description/homo_sapiens.tsv) | sed 's/ /\t/g') > "$OUT2"
fi

KEPT=$(join -11 -21 output/01_Description/DCRs_with_2or3_frames.txt <(sort -k1,1 output/01_Description/homo_sapiens.tsv) | sed 's/ /\t/g' | cut -f 1 | sort | uniq | wc -l | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] After filtering @${IDTY}%, ${KEPT}/${DCRIN} DCR are kept"

# 4. filter to retain DCR at non-ambiguous genes only
OUT2_2="output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv"
if [ "$FLAG" = true ] && [ -s "$OUT2_2" ]; then
     echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping creation of $OUT2_2 (already exists and not empty)"
else
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] removing genes with ambiguous chromosome mapping"
    cat data/Species-MSAs/homo_sapiens/mappings/*.out |\
    awk '/Sequence ID:/ {print $0; next} /Aminos/{getline coord; n=split(coord, Scoord, ","); print $0"\t"Scoord[1]}' |\
    sed 's/\.\./-/g' | awk '/Sequence ID/ {sequenceID = $3; next} /\*/{aminoPos = $NF; split($4, coord, "[:-]"); print coord[1]"\t"coord[2]"\t"coord[3]"\t"sequenceID"\t"aminoPos}' |\
    awk '/revcomp/ {strand = "-"; sub(/\[revcomp\]/,""); print $1"\t"$3-1"\t"$2"\t"$4"\t1000\t"strand"\t"$2-$NF; next} !/revcomp/ {chrom=$1; strand = "+"; print $1"\t"$2-1"\t"$3"\t"$4"\t1000\t"strand"\t"$NF-$2}' > output/01_Description/UniprotKB_isoforms_mapping.bed
    # check if each gene has consistent chromosome mapping
    awk '{split($4, gene, "[\\|/]"); print $1"\t"gene[2]}' output/01_Description/UniprotKB_isoforms_mapping.bed | sort | uniq | cut -f 2 | sort | uniq -d > output/01_Description/isoTrack_ambiguous_gene_mapping.txt
    misses=$(wc -l output/01_Description/isoTrack_ambiguous_gene_mapping.txt | awk '{print $1}')
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] ${misses} genes have ambiguous chromosome mapping"
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Removing associated DCRS..."
    grep -vwf output/01_Description/isoTrack_ambiguous_gene_mapping.txt output/01_Description/homo_sapiens_filter.${IDTY}pcnt.tsv > output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv
    # store retained DCR
    RETAINED=$(awk 'NR > 1 {print $1}' output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv | sort | uniq)
fi
finalDCRcnt=$(awk 'NR > 1 {print $1}' output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv | sort | uniq | wc -l | awk '{print $1}')
finalGenecnt=$(awk 'NR > 1 {print $2}' output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv | sort | uniq | wc -l | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] After removing ambiguously mapped genes, the filtered DCR table has ${finalDCRcnt}/${DCRIN} DCRs for ${finalGenecnt} genes"

# Addition 09-26-25: remove DCR with ambiguous isoforms mappings (wrong start or end due to small first or last exon)
# get the problematic isoforms mapping
bedtools getfasta -s -name -fi data/hg38.fa -bed <(cut -f 1-6 output/01_Description/Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguous.extBed.bed | awk '{if($6 == "+"){print $1"\t"$2"\t"$2+3"\t"$4"|start\t\.\t"$6; print $1"\t"$3"\t"$3+3"\t"$4"|end\t\.\t"$6} else {print $1"\t"$2-3"\t"$2"\t"$4"|end\t\.\t"$6; print $1"\t"$3-3"\t"$3"\t"$4"|start\t\.\t"$6}}') | awk 'getline seq {print $0"\t"toupper(seq)}' | awk '/start/{if($2 == "ATG"){print $0"\tPASS"} else {print $0"\tFAIL"}; next} /end/{if($2 == "TAA" || $2 == "TAG" || $2 == "TGA"){print $0"\tPASS"} else {print $0"\tFAIL"}}' | grep 'FAIL' | awk '{split($1, name, "::"); print name[2]"\t"name[1]"\t"$2"\t"$3}' | sed 's/(/\t/g;s/)//g' | awk '{split($1,chrpos,":"); print chrpos[1]"\t"chrpos[2]"\t"$3"\t"$2"\t"$4"\t"$5}' | awk '{split($2,pos,"-"); print $1"\t"pos[1]"\t"pos[2]"\t"$3"|"$5"|"$6"\t\.\t"$4}' > problematic_end_regions.bed
# find the DCRs to remove
join -11 -21 <(cut -f 4 problematic_end_regions.bed | sort | uniq | sed 's/|end/\t|end/g;s/|start/\t|start/g;s/>//g' | sort -k1,1) <(sort -k1,1 output/01_Description/isoforms_to_DCR_metadata_canon.tsv) | awk '{print $5}' | sort | uniq > DCR_to_remove_due_to_isoforms_problematic_end_regions.txt
# remove them
cat <(head -n1 output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv) <(grep -vf <(awk '{print "^"$0"\\b"}' DCR_to_remove_due_to_isoforms_problematic_end_regions.txt) <(awk 'NR > 1' output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv)) > output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes_fixed.tsv
mv output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes_fixed.tsv output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv

# classify with python script ###################

OUT3_1="output/01_Description/unique_series_counts.csv"
OUT3_2="output/01_Description/DCR_pattern_matrix.csv"
OUT3_3="output/01_Description/pattern_correspondence.csv"
LOG3="logs/01_Description/classify.DCR.log"
if [ "$FLAG" = true ] && [ -s "$OUT3_1" ] && [ -s "$OUT3_2" ] && [ -s "$OUT3_3" ]; then
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping classify.DCR.py (log exists and not empty)"
else
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Classifying DCR 1/2"
    # edit new enhanced classification method
    python3.11 scripts/01_Description/categorize_DCR_with_up_down_exons_edges_k16.py --input output/01_Description/homo_sapiens_filter.80pcnt.noAmbigGenes.tsv --out-detailed output/01_Description/NEW_DCR_architectures_exon_edges.tsv --out-upset output/01_Description/NEW_DCR_architectures_exon_edges_upset.tsv
    #python3.11 scripts/01_Description/categorize_DCR_with_up_down_exons.py -i output/01_Description/homo_sapiens_filter.${IDTY}pcnt.noAmbigGenes.tsv > "$LOG3" 2>&1
fi




# DCR exons count and distribution ##########

# NOTE/TO DO: for now I let the error print in stdout: the # of exon and distribution will need to be adjusted manually in the R script for aestetic purposes.
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Running DCR exon count and distribution R script"
Rscript scripts/01_Description/DCR_exon_length_distribution.r output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb.bed
if ! [ -s "output/01_Description/DCR_exon_length_distribution.pdf" ]; then
    echo "[warning $(date '+%Y-%m-%d %H:%M:%S')] The exon length distribution graph was not produced. Please check the R script for errors."
fi

## Create isoform track #####################
OUT8="output/01_Description/UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed"
if [ "$FLAG" = true ] && [ -s "$OUT8" ]; then
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping isoform tracks creation (already exists and not empty)"
else
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Creating isoform tracks"
    # then (this can probably be simplified to match required input)
    RAWISOCNT=$(cut -f 4 output/01_Description/UniprotKB_isoforms_mapping.bed | sort | uniq | wc -l | awk '{print $1}')
    RAWISOGENECNT=$(cut -f 4 output/01_Description/UniprotKB_isoforms_mapping.bed | sed 's/\|/\t/g' | cut -f 2 | sort | uniq | wc -l | awk '{print $1}')
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] ${RAWISOCNT} isoforms from ${RAWISOGENECNT} genes in input"
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Convert tracks for visualization and NMD evaluation"
    # convert the bed track for the NMD script
    awk '{split($4, iso, "[\\|/]"); print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"iso[2]}' output/01_Description/UniprotKB_isoforms_mapping.bed > output/01_Description/UniprotKB_isoforms_mapping_forNMD.bed
    # remove ambiguous genes
    grep -vwf output/01_Description/isoTrack_ambiguous_gene_mapping.txt output/01_Description/UniprotKB_isoforms_mapping_forNMD.bed > output/01_Description/UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed
    FISOCNT=$(cut -f 4 output/01_Description/UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed | sort | uniq | wc -l | awk '{print $1}')
    FISOGENECNT=$(cut -f 4 output/01_Description/UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed | sed 's/\|/\t/g' | cut -f 2 | sort | uniq | wc -l | awk '{print $1}')
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] ${FISOCNT} isoforms from ${FISOGENECNT} genes are ready for NMD evaluation"
    # convert the bed track to extended bed for visualization, only for isoforms in retained DCRs (it's long)
    OUT9="output/01_Description/Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguous.extBed.bed"
    if [ "$FLAG" = true ] && [ -s "$OUT9" ]; then
        echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping isoform extended bed track creation (already exists and not empty)"
    else
        python3.11 scripts/01_Description/convert_bed_to_ExtendedBed.py -i output/01_Description/UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed |\
         awk '{$9 = "255,0,0"; print $0}' | sed 's/ /\t/g' > output/01_Description/Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguous.extBed.bed
    fi
fi  

# Addition 09-26-25: remove DCR with ambiguous isoforms mappings (wrong start or end due to small first or last exon)
grep -vf <(awk '{print "^"$0"\\b"}' DCR_to_remove_due_to_isoforms_problematic_end_regions.txt) <(awk '{print $4"\t"$0}' output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb.bed) > output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb_fixed.bed

# NMD susceptibility ########################
# extract filtered DCRs IDs
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Testing for NMD susceptibility..."
cut -f 4 output/01_Description/filtered_DCE_overlapping_regions_with_gene_and_dceID_and_frameNb.bed | sort | uniq > output/01_Description/filtered.DCR.id.txt
NMDIN=$(wc -l output/01_Description/filtered.DCR.id.txt | awk '{print $1}')
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] ${NMDIN} DCRs will be tested for NMD susceptibility"
# run the script
if [ "$FLAG" = true ] && [ -s "output/01_Description/DCR_NMD_analysis.out" ]; then
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Skipping NMD analysis (already exists and not empty)"
else
    # create isoform metadata file 
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Running NMD analysis..."
    ./scripts/01_Description/DCR_isoform_group_to_frame.sh ${INDIR}/homo_sapiens > output/01_Description/isoforms_to_DCR_metadata.tsv
    # run NMD script
    python3.11 scripts/01_Description/isoform_NMD.py
fi  

# add uncertainty to per-isoform DCR assessment
awk '{if ($5 == "no_NMD" && $4 <= 50) print $0 "\tnot_enough_evidence"; else if ($5 == "no_NMD" && $4 == "last_exon") print $0 "\tno_NMD"; else print $0 "\tNMD"}' output/01_Description/NMDtable.tsv > output/01_Description/NMDtable_updated.tsv

# Addition 09-26-25: remove DCR with ambiguous isoforms mappings (wrong start or end due to small first or last exon)
grep -vf <(awk '{print "^"$0"\\b"}' DCR_to_remove_due_to_isoforms_problematic_end_regions.txt) output/01_Description/DCR_NMD_analysis.out > output/01_Description/DCR_NMD_analysis.out_fixed
mv output/01_Description/DCR_NMD_analysis.out_fixed output/01_Description/DCR_NMD_analysis.out

# Splice site analysis ######################

# DCR
./scripts/01_Description/extract_DCE_splice_sites.sh ${INDIR}/homo_sapiens/ > output/01_Description/EXPERIMENTAL_filtered_DCE_splice_sites.bed

# Addition 09-26-25: remove DCR with ambiguous isoforms mappings (wrong start or end due to small first or last exon)
grep -vf <(awk '{print "DCEid="$0";"}' DCR_to_remove_due_to_isoforms_problematic_end_regions.txt) output/01_Description/EXPERIMENTAL_filtered_DCE_splice_sites.bed > output/01_Description/EXPERIMENTAL_filtered_DCE_splice_sites.bed_fixed
mv output/01_Description/EXPERIMENTAL_filtered_DCE_splice_sites.bed_fixed output/01_Description/EXPERIMENTAL_filtered_DCE_splice_sites.bed

# back to normal
bedtools getfasta -s -bed <(grep 'junction=upstream;site=donor' output/01_Description/EXPERIMENTAL_filtered_DCE_splice_sites.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t\.\t"$5}'|\
 awk '{print $1"_"$2"_"$3"\t"$0}' | sort -u -k1,1 | cut -f 2-) -fi data/hg38.fa -name > output/01_Description/EXPERIMENTAL_UPSTREAM_DONOR_filtered_DCE_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=upstream;site=acceptor' output/01_Description/EXPERIMENTAL_filtered_DCE_splice_sites.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t\.\t"$5}'|\
 awk '{print $1"_"$2"_"$3"\t"$0}' | sort -u -k1,1 | cut -f 2-) -fi data/hg38.fa -name > output/01_Description/EXPERIMENTAL_UPSTREAM_ACCEPTOR_filtered_DCE_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=downstream;site=donor' output/01_Description/EXPERIMENTAL_filtered_DCE_splice_sites.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t\.\t"$5}'|\
 awk '{print $1"_"$2"_"$3"\t"$0}' | sort -u -k1,1 | cut -f 2-) -fi data/hg38.fa -name > output/01_Description/EXPERIMENTAL_DOWNSTREAM_DONOR_filtered_DCE_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=downstream;site=acceptor' output/01_Description/EXPERIMENTAL_filtered_DCE_splice_sites.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t\.\t"$5}'|\
 awk '{print $1"_"$2"_"$3"\t"$0}' | sort -u -k1,1 | cut -f 2-) -fi data/hg38.fa -name > output/01_Description/EXPERIMENTAL_DOWNSTREAM_ACCEPTOR_filtered_DCE_splice_sites.fa
# print logos
/opt/homebrew/bin/weblogo -f output/01_Description/EXPERIMENTAL_UPSTREAM_DONOR_filtered_DCE_splice_sites.fa -o output/01_Description/LOGO_DCR_upstream_donor.pdf -A DNA -U probability -F PDF -t "upstream donor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/EXPERIMENTAL_UPSTREAM_ACCEPTOR_filtered_DCE_splice_sites.fa -o output/01_Description/LOGO_DCR_upstream_acceptor.pdf -A DNA -U probability -F PDF -t "upstream acceptor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/EXPERIMENTAL_DOWNSTREAM_DONOR_filtered_DCE_splice_sites.fa -o output/01_Description/LOGO_DCR_downstream_donor.pdf -A DNA -U probability -F PDF -t "downstream donor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/EXPERIMENTAL_DOWNSTREAM_ACCEPTOR_filtered_DCE_splice_sites.fa -o output/01_Description/LOGO_DCR_downstream_acceptor.pdf -A DNA -U probability -F PDF -t "downstream acceptor site" -s large -c classic

# controls
python3.11 scripts/01_Description/extract_control_exons_splice_site.py
bedtools getfasta -s -bed <(grep 'junction=upstream;site=donor' output/01_Description/EXPERIMENTAL_control_splice_sites.bed |\
 sort | uniq) -fi data/hg38.fa -name > output/01_Description/EXPERIMENTAL_control_UPSTREAM_DONOR_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=upstream;site=acceptor' output/01_Description/EXPERIMENTAL_control_splice_sites.bed |\
 sort | uniq) -fi data/hg38.fa -name > output/01_Description/EXPERIMENTAL_control_UPSTREAM_ACCEPTOR_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=downstream;site=donor' output/01_Description/EXPERIMENTAL_control_splice_sites.bed |\
 sort | uniq) -fi data/hg38.fa -name > output/01_Description/EXPERIMENTAL_control_DOWNSTREAM_DONOR_splice_sites.fa
bedtools getfasta -s -bed <(grep 'junction=downstream;site=acceptor' output/01_Description/EXPERIMENTAL_control_splice_sites.bed |\
    sort | uniq) -fi data/hg38.fa -name > output/01_Description/EXPERIMENTAL_control_DOWNSTREAM_ACCEPTOR_splice_sites.fa
/opt/homebrew/bin/weblogo -f output/01_Description/EXPERIMENTAL_control_UPSTREAM_DONOR_splice_sites.fa -o output/01_Description/LOGO_control_upstream_donor.pdf -A DNA -U probability -F PDF -t "upstream donor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/EXPERIMENTAL_control_UPSTREAM_ACCEPTOR_splice_sites.fa -o output/01_Description/LOGO_control_upstream_acceptor.pdf -A DNA -U probability -F PDF -t "upstream acceptor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/EXPERIMENTAL_control_DOWNSTREAM_DONOR_splice_sites.fa -o output/01_Description/LOGO_control_downstream_donor.pdf -A DNA -U probability -F PDF -t "downstream donor site" -s large -c classic
/opt/homebrew/bin/weblogo -f output/01_Description/EXPERIMENTAL_control_DOWNSTREAM_ACCEPTOR_splice_sites.fa -o output/01_Description/LOGO_control_downstream_acceptor.pdf -A DNA -U probability -F PDF -t "downstream acceptor site" -s large -c classic
echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Pipeline step 01 (DCR Description) completed"

###### the code below was required to generate the list of DCR to exclude due to incorrect mapping which generate the file "DCR_to_remove_due_to_isoforms_problematic_end_regions.txt" #########

# bedtools getfasta -s -name -fi data/hg38.fa -bed <(cut -f 1-6 output/01_Description/Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguous.extBed.bed | awk '{if($6 == "+"){print $1"\t"$2"\t"$2+3"\t"$4"|start\t\.\t"$6; print $1"\t"$3"\t"$3+3"\t"$4"|end\t\.\t"$6} else {print $1"\t"$2-3"\t"$2"\t"$4"|end\t\.\t"$6; print $1"\t"$3-3"\t"$3"\t"$4"|start\t\.\t"$6}}') | awk 'getline seq {print $0"\t"toupper(seq)}' | awk '/start/{if($2 == "ATG"){print $0"\tPASS"} else {print $0"\tFAIL"}; next} /end/{if($2 == "TAA" || $2 == "TAG" || $2 == "TGA"){print $0"\tPASS"} else {print $0"\tFAIL"}}' | grep 'FAIL' | awk '{split($1, name, "::"); print name[2]"\t"name[1]"\t"$2"\t"$3}' | sed 's/(/\t/g;s/)//g' | awk '{split($1,chrpos,":"); print chrpos[1]"\t"chrpos[2]"\t"$3"\t"$2"\t"$4"\t"$5}' | awk '{split($2,pos,"-"); print $1"\t"pos[1]"\t"pos[2]"\t"$3"|"$5"|"$6"\t\.\t"$4}' > problematic_end_regions.bed
# join -11 -21 <(cut -f 4 problematic_end_regions.bed | sort | uniq | sed 's/|end/\t|end/g;s/|start/\t|start/g;s/>//g' | sort -k1,1) <(sort -k1,1 output/01_Description/isoforms_to_DCR_metadata_canon.tsv) |sed 's/ /\t/g' > problematic_end_regions_and_DCR_metadata.tsv
# join -11 -21 <(cut -f 4 problematic_end_regions.bed | sort | uniq | sed 's/|end/\t|end/g;s/|start/\t|start/g;s/>//g' | sort -k1,1) <(sort -k1,1 output/01_Description/isoforms_to_DCR_metadata_canon.tsv) | awk '{print $5}' | sort | uniq > DCR_to_remove_due_to_isoforms_problematic_end_regions.txt