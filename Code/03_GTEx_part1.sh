#!/usr/bin/env bash

###################################################
# 03 GTEx analyses - PART 1: junction preparation #
###################################################

# Parse arguments
usage() {
    echo "Usage: $0 -i <Extracted_DCR_directory> [-s]"
    echo "  -s    Skip command if output file exists and is not empty (default = false)"
    echo "  -i    Input dir (required)"
    exit 1
}

FLAG=false
INPUT_DIR="ata/Extracted-DCEs_06202025"

while getopts ":si:I:" opt; do
  case $opt in
    s)
      FLAG=true
      ;;
    i)
      INPUT_DIR="$OPTARG"
      ;;
    \?)
      usage
      ;;
  esac
done

# Check if input directory is provided
if [ -z "$INPUT_DIR" ]; then
    echo "Error: Input directory is required."
    usage
fi  

OUT1="output/03_GTEx/unique_junctions.fa"
if [ "$FLAG" = true ] && [ -s "$OUT1" ]; then
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Unique junctions file $OUT1 already exists and is not empty. Skipping command."
else
    echo "[info $(date '+%Y-%m-%d %H:%M:%S')] Creating junction file to quantify"
    # Linearize needles in order to find duplicates
    awk 'getline seq {print $0"\t"toupper(seq)}' data/Extracted-DCEs_06202025/homo_sapiens.IREs-removed.fa > output/03_GTEx/homo_sapiens.IREs-removed.fa.tsv
    # count occurence of each needle
    cut -f 2 output/03_GTEx/homo_sapiens.IREs-removed.fa.tsv | sort | uniq -c | awk '{print "junction_"NR"_"$1"\t"$2"\t"$1}' > output/03_GTEx/junctions_uniq_name_and_count.tsv
    # create a correspondance table for demultiplexing
    join -12 -22 <(sort -k2,2 output/03_GTEx/homo_sapiens.IREs-removed.fa.tsv) <(sort -k2,2 output/03_GTEx/junctions_uniq_name_and_count.tsv) | awk '{print $2"\t"$1"\t"$3}' > output/03_GTEx/junctions_name_correspondance_for_demultiplexing.tsv
    # create a fasta of unique needles
    awk '{print ">"$1; print $2}' output/03_GTEx/junctions_uniq_name_and_count.tsv > output/03_GTEx/unique_junctions.fa
fi

### end of Part 1 ###