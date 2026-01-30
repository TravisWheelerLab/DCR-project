#!/bin/bash

mkdir -p models/dssp_out/parsed
# get the list of isoforms to process 
isoforms=$(awk '{print $1}' DCRpep.DCRiso.start.end.tsv)

for iso in $isoforms; do
    # get the regex for the corresponding dssp cif file
    isoregex=$(grep -w "${iso}" DCRpep.DCRiso.start.end.tsv | cut -f 2 | sed 's/|/\t/g;s/DCR=/dcr_/g;s/frame=/frame_/g' | awk '{print $3"_"$4"_"$5}')
    form=$(grep -w "${iso}" DCRpep.DCRiso.start.end.tsv | sed 's/|/\t/g;s/DCR=/dcr_/g;s/frame=/frame_/g' | cut -f 6)
    DCR=$(echo "${iso}" | sed 's/_/\t/g;s/iso//g;s/f//g' | cut -f 1)
    isoNB=$(echo "${iso}" | sed 's/_/\t/g;s/iso//g;s/f//g' | cut -f 2)
    frameNB=$(echo "${iso}" | sed 's/_/\t/g;s/iso//g;s/f//g' | cut -f 3)
    dssp_file=$(ls models/dssp_out | grep "${isoregex}")
    # get start and end withing the isoform
    start=$(grep -w "${iso}" DCRpep.DCRiso.start.end.tsv | cut -f 3)
    end=$(grep -w "${iso}" DCRpep.DCRiso.start.end.tsv | cut -f 4)
    # identify dssp cif file region with the structure summary
    start_line=$(grep -n '_dssp_struct_summary.z_ca' models/dssp_out/${dssp_file} | cut -d: -f1 | awk '{print $1+1}')
    start_line=$((start_line + start - 1))
    end_line=$((start_line + end - start))
       
    echo "${iso}: parsing lines ${start_line} to ${end_line} from file models/dssp_out/${dssp_file} for residues ${start} to ${end}"
    # extract in a parsable tsv
    awk -v start=${start_line} -v end=${end_line} -v dcr=${DCR} -v iso=${isoNB} -v frame=${frameNB} -v form=${form} 'NR >= start && NR <= end {print dcr"\t"iso"\t"frame"\t"$3"\t"$4"\t"$5"\t"form}' models/dssp_out/${dssp_file} > models/dssp_out/parsed/${iso}_dssp_parsed.tsv
done

cat  models/dssp_out/parsed/* > DCR_peps_structure_annotation.tsv