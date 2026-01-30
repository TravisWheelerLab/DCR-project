#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#Created on Tue Jun 18 10:00:00 2024
#extract plddts at DCR regions from AF3 structures of DCR isoforms

DCRpepCoords="../../DCRpep.DCRiso.start.end.tsv" # output from 03_extract_DCR_peptides.py

rm -f extract_DCR_plddts.tsv
touch extract_DCR_plddts.tsv

filelist="*.CA.tsv"
for file in ${filelist}
do
# loop starts here
DCR=$(echo $file | sed 's/_/\t/g' | awk '{print $6 }') # DCR ID
DCRregex=$(echo ${DCR} | awk -v dcrid=${DCR} '{print "|DCR="dcrid"|"}') # escape dot for regex
isoform=$(echo $file | sed 's/_/\t/g' | awk '{print $3 }') # isoform ID (number)
gene=$(echo $file | sed 's/_/\t/g' | awk '{print $2 }') # gene name
version=$(echo $file | sed 's/_/\t/g' | awk '{if(NF == 10){print $7} else {print $9}}') # canonical or alt

DCRisoCoords=$(grep "${DCRregex}" ${DCRpepCoords} | grep "${version}") # get DCR isoform coords
start=$(echo ${DCRisoCoords} | awk '{print $3 }') # DCR start
end=$(echo ${DCRisoCoords} | awk '{print $4 }') # DCR end   

# write the output
awk -v dcr=${DCR} -v gene=${gene} -v iso=${isoform} -v type=${version} -v start=${start} -v end=${end} '$3 >= start && $3 <= end {print dcr"\t"gene"\t"iso"\t"$3"\t"$4"\t"$5"\t"type}' ${file} >> extract_DCR_plddts.tsv
done