#!/usr/bin/env bash
filelist=${1}/*

for file in ${filelist[@]}
do
	#echo $file
	CHR=$(grep 'Chromosome' ${file} | awk '{split($3,chr,"["); if(chr[2]~"revcomp") {print chr[1]"\t-"} else {print chr[1]"\t+"}}')
	DCE=$(grep 'DCE Index' ${file} | awk '{print $4}') #echo "$CHR"
	if grep -qw ${DCE} output/01_Description/filtered.DCR.id.txt; then
		STRND=$(echo ${CHR} | awk '{print $2}')
		grep '> Group\|Genome Ranges:' $file | awk 'getline range {split(range,coord,"|"); print $2"_"$3"\t"$NF" "coord[1]"\t"coord[2]"\t"coord[3]}'| sed -E 's/ \/ /\//g;s/://g;s/\.\.[0-9]+\/[0-9]+\.\./\.\./g' |\
		 awk -v chr="${CHR}" -v dce="${DCE}" '{split($5,left,"\\.\\."); split($6,main,"\\.\\."); split($7,right,"\\.\\."); isoform = $2; if(left[1]~/N-TERM/){print chr"\t"left[1]"\t"main[1]"\t"isoform";DCEid="dce";junction=upstream"; print chr"\t"main[length(main)]"\t"right[1]"\t"isoform"_DCEid:"dce";junction=downstream"} else {print chr"\t"left[2]"\t"main[1]"\t"isoform";DCEid="dce";junction=upstream"; print chr"\t"main[length(main)]"\t"right[1]"\t"isoform"_DCEid:"dce";junction=downstream"}}' |\
		 grep -v 'C-TERM\|N-TERM' |\
		 awk '{
				    if ($2 == "+") {
				    	if($5~/upstream/){
					        print $1"\t"($3)"\t"($3+10)"\t"$5";site=donor\t+";
					        print $1"\t"($4-11)"\t"($4-1)"\t"$5";site=acceptor\t+";
				        } else {
					        print $1"\t"($3)"\t"($3+10)"\t"$5";site=donor\t+";
					        print $1"\t"($4-11)"\t"($4-1)"\t"$5";site=acceptor\t+";        	
				        }
				    } else {
				    	if($5~/upstream/){
					        print $1"\t"($4)"\t"($4+10)"\t"$5";site=acceptor\t-";
					        print $1"\t"($3-11)"\t"($3-1)"\t"$5";site=donor\t-";
					        } else {
					        print $1"\t"($4)"\t"($4+10)"\t"$5";site=acceptor\t-";
					        print $1"\t"($3-11)"\t"($3-1)"\t"$5";site=donor\t-";        	
				        }
				    }
				}'
	fi
done
