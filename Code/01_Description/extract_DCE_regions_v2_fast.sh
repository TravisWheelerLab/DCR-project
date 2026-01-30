#!/usr/bin/env bash
# Description: This script parse individual DCE .out file to extract the dual coding region of each DCE.
# Usage: ls Extracted-DCEs_070124/homo_sapiens/* | parallel ./extract_DCE_regions_v2_fast.sh

filelist=$1 #(ls Extracted-DCEs_070124/homo_sapiens/*)
#filelist="Extracted-DCEs_070124/homo_sapiens/98.myb.out"

for file in "${filelist}"
do
	CHR=$(grep 'Chromosome' ${file} | awk '{split($3,chr,"["); if(chr[2]~"revcomp") {print chr[1]"\t-"} else {print chr[1]"\t+"}}')
	DCE=$(grep 'DCE Index' ${file} | awk '{print $4}') #echo "$CHR"
	STRND=$(echo ${CHR} | awk '{print $2}')
	#grouping=$(grep '> Group\|Genome Ranges:' $file |  awk 'getline range {split(range,coord,"|"); print $2"_"$3"\t"$NF" "coord[2]}' | sed 's/ \/ /\//g;s/://g' | awk -v chr="${CHR}" '{split($3,coords,"/"); for(i = 1; i <= length(coords); i++){print $1"\t"$2"\t"chr"\t"coords[i]"\n"}}') #| sed 's/Group_//g')
	grouping=$(grep '> Group\|Genome Ranges:' $file |  awk 'getline range {split(range,coord,"|"); print $2"_"$3"\t"$NF" "coord[2]}' | sed 's/ \/ /\//g;s/://g' | awk -v chr="${CHR}" '{split($3,coords,"/"); for(i = 1; i <= length(coords); i++){print $1"\t"$2".."coords[i]"\t"chr"\n"}}') #| sed 's/Group_//g')
	echo -e "$grouping"
	#frames=$(grep '   Group' $file | grep 'frame' | awk '{print $1"_"$2"\t"$NF"\n"}' | sed 's/\[//g;s/\]//g;s/frame://g') # > frames.tmp
	#frames=$(grep '   Group' $file | grep 'frame' | awk '{print $1"_"$2"\t"$NF"\n"}' | sed 's/\[//g;s/\]//g;s/frame://g' | awk -F '\t' '{split($1, groups, ","); count = $2; for (i in groups) {gsub(/^Group[s]?_/, "", groups[i]); print "Group_"groups[i] "\t" count}}')
	#echo -e "$frames"
	# now join group to frame and correct coordinates by doing -1 to start position (bed format from DCE .out format)
	#GroupedFrames=$(join -11 -21 <(echo "$grouping" | grep 'Group' | sort -k1,1 ) <(echo "$frames" | grep 'Group' | sort -k1,1) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6}' | sed 's/\.\./\t/g' | awk '{if($3 == "+") {print $2"\t"$4-1"\t"$5"\t"$1"\t\.\t"$3"\t"$NF} else {print $2"\t"$5-1"\t"$4"\t"$1"\t\.\t"$3"\t"$NF}}' | sort -k1,1 -k2,2n)
	GroupedFrames=$(join -12 -21 <(echo "$grouping" | grep 'Group' | sort -k2,2 ) <(sort -k1,1 output/01_Description/DCE_exons_with_phase.tsv) | sed 's/\.\./\t/g' | awk '{if($6 == "+") {print $5"\t"$2-1"\t"$3"\t"$1"\t\.\t"$6"\t"$7} else {print $5"\t"$3-1"\t"$2"\t"$1"\t\.\t"$6"\t"$NF}}' | sort -k1,1 -k2,2n)
	echo "$GroupedFrames"
	#echo ""
	# we need to keep the longest per frame	
	w_length=""
	FrameNumbers=$(echo "$GroupedFrames" | awk '{print $NF}' | sort | uniq)
	#echo "$GroupedFrames" | awk '{print $NF}' | sort | uniq
	#echo "$FrameNumbers"
	for i in ${FrameNumbers[@]} #`seq 1 3` # 3 frames at best
		do
		# echo "${i}"
		# echo ""
		#fm=$(echo "${i}" | tr -d '\n')
		#echo "$fm"
		#echo "$GroupedFrames" | awk -v var="${fm}" '$NF==var' # <(echo "$GroupedFrames")
		#echo "$GroupedFrames" | awk -v var=${i} '$NF == var' | sort -k1,1 -k2,2n | bedtools merge -s -c 4,5,6,7 -o distinct -i -
		#echo "---"
		w_length+="$(echo "$GroupedFrames" | awk -v var=${i} '$NF == var' | sort -k1,1 -k2,2n | bedtools merge -s -c 4,5,6,7 -o distinct -i -)\n"
		done
	# and then bedtools genomecov to get the DCE region
	w_len_sorted=$(echo -e "$w_length" | sed '/^$/d' | awk '{print $1"::"$2"::"$3"\t"$0}' | sort -u -k1,1 | cut -f 2-)
	#echo "---"
	#echo -e "$w_len_sorted" # | sort -k1,1 -k2,2n | bedtools genomecov -i - -g hg38.genome -bg | sort -k4,4nr | awk -v DCE="${DCE}" -v STRND="${STRND}" '$4 >= 2 {print $1"\t"$2"\t"$3"\t"DCE"\t"$4"\t"STRND}' | sort -k1,1 -k2,2n | bedtools merge -s -i - | awk -v DCE="${DCE}" -v STRND="${STRND}" '{print $0"\t"DCE"\t\.\t"STRND}' #>> DCE_overlapping_regions.bed
	#echo -e "$w_length"
	echo -e "$w_length" | sort -k1,1 -k2,2n | bedtools genomecov -i - -g data/hg38.genome -bg | sort -k4,4nr | awk -v DCE="${DCE}" -v STRND="${STRND}" '$4 >= 2 {print $1"\t"$2"\t"$3"\t"DCE"\t"$4"\t"STRND}' | sort -k1,1 -k2,2n | bedtools merge -s -c 4,5,6 -o distinct -i - | awk '{print $1"\t"$2"\t"$3"\tDCEid="$4";nFrames="$5"\t\.\t"$6}' # >> DCE_overlapping_regions.bed
	#| awk -v DCE=${DCE} -v STRND=${STRND} 'BEGIN {max=0} {if($4>max) {print $1"\t"$2"\t"$3"\t"DCE"\t"$4"\t"STRND; max=$4; next}}' #>> DCE_overlapping_regions.bed
done
