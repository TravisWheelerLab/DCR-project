#!/usr/bin/env bash
DCRDIR=$1 #data/Extracted-DCEs_06202025/homo_sapiens/*
filelist=$(ls ${DCRDIR}/*)
#filelist=($(ls Extracted-DCEs_070124/homo_sapiens/ | grep -wf DCE_NMD_Candidates | awk '{print "Extracted-DCEs_070124/homo_sapiens/"$0}'))

#echo $filelist
for file in ${filelist}
do
	#echo "$file"
	CHR=$(grep 'Chromosome' ${file} | awk '{split($3,chr,"["); if(chr[2]~"revcomp") {print chr[1]"\t-"} else {print chr[1]"\t+"}}')
	DCE=$(grep 'DCE Index' ${file} | awk '{print $4}') #echo "$CHR"
	STRND=$(echo ${CHR} | awk '{print $2}')
	#grouping=$(grep '> Group\|Genome Ranges:' $file |  awk 'getline range {split(range,coord,"|"); print $2"_"$3"\t"$NF" "coord[2]"\t"coord[3]}' | sed 's/ \/ /\//g;s/://g' | awk -v chr="${CHR}" '{split($3,coords,"/"); for(i = 1; i <= length(coords); i++){print $1"\t"$2"\t"chr"\t"coords[i]"\t"$NF"\n"}}') #| sed 's/Group_//g')
	grouping=$(cat $file | grep 'homo_sapiens|' | awk '{print $3"\t"$NF}' |  awk -F'\t' 'BEGIN { OFS=FS } /^[0-9]+:/ { key=$1 } /^[^0-9]/ { $0 = key FS $2 } 1' | awk '{sub(/:/,"",$1); print "Group_"$0}')
	#echo -e "$grouping"
	#frames=$(grep '   Group' $file | grep 'frame' | awk '{print $1"_"$2"\t"$NF"\n"}' | sed 's/\[//g;s/\]//g;s/frame://g') # > frames.tmp
	frames=$(grep '   Group' $file | grep 'frame' | awk '{print $1"_"$2"\t"$NF"\n"}' | sed 's/\[//g;s/\]//g;s/frame://g' | awk -F '\t' '{split($1, groups, ","); count = $2; for (i in groups) {gsub(/^Group[s]?_/, "", groups[i]); print "Group_"groups[i] "\t" count}}')
	#echo -e "$frames"
	# now join group to frame 
	GroupedFrames=$(join -11 -21 <(echo "$grouping" | sort -k1,1 ) <(echo "$frames" | sort -k1,1)) #| awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sed 's/\.\./\t/g' | awk '{if($3 == "+") {print $2"\t"$4-1"\t"$5"\t"$1"\t\.\t"$3"\t"$(NF-1)"\t"$NF} else {print $2"\t"$5-1"\t"$4"\t"$1"\t\.\t"$3"\t"$(NF-1)"\t"$NF}}' | sort -k1,1 -k2,2n)
	# get the termini of each group
	Gend=$(cat $file | grep '> Group\|Genome Ranges' | awk 'getline range {print $3"\t"range}' | awk '{sub(/:/, "", $1); print "Group_"$1"\t"$NF}')
	# joint it to group
	# echo "$GroupedFrames"
	# echo "$Gend"
	Final=$(join -11 -21 <(echo "$GroupedFrames" | sort -k1,1 ) <(echo "$Gend" | sort -k1,1) | sed 's/ /\t/g')
	# filter to pick one representative isoform per frame
	#fmtd_dce=$(echo "$GroupedFrames" | sort -k3,3n | sort -u -k8,8)
	if echo ${DCE} | grep -q -wf output/01_Description/filtered.DCR.id.txt; then
		echo -e "$Final" | awk -v dce=${DCE} '{print $2"\t"$3"\t"$4"\t"dce}' #| sed 's/frames://g'
	fi
done