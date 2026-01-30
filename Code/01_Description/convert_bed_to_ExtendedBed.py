#!/usr/bin/env python3
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Convert BED to Extended BED format.")
parser.add_argument('-i', '--input', required=True, help='Input BED file')
args = parser.parse_args()
input = args.input

def main(input=input):
	bed_in = pd.read_csv(input, sep = "\t")
	isoforms = bed_in.iloc[:, 3].unique()
	for form in isoforms:
		current_isoform_lines = bed_in[bed_in.iloc[:, 3] == form]
		isoform_starts = current_isoform_lines.iloc[:, 1].tolist()
		isoform_starts.sort()
		isoform_ends = current_isoform_lines.iloc[:, 2].tolist()
		isoform_ends.sort()
		isoform_start = min(isoform_starts)
		isoform_block_starts = [pos - isoform_start for pos in isoform_starts]
		isoform_end = max(isoform_ends)
		isoform_chr = current_isoform_lines.iloc[:, 0].unique()
		isoform_strand = current_isoform_lines.iloc[:, 4].unique()
		isoform_block_count = len(isoform_starts)
		isoform_block_length = [isoform_ends_i - isoform_starts_i for isoform_starts_i, isoform_ends_i in zip(isoform_starts, isoform_ends)] # [a_i - b_i for a_i, b_i in zip(a, b)]
		print(str(isoform_chr[0]) + "\t",			# chr
			str(isoform_start) + "\t",				# start
			str(isoform_end) + "\t",				# end
			str(form) + "\t",						# name
			"1000\t",									# score
			str(isoform_strand[0]) + "\t",				# strand
			str(isoform_start) + "\t",				# thickStart
			str(isoform_end) + "\t",				# thickEnd
			"0\t",									# itemRGB
			str(isoform_block_count),				# block count
			','.join(map(str, isoform_block_length)),		# block sizes
			','.join(map(str, isoform_block_starts)))		# block starts

		# with open('isoforms_mapping_coordinates_withDCEid.bed') as input_bed:
		# 	lines = input_bed.readlines()  # Read all lines into a list
		# 	isoforms = list(set([line.split()[3] for line in lines]))
		# 	for form in isoforms:
		# 		iso_lines = []
		# 		for line in lines:
		# 			current_isoform = line.split()[3]
		# 			if current_isoform == form:
		# 				iso_lines += [line.strip()]
		# 		print("---")
		# 		print(iso_lines)

main()