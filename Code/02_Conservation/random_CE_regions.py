#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
import random
from alive_progress import alive_bar

# argument parser
parser = argparse.ArgumentParser(description='Create random exon intervals according to DCE regions length distribution')
parser.add_argument('-n', '--number', type = str, metavar='<int>', help='number of intervals to generate', dest = 'n_intervals', default = '')
parser.add_argument('-e', '--exons', type = str, metavar='<.bed>', help='mouse exon coordinates', dest = 'exons', default = '')
parser.add_argument('-d', '--dce', type = str, metavar='<.bed>', help='human DCE regions coordinates', dest = 'dces', default = '')
# parse the arguments
args = parser.parse_args()

# check we have the inputs
if args.n_intervals == '': # number of intervals to generate
    print("ERROR: missing input number of intervals (-n <int>)")
    parser.print_help()
    sys.exit()
print(args)
if args.exons == '': # list of human exons in .bed
    print("ERROR: missing input exons list (-e <exons.bed>)")
    parser.print_help()
    sys.exit()
print(args)
if args.dces == '': # list of human DCE regions in .bed
    print("ERROR: missing input dce regions (-e <dces.bed>)")
    parser.print_help()
    sys.exit()
print(args)

def generate_intervals(n_ints = args.n_intervals, exons = args.exons, dces = args.dces):
	# read the exon table
	read_exons = pd.read_csv(exons, sep='\t', header=None)
	read_exons = read_exons.rename(columns={0: "chr", 1: "start", 2: "end", 3: "name", 4: "score", 5: "strand"})
	# read the DCE table
	read_dces = pd.read_csv(dces, sep='\t', header=None)
	read_dces = read_dces.rename(columns={0: "chr", 1: "start", 2: "end", 3: "name", 4: "max_overlap", 5: "strand"})
	# calculate the length of regions and exons
	read_dces['length'] = read_dces['end'] - read_dces['start']
	read_exons['length'] = read_exons['end'] - read_exons['start']
	
	# create an empty list to host the chosen intervals 
	chosen_ints = []
	# here start the main loop:
	with alive_bar(int(args.n_intervals), bar = 'circles', spinner = 'classic') as bar:
		for i in range(int(args.n_intervals)):
			# pick a region length
			reg_serie = read_dces['length'].sample()
			reg_len = int(reg_serie.iloc[0])
			#print(reg_len)
			# pick an exon
			exon_picked = read_exons.sample()
			# test if it is long enough
			exon_picked_len = exon_picked['length'].iloc[0]
			while exon_picked_len < reg_len:
				# we try another one
				exon_picked = read_exons.sample()
				exon_picked_len = exon_picked['length'].iloc[0]
			# now define the range where we can safely pick a start
			safe_start = exon_picked['end'].iloc[0] - reg_len
			new_start = random.randint(exon_picked['start'].iloc[0], safe_start)
			new_end = new_start + reg_len
			new_record = [exon_picked['chr'].iloc[0], new_start, new_end, exon_picked['name'].iloc[0], exon_picked['score'].iloc[0], exon_picked['strand'].iloc[0] ]
			#print(new_record)
			chosen_ints.append(new_record)
			bar()
		# turn into df
		output_df = pd.DataFrame(chosen_ints, columns = ['chr', 'start', 'end', 'name', 'score', 'strand']) 
		output_df.to_csv("output/02_Conservation/fake_DCE_intervals.bed", sep = "\t", index=False)

# execute code
generate_intervals()