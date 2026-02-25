#!/usr/bin/env python3
import argparse
import sys
import os
import numpy as np
import pandas as pd
import random
from alive_progress import alive_bar

# argument parser
parser = argparse.ArgumentParser(description='Create random exon intervals according to DCE regions length distribution')
parser.add_argument('-n', '--number', type = str, metavar='<int>', help='number of intervals to generate', dest = 'n_intervals', default = '')
parser.add_argument('-e', '--exons', type = str, metavar='<.bed>', help='mouse exon coordinates', dest = 'exons', default = '')
parser.add_argument('-d', '--dce', type = str, metavar='<.bed>', help='human DCE regions coordinates', dest = 'dces', default = '')
parser.add_argument('-f', '--fasta', type = str, metavar='<.fa>', help='mm39 fasta file (used to filter valid chromosomes via .fai index)', dest = 'fasta', default = '')
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
if args.fasta == '': # mm39 fasta file
    print("ERROR: missing input fasta file (-f <mm39.fa>)")
    parser.print_help()
    sys.exit()
print(args)

def get_valid_chromosomes(fasta_path):
	"""Read valid chromosome names from the .fai index of the fasta file."""
	fai_path = fasta_path + ".fai"
	if os.path.exists(fai_path):
		valid_chroms = set()
		with open(fai_path) as f:
			for line in f:
				valid_chroms.add(line.split('\t')[0])
		return valid_chroms
	else:
		# fallback: parse fasta headers directly
		print(f"[warning] .fai index not found at {fai_path}, parsing fasta headers instead (slow)...")
		valid_chroms = set()
		with open(fasta_path) as f:
			for line in f:
				if line.startswith('>'):
					valid_chroms.add(line[1:].split()[0])
		return valid_chroms

def generate_intervals(n_ints = args.n_intervals, exons = args.exons, dces = args.dces, fasta = args.fasta):
	# get valid chromosomes from the fasta index
	valid_chroms = get_valid_chromosomes(fasta)
	print(f"[info] Found {len(valid_chroms)} valid chromosomes in {fasta}")

	# read the exon table
	read_exons = pd.read_csv(exons, sep='\t', header=None)
	read_exons = read_exons.rename(columns={0: "chr", 1: "start", 2: "end", 3: "name", 4: "score", 5: "strand"})
	# read the DCE table
	read_dces = pd.read_csv(dces, sep='\t', header=None)
	read_dces = read_dces.rename(columns={0: "chr", 1: "start", 2: "end", 3: "name", 4: "max_overlap", 5: "strand"})
	# calculate the length of regions and exons
	read_dces['length'] = read_dces['end'] - read_dces['start']
	read_exons['length'] = read_exons['end'] - read_exons['start']

	# filter exons to only chromosomes present in the fasta
	n_before = len(read_exons)
	read_exons = read_exons[read_exons['chr'].isin(valid_chroms)].reset_index(drop=True)
	n_after = len(read_exons)
	if n_before != n_after:
		print(f"[info] Filtered out {n_before - n_after} exons on chromosomes not in {fasta} ({n_after} remaining)")

	# pre-filter DCE lengths to >= 21 nt (avoids rejection sampling in the loop)
	dce_lengths = read_dces['length'].values
	dce_lengths = dce_lengths[dce_lengths >= 21]

	# sort exons by length and extract columns as numpy arrays for fast access
	sort_idx = np.argsort(read_exons['length'].values)
	exon_chr = read_exons['chr'].values[sort_idx]
	exon_start = read_exons['start'].values[sort_idx]
	exon_end = read_exons['end'].values[sort_idx]
	exon_name = read_exons['name'].values[sort_idx]
	exon_score = read_exons['score'].values[sort_idx]
	exon_strand = read_exons['strand'].values[sort_idx]
	exon_length = read_exons['length'].values[sort_idx]
	n_exons = len(exon_length)

	# create an empty list to host the chosen intervals
	chosen_ints = []
	n = int(args.n_intervals)
	# here start the main loop:
	with alive_bar(n, bar = 'circles', spinner = 'classic') as bar:
		for i in range(n):
			# pick a region length (already >= 21 nt)
			reg_len = int(random.choice(dce_lengths))
			# find the first exon long enough using binary search on sorted lengths
			min_idx = np.searchsorted(exon_length, reg_len, side='left')
			# pick a random exon from the eligible range [min_idx, n_exons)
			idx = random.randint(min_idx, n_exons - 1)
			# define the range where we can safely pick a start
			safe_start = exon_end[idx] - reg_len
			new_start = random.randint(exon_start[idx], safe_start)
			new_end = new_start + reg_len
			chosen_ints.append([exon_chr[idx], new_start, new_end, exon_name[idx], exon_score[idx], exon_strand[idx]])
			bar()
		# turn into df
		output_df = pd.DataFrame(chosen_ints, columns = ['chr', 'start', 'end', 'name', 'score', 'strand'])
		output_df.to_csv("output/02_Conservation/fake_DCE_intervals.bed", sep = "\t", index=False)

# execute code
generate_intervals()