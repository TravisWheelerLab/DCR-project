#!/usr/bin/env python3.11

import pybedtools
import random
from alive_progress import alive_bar

# inputs
input_bed = open("data/mm39_exons_coding_regions_phased_4bpMin_noRandom.bed", 'r')
fasta = "data/mm39.fa"

# set random seed for reproducibility
random.seed(42)

# run getfasta one by one, save DNA in string
total = input_bed.readlines() # check file length for progress bar
output = [] # create empty list to host the output
with alive_bar(len(total), bar = 'smooth', spinner = 'dots_waves') as bar:
	for interval in total:
		# check if we have at least one codon
		start = interval.split("\t")[1]
		stop  = interval.split("\t")[2]
		exon_name = interval.split("\t")[3]
		if int(stop) - int(start) > 2:
			bedline = pybedtools.BedTool(interval, from_string=True) # make the line a bed object
			a = bedline.sequence(fi=fasta, s=True) # get the fasta record
			a_string = open(a.seqfn).read() # read it as a string
			a_list = [a_string.split("\n")] # make it a list
			seq = a_list[0][1]  # get the sequence from the list
			split_seq = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3] # split into codons, drop 2 or 1 bp leftovers
			# remove STOP codons
			no_stop_split_seq = [x for x in split_seq if not x.startswith('TAA') and not x.startswith('TAG') and not x.startswith('TGA')] # remove stops
			# shuffle the codons
			shuf_split_seq = random.sample(no_stop_split_seq, len(no_stop_split_seq))
			# accolate the codons into a sequence
			shuf_seq = ''.join(shuf_split_seq)
			# create the fasta entry
			out_line = '>' + exon_name + "_shuffled" + "\n" + str(shuf_seq)
			# append to output
			output.append(out_line)
			bar()
# write output
outfile = open("output/02_Conservation/fake_mm39_exons.fa", "w")
for line in output:
	outfile.write(str(line) + "\n")
outfile.close()
