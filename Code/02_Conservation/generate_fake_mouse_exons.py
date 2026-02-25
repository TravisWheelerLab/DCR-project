#!/usr/bin/env python3.11

import pybedtools
import random

# inputs
input_bed = "output/02_Conservation/mm39_CDS_codon_phased_21bpMin_stdChrom_500less.bed"
fasta = "data/mm39.fa"

# set random seed for reproducibility
random.seed(42)

# Filter intervals with at least one codon (length > 2) and extract all sequences at once
bed = pybedtools.BedTool(input_bed).filter(lambda x: len(x) > 2).saveas()
seqs = bed.sequence(fi=fasta, s=True)

# Collect exon names from the filtered BED (preserving order)
exon_names = [interval.name for interval in bed]

# Parse all sequences from the single FASTA output (order matches BED)
sequences = []
with open(seqs.seqfn) as f:
	for line in f:
		line = line.strip()
		if not line.startswith('>'):
			sequences.append(line)

# Process each exon: split into codons, remove stops, shuffle, write output
stop_codons = {'TAA', 'TAG', 'TGA', 'taa', 'tag', 'tga'}
with open("output/02_Conservation/fake_mm39_exons.fa", "w") as outfile:
	for exon_name, seq in zip(exon_names, sequences):
		# split into codons, drop incomplete trailing codons
		codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
		# remove STOP codons
		codons = [c for c in codons if c.upper() not in stop_codons]
		# shuffle the codons
		codons = random.sample(codons, len(codons))
		# write the fasta entry
		outfile.write('>' + exon_name + "_shuffled\n" + ''.join(codons) + "\n")
