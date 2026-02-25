#!/usr/bin/env python3
"""
Extract codon-phased CDS regions from Ensembl GTF and output as BED.

For each CDS entry in the GTF:
  - Use the frame field (column 8) to trim incomplete codons at the 5' end
  - Trim any remaining incomplete codon at the 3' end
  - Filter to standard chromosomes (1-19, X, Y)
  - Filter to regions >= 21 bp
  - Add "chr" prefix for UCSC/mm39.fa compatibility

Output: 6-column BED file (chr, start, end, name, score, strand)
compatible with generate_fake_mouse_exons.py (codon shuffling).

GTF frame field:
  0 = feature starts at the first base of a codon
  1 = 1 base before the first complete codon (skip 1 from 5' end)
  2 = 2 bases before the first complete codon (skip 2 from 5' end)
"""

import gzip
import sys
import os

# --- Configuration ---
INPUT_GTF = "data/Mus_musculus.GRCm39.115.gtf.gz"
OUTPUT_BED = "output/02_Conservation/mm39_CDS_codon_phased_21bpMin_stdChrom.bed"
MIN_LENGTH = 21
STANDARD_CHROMS = set([str(i) for i in range(1, 20)] + ["X", "Y"])

# --- Main ---
n_total_cds = 0
n_kept = 0
n_filtered_chrom = 0
n_filtered_length = 0

with gzip.open(INPUT_GTF, "rt") as fin, open(OUTPUT_BED, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue

        chrom = fields[0]
        feature = fields[2]
        gtf_start = int(fields[3])  # 1-based inclusive
        gtf_end = int(fields[4])    # 1-based inclusive
        strand = fields[6]
        frame = fields[7]

        if feature != "CDS":
            continue

        n_total_cds += 1

        # Filter chromosomes
        if chrom not in STANDARD_CHROMS:
            n_filtered_chrom += 1
            continue

        # Frame must be valid
        if frame not in ("0", "1", "2"):
            continue

        frame = int(frame)

        # Convert GTF (1-based inclusive) to BED (0-based half-open)
        bed_start = gtf_start - 1
        bed_end = gtf_end

        # Adjust for codon phase
        if strand == "+":
            # Frame = bases to skip from 5' end (= genomic left on + strand)
            bed_start += frame
            # Trim 3' end (= genomic right on + strand) to complete codons
            length = bed_end - bed_start
            bed_end -= (length % 3)
        elif strand == "-":
            # Frame = bases to skip from 5' end (= genomic right on - strand)
            bed_end -= frame
            # Trim 3' end (= genomic left on - strand) to complete codons
            length = bed_end - bed_start
            bed_start += (length % 3)

        # Length filter
        final_length = bed_end - bed_start
        if final_length < MIN_LENGTH:
            n_filtered_length += 1
            continue

        # Sanity check: length must be divisible by 3
        assert final_length % 3 == 0, (
            f"Length not divisible by 3: {chrom}:{bed_start}-{bed_end} "
            f"({final_length} bp, strand={strand}, frame={frame})"
        )

        # Parse attributes for a descriptive name
        attrs = fields[8]
        gene_name = ""
        transcript_id = ""
        exon_number = ""
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_name"):
                gene_name = attr.split('"')[1]
            elif attr.startswith("transcript_id"):
                transcript_id = attr.split('"')[1]
            elif attr.startswith("exon_number"):
                exon_number = attr.split('"')[1]

        name = f"{gene_name}_{transcript_id}_exon{exon_number}"

        # Write BED line with "chr" prefix
        fout.write(f"chr{chrom}\t{bed_start}\t{bed_end}\t{name}\t0\t{strand}\n")
        n_kept += 1

print(f"Total CDS entries in GTF:        {n_total_cds}")
print(f"Filtered (non-standard chrom):   {n_filtered_chrom}")
print(f"Filtered (< {MIN_LENGTH} bp after phasing): {n_filtered_length}")
print(f"Kept:                            {n_kept}")
print(f"Output written to: {OUTPUT_BED}")
