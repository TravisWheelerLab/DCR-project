import random
from collections import defaultdict
import sys

# Input and output filenames
input_bed = "output/01_Description/UniprotKB_isoforms_mapping.bed"
output_bed = "output/01_Description/EXPERIMENTAL_control_splice_sites.bed"

# Optional reproducibility
random.seed(42)

def parse_gene(name):
    return name.split("|")[1]

def parse_isoform(name):
    return name

genes = defaultdict(lambda: defaultdict(list))

# Load BED input
with open(input_bed) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 6:
            continue
        chrom, start, end, name, score, strand = parts[:6]
        start, end = int(start), int(end)
        gene = parse_gene(name)
        isoform = parse_isoform(name)
        genes[gene][isoform].append((start, end, chrom, strand))

# Store output with deduplication
output_coords = defaultdict(set)

for gene, isoforms in genes.items():
    # Get all unique exons across isoforms
    all_exons = set()
    for exons in isoforms.values():
        all_exons.update((start, end) for start, end, *_ in exons)

    if len(all_exons) <= 1:
        continue  # Skip single-exon genes

    picked_exon = random.choice(list(all_exons))

    # Find isoforms using the picked exon
    matching_isoforms = [
        (isoform, sorted(exons, key=lambda x: x[0] if x[3] == '+' else -x[1]))
        for isoform, exons in isoforms.items()
        if picked_exon in [(start, end) for start, end, *_ in exons]
    ]

    for isoform, exons in matching_isoforms:
        exon_coords = [(start, end) for start, end, *_ in exons]
        chrom = exons[0][2]
        strand = exons[0][3]

        try:
            idx = exon_coords.index(picked_exon)
        except ValueError:
            continue

        def add_region(start, end, junction, site):
            if end <= start:
                return
            annotation = f"{gene}|junction={junction};site={site}"
            output_coords[gene].add((chrom, start, end, annotation, strand))

        ### Upstream intron ###
        if idx > 0:
            up_start, up_end = exon_coords[idx - 1]
            pk_start, pk_end = picked_exon

            if strand == '+':
                # Donor = end of upstream exon
                add_region(up_end, up_end + 10, "upstream", "donor")
                # Acceptor = start of picked exon
                add_region(pk_start - 10, pk_start, "upstream", "acceptor")
            else:
                # Donor = start of upstream exon
                add_region(up_start - 10, up_start, "upstream", "donor")
                # Acceptor = end of picked exon
                add_region(pk_end, pk_end + 10, "upstream", "acceptor")

        ### Downstream intron ###
        if idx < len(exon_coords) - 1:
            pk_start, pk_end = picked_exon
            down_start, down_end = exon_coords[idx + 1]

            if strand == '+':
                # Donor = end of picked exon
                add_region(pk_end, pk_end + 10, "downstream", "donor")
                # Acceptor = start of downstream exon
                add_region(down_start - 10, down_start, "downstream", "acceptor")
            else:
                # Donor = start of picked exon
                add_region(pk_start - 10, pk_start, "downstream", "donor")
                # Acceptor = end of downstream exon
                add_region(down_end, down_end + 10, "downstream", "acceptor")

# Write output BED
with open(output_bed, "w") as out:
    for gene, coords in output_coords.items():
        for chrom, start, end, annotation, strand in sorted(coords):
            out.write(f"{chrom}\t{start}\t{end}\t{annotation}\t.\t{strand}\n")