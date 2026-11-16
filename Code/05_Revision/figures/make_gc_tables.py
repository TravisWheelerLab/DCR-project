#!/usr/bin/env python3
"""make_gc_tables.py -- regenerate per-exon / per-region GC tables for the GC-controls
supplementary figure from the executed DOOVER outputs.

  <datadir>/gc_per_exon_human.tsv  (no header) cls, chrom, start, end, GC%, length
      bedtools nuc on the three human exon beds of output/02_Conservation:
      DCR_exons_DCR_genes_UNIQ.bed (D), nonDCR_exon_DCR_gene_UNIQ.bed (D*),
      nonDCR_exon_nonDCR_gene_UNIQ.bed (ND)
  <datadir>/gc_per_region_mouse.tsv  cls, name, gc, len
      GC of DCE_overlapping_regions_LO.to.mm39.fa (mouse DCR orthologs) and
      fake_DCE_intervals.fa (random mouse control), GC = (G+C)/(A+C+G+T)

Usage: python make_gc_tables.py [--doover DIR] [--hg38 hg38.fa] [--datadir DIR]   (bedtools on PATH)
Reproduces gc_audit_summary.csv medians exactly (checked in the R figure script).
Runtime ~5 s with a local hg38.fa. Manuscript: Supplementary Figure S8 panels A-B.
"""
import os, subprocess, sys
import argparse as _argparse
_doover = os.environ.get("DOOVER", ".")
_ap = _argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
_ap.add_argument("--doover", default=_doover, help="DOOVER executed-pipeline directory (default $DOOVER)")
_ap.add_argument("--hg38", default=os.environ.get("HG38", os.path.join(os.environ.get("REFDIR", os.path.join(_doover, "data")), "hg38.fa")))
_ap.add_argument("--datadir", default=os.path.join(os.environ.get("WORK", "work"), "figure_data"),
                 help="output directory for the derived tables (the former data/)")
_args = _ap.parse_args()

HG38 = _args.hg38
C = os.path.join(_args.doover, "output", "02_Conservation")
os.makedirs(_args.datadir, exist_ok=True)

with open(os.path.join(_args.datadir, "gc_per_exon_human.tsv"), "w") as out:
    for cls, f in [("DCR_exon_human", "DCR_exons_DCR_genes_UNIQ.bed"),
                   ("nonDCR_exon_DCR_gene", "nonDCR_exon_DCR_gene_UNIQ.bed"),
                   ("nonDCR_exon_nonDCR_gene", "nonDCR_exon_nonDCR_gene_UNIQ.bed")]:
        bed3 = subprocess.run(["cut", "-f1-3", f"{C}/{f}"], capture_output=True, text=True,
                              check=True).stdout
        nuc = subprocess.run(["bedtools", "nuc", "-fi", HG38, "-bed", "-"], input=bed3,
                             capture_output=True, text=True, check=True).stdout
        for line in nuc.splitlines()[1:]:
            p = line.split("\t")      # chrom start end pct_at pct_gc A C G T N other len
            out.write(f"{cls}\t{p[0]}\t{p[1]}\t{p[2]}\t{float(p[4]) * 100:g}\t{p[11]}\n")


def fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip().upper())
    if name is not None:
        yield name, "".join(buf)


def gc(s):
    acgt = sum(s.count(c) for c in "ACGT")
    return 100 * (s.count("G") + s.count("C")) / acgt if acgt else float("nan")


with open(os.path.join(_args.datadir, "gc_per_region_mouse.tsv"), "w") as out:
    out.write("cls\tname\tgc\tlen\n")
    for cls, f in [("mouse_DCR_ortholog", "DCE_overlapping_regions_LO.to.mm39.fa"),
                   ("random_ctrl_mouse", "fake_DCE_intervals.fa")]:
        for n, s in fasta(f"{C}/{f}"):
            out.write(f"{cls}\t{n}\t{gc(s)}\t{len(s)}\n")
