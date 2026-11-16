#!/usr/bin/env python3
"""verify_latent_exon_fraction.py -- exon-level check of the Figure 2B(i) fractions.

For each exon (strand-oriented sequence from hg38), find stop-free runs in each of the
3 forward frames (esl-translate --watson semantics: stop-to-stop, leading and trailing
partial runs included, coding nt only). The exon's "max dual-coding overlap" is the
longest genomic overlap between a run in one frame and a run in another frame. The
permissive sub-interval rule passes an exon when this overlap is >= 21 nt.

Validation target (definition_reconciliation.csv): DCR exons 98.0% (n = 3,324), median
max overlap 79 nt. Output: data/latent_exon_fraction.csv (per class: n, n_pass, pct,
median max overlap) -- tells whether N = 289,416 is an exon count.
Output path: <datadir>/latent_exon_fraction.csv.
Usage: python verify_latent_exon_fraction.py [--doover DIR] [--hg38 hg38.fa] [--datadir DIR]
Runtime ~20 s with a local hg38.fa (two passes over ~300k exons). Manuscript: Figure 2B(i) (check of the exon-level
fractions; not plotted).
"""
import os, subprocess, sys, statistics
import pandas as pd
import argparse as _argparse
_doover = os.environ.get("DOOVER", ".")
_ap = _argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
_ap.add_argument("--doover", default=_doover, help="DOOVER executed-pipeline directory (default $DOOVER)")
_ap.add_argument("--hg38", default=os.environ.get("HG38", os.path.join(os.environ.get("REFDIR", os.path.join(_doover, "data")), "hg38.fa")))
_ap.add_argument("--datadir", default=os.path.join(os.environ.get("WORK", "work"), "figure_data"),
                 help="output directory for the derived tables (the former data/)")
_args = _ap.parse_args()

HG38 = _args.hg38
# The original analysis read exon sequences on the genomic + strand (bedtools getfasta
# without -s); that variant reproduces the published 98.0% / median 79 nt exactly.
# Both orientations are reported.
C = os.path.join(_args.doover, "output", "02_Conservation")
STOPS = {"TAA", "TAG", "TGA"}


def runs(seq, off):
    out, start, i = [], off, off
    while i + 3 <= len(seq):
        if seq[i:i + 3] in STOPS:
            if i > start:
                out.append((start, i))
            start = i + 3
        i += 3
    end = off + 3 * ((len(seq) - off) // 3)
    if end > start:
        out.append((start, end))
    return out


def max_overlap(seq):
    r = [runs(seq, f) for f in range(3)]
    best = 0
    for a in range(3):
        for b in range(a + 1, 3):
            for s1, e1 in r[a]:
                for s2, e2 in r[b]:
                    ov = min(e1, e2) - max(s1, s2)
                    if ov > best:
                        best = ov
    return best


rows = []
for orient in ("genomic_plus", "transcript"):
  for cls, f in [("DCR_exon_human", "DCR_exons_DCR_genes_UNIQ.bed"),
                 ("nonDCR_exon_DCR_gene", "nonDCR_exon_DCR_gene_UNIQ.bed"),
                 ("nonDCR_exon_nonDCR_gene", "nonDCR_exon_nonDCR_gene_UNIQ.bed")]:
      bed6 = []
      for line in open(f"{C}/{f}"):
          p = line.rstrip("\n").split("\t")
          strand = "+" if p[3].endswith("_f") else "-"
          bed6.append(f"{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\t0\t{strand}")
      fa = subprocess.run(["bedtools", "getfasta", "-fi", HG38, "-bed", "-", "-tab"] + (["-s"] if orient == "transcript" else []),
                          input="\n".join(bed6) + "\n", capture_output=True, text=True,
                          check=True).stdout.splitlines()
      mo = [max_overlap(l.split("\t")[1].upper()) for l in fa]
      n_pass = sum(m >= 21 for m in mo)
      rows.append(dict(orientation=orient, cls=cls, n=len(mo), n_pass=n_pass, pct=round(100 * n_pass / len(mo), 2),
                       median_max_overlap=statistics.median(mo)))
      print(rows[-1], flush=True)
df = pd.DataFrame(rows)
for orient, d in list(df.groupby("orientation")):
    non = d[d.cls != "DCR_exon_human"]
    df.loc[len(df)] = dict(orientation=orient, cls="nonDCR_all", n=non.n.sum(),
                           n_pass=non.n_pass.sum(),
                           pct=round(100 * non.n_pass.sum() / non.n.sum(), 2),
                           median_max_overlap=float("nan"))
os.makedirs(_args.datadir, exist_ok=True)
df.to_csv(os.path.join(_args.datadir, "latent_exon_fraction.csv"), index=False)
print(df.to_string())
