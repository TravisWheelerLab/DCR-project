#!/usr/bin/env python3
"""
A2 step 2 - assemble the comparison arms and annotate length, frame preservation, GC and mean phyloP100way.
Arms
  DCR                : 1,407 catalogue DCRs (1,741 intervals; filtered_DCR_overlapping_regions_with_gene_and_dceID.bed)
  regAS_DCRgenes     : regular_AS_internal segments in the 1,296 catalogue DCR genes
  regAS_nonDCRgenes  : regular_AS_internal segments in genes with NO dual-coding anywhere (no derived
                       dual-coding atom, no raw MIRAGE DCR call in homo_sapiens.tsv, not in catalogue)
  const_DCRgenes / const_nonDCRgenes : constitutive coding segments (context only)
Genes with derived/raw dual coding that did not reach the final catalogue ('grey-zone') are excluded from all
regular-AS arms.
Paths   in/  = --indir  (default $DOOVER/output/01_Description)
        out/ = --outdir (default $WORK/splicing_vs_AS)
Manuscript  Supplementary Figure S10 (DCRs vs ordinary alternative splicing) and its Results/Methods
            text; full write-up in docs/A2_results.md.
"""

import argparse as _argparse, os as _os
_doover = _os.environ.get("DOOVER", ".")
_refdir = _os.environ.get("REFDIR", _os.path.join(_doover, "data"))
_ap = _argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
_ap.add_argument("--indir", default=_os.path.join(_doover, "output", "01_Description"),
                 help="DOOVER output/01_Description (pipeline outputs read by this step)")
_ap.add_argument("--outdir", default=_os.path.join(_os.environ.get("WORK", "work"), "splicing_vs_AS"),
                 help="A2 working/output directory (the former out/)")
_ap.add_argument("--hg38", default=_os.environ.get("HG38", _os.path.join(_refdir, "hg38.fa")))
_ap.add_argument("--phylop", default=_os.environ.get("PHYLOP_BW", _os.path.join(_refdir, "hg38.phyloP100way.bw")))
_args = _ap.parse_args()
IN = lambda f: _os.path.join(_args.indir, f)
OUT = lambda f: _os.path.join(_args.outdir, f)
_os.makedirs(_args.outdir, exist_ok=True)
import pandas as pd, numpy as np, pyBigWig
from pyfaidx import Fasta
MINLEN = 6   # = shortest catalogue DCR interval
S = pd.read_csv(OUT("segments.tsv"), sep="\t", keep_default_na=False)
dcr = pd.read_csv(IN("filtered_DCR_overlapping_regions_with_gene_and_dceID.bed"), sep="\t", header=None,
                  names=["chr", "start", "end", "DCR", "strand", "gene", "frames"])
raw = pd.read_csv(IN("homo_sapiens.tsv"), sep="\t")
cat_genes = set(dcr.gene)
dual_genes = set(S.gene[S.cls == "dual_coding"])
raw_genes = set(raw.Gene.str.lower())
non_dcr_genes = set(S.gene) - cat_genes - dual_genes - raw_genes
grey = (dual_genes | raw_genes) - cat_genes

def arm(r):
    if r.gene in cat_genes: tag = "DCRgenes"
    elif r.gene in non_dcr_genes: tag = "nonDCRgenes"
    else: return None
    if r.cls == "regular_AS_internal": return "regAS_" + tag
    if r.cls == "constitutive": return "const_" + tag
    return None
S["arm"] = S.apply(arm, axis=1)
E = S[S.arm.notna() & (S.length >= MINLEN)].copy()
E["event_id"] = E.gene + ":" + E.chr + ":" + E.start.astype(str) + "-" + E.end.astype(str)
E["intervals"] = E.chr + ":" + E.start.astype(str) + "-" + E.end.astype(str)

# DCR arm, one row per DCR (intervals pooled)
d = dcr.groupby("DCR").agg(gene=("gene", "first"), chr=("chr", "first"), strand=("strand", "first"),
                          start=("start", "min"), end=("end", "max"),
                          length=("end", lambda x: 0), n_int=("start", "size")).reset_index()
d["length"] = dcr.assign(L=dcr.end - dcr.start).groupby("DCR").L.sum().values
d["intervals"] = dcr.groupby("DCR").apply(lambda t: ";".join(f"{c}:{a}-{b}" for c, a, b in zip(t.chr, t.start, t.end))).values
d["arm"] = "DCR"; d["event_id"] = "DCR" + d.DCR.astype(str)
cols = ["event_id", "arm", "gene", "chr", "strand", "start", "end", "length", "intervals"]
E = pd.concat([d[cols + ["DCR"]], E[cols + ["cov_isoforms", "spanning_excluders", "n_iso_gene", "n_cov",
                                            "all_cov_iso_clean", "all_excl_iso_clean"]]], ignore_index=True)
E["mult3"] = (E.length % 3 == 0)

bw = pyBigWig.open(_args.phylop); fa = Fasta(_args.hg38, sequence_always_upper=True)
chroms = bw.chroms()
pp, gc, ppcov = [], [], []
for iv in E.intervals:
    tot = 0.0; n = 0; g_ = 0; L = 0
    for x in iv.split(";"):
        c, ab = x.rsplit(":", 1); a, b = map(int, ab.split("-"))
        s = str(fa[c][a:b]); L += len(s); g_ += s.count("G") + s.count("C")
        if c in chroms:
            v = np.array(bw.values(c, a, b), dtype=float); v = v[~np.isnan(v)]
            tot += v.sum(); n += len(v)
    pp.append(tot / n if n else np.nan); ppcov.append(n / L); gc.append(g_ / L)
E["phyloP_mean"] = pp; E["phyloP_bp_frac"] = ppcov; E["GC"] = gc
E.to_csv(OUT("A2_events.tsv"), sep="\t", index=False)
print(E.arm.value_counts().to_string())
print("genes: catalogue", len(cat_genes), "nonDCR", len(non_dcr_genes), "grey-zone excluded", len(grey))
print(E.groupby("arm").gene.nunique().to_string())
