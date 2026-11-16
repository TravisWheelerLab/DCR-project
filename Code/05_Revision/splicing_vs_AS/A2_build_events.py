#!/usr/bin/env python3
"""
A2 step 1 - build atomic coding segments from the SwissProt isoform models (MIRAGE2 mappings)
and classify them as dual-coding / regular alternative splicing / constitutive.

Input : --indir/UniprotKB_isoforms_mapping.bed (DOOVER output/01_Description; protein-derived CDS exons,
        each isoform model starts at its first codon and excludes the stop codon)
        --indir/isoTrack_ambiguous_gene_mapping.txt  (genes removed by the pipeline)
        --indir/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed (final DCR catalogue intervals)
        --hg38 hg38.fa (local copy)
Frames: the codon position of every coding base is computed from the protein model itself
        (cumulative CDS length from the first codon, offset 0). No GENCODE annotation is used.
Output: --outdir/atoms.tsv.gz       one row per atomic interval x gene
        --outdir/segments.tsv       merged segments (maximal runs of adjacent atoms with identical
                               isoform-coverage set and class)
        --outdir/isoform_frames.tsv per-isoform QC (internal stops at offset 0)
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
_args = _ap.parse_args()
IN = lambda f: _os.path.join(_args.indir, f)
OUT = lambda f: _os.path.join(_args.outdir, f)
_os.makedirs(_args.outdir, exist_ok=True)
import pandas as pd, numpy as np, collections, gzip
from pyfaidx import Fasta

iso = pd.read_csv(IN("UniprotKB_isoforms_mapping.bed"), sep="\t", header=None,
                  names=["chr", "start", "end", "isoform", "score", "strand", "phase"])
amb = set(open(IN("isoTrack_ambiguous_gene_mapping.txt")).read().split())
iso["gene"] = iso.isoform.str.split(r"[|/]").str[1]
iso = iso[~iso.gene.isin(amb)].copy()

# ---- per-isoform QC: internal stops when read from offset 0 (protein-model frame)
g = Fasta(_args.hg38, sequence_always_upper=True)
comp = str.maketrans("ACGTN", "TGCAN")
STOP = {"TAA", "TAG", "TGA"}
qc = []
for name, t in iso.groupby("isoform", sort=False):
    st = t.strand.iloc[0]; t = t.sort_values("start"); c = t.chr.iloc[0]
    seq = "".join(str(g[c][a:b]) for a, b in zip(t.start, t.end))
    if st == "-": seq = seq.translate(comp)[::-1]
    n = sum(1 for i in range(0, len(seq) - 2, 3) if seq[i:i + 3] in STOP)
    qc.append((name, t.gene.iloc[0], st, len(seq), len(t), n, seq[:3] == "ATG"))
qc = pd.DataFrame(qc, columns=["isoform", "gene", "strand", "cds_len", "n_exons", "internal_stops_off0", "starts_ATG"])
qc.to_csv(OUT("isoform_frames.tsv"), sep="\t", index=False)
iso_ok = set(qc.isoform[qc.internal_stops_off0 == 0])

# ---- per-isoform exon list with transcript offset at exon 5' end
iso = iso.sort_values(["isoform", "start"])
ex = []  # (gene, chr, strand, isoform, start, end, c0)
for name, t in iso.groupby("isoform", sort=False):
    st = t.strand.iloc[0]
    rows = list(zip(t.start, t.end))
    if st == "-": rows = rows[::-1]
    c0 = 0
    for a, b in rows:
        ex.append((t.gene.iloc[0], t.chr.iloc[0], st, name, a, b, c0)); c0 += b - a
ex = pd.DataFrame(ex, columns=["gene", "chr", "strand", "isoform", "start", "end", "c0"])

def cp_at(row, pos):
    """codon position (0,1,2) of genomic base pos in this isoform exon"""
    if row.strand == "+": return (row.c0 + pos - row.start) % 3
    return (row.c0 + (row.end - 1 - pos)) % 3

atoms = []
for gene, t in ex.groupby("gene", sort=False):
    if t.chr.nunique() > 1 or t.strand.nunique() > 1: continue
    isos = sorted(t.isoform.unique())
    span = {i: (s.start.min(), s.end.max()) for i, s in t.groupby("isoform")}
    bps = sorted(set(t.start) | set(t.end))
    for a, b in zip(bps[:-1], bps[1:]):
        cov = t[(t.start <= a) & (t.end >= b)]
        if len(cov) == 0: continue  # intron of every isoform
        fr = {r.isoform: cp_at(r, a) for r in cov.itertuples()}
        covset = tuple(sorted(fr))
        excl = [i for i in isos if i not in fr]
        spanning_excl = [i for i in excl if span[i][0] < a and span[i][1] > b]
        nframes = len(set(fr.values()))
        atoms.append((gene, t.chr.iloc[0], t.strand.iloc[0], a, b, len(isos), len(covset), nframes,
                      ",".join(covset), ",".join(f"{fr[i]}" for i in covset), ",".join(spanning_excl), len(excl)))
A = pd.DataFrame(atoms, columns=["gene", "chr", "strand", "start", "end", "n_iso_gene", "n_cov", "n_frames",
                                 "cov_isoforms", "cov_cp", "spanning_excluders", "n_excl"])
A["n_span_excl"] = A.spanning_excluders.apply(lambda s: 0 if s == "" else len(s.split(",")))

def cls(r):
    if r.n_frames > 1: return "dual_coding"
    if r.n_excl == 0: return "constitutive"
    if r.n_span_excl > 0: return "regular_AS_internal"
    return "variable_terminal"   # absent only from isoforms that start/end before it (alt. termini, truncations)
A["class"] = A.apply(cls, axis=1)
A.to_csv(OUT("atoms.tsv.gz"), sep="\t", index=False, compression="gzip")

# ---- merge adjacent atoms with identical coverage set + class (+ for regular AS, identical spanning-excluder set)
segs = []
for gene, t in A.groupby("gene", sort=False):
    t = t.sort_values("start"); cur = None
    for r in t.itertuples():
        key = (r[-1], r.cov_isoforms, r.spanning_excluders)   # r[-1] = class
        if cur is not None and r.start == cur["end"] and key == cur["key"]:
            cur["end"] = r.end; cur["n_atoms"] += 1; cur["cps"].append(r.cov_cp)
        else:
            if cur is not None: segs.append(cur)
            cur = dict(gene=gene, chr=r.chr, strand=r.strand, start=r.start, end=r.end, key=key, n_atoms=1,
                       cls=r[-1], cov_isoforms=r.cov_isoforms, spanning_excluders=r.spanning_excluders,
                       n_iso_gene=r.n_iso_gene, n_cov=r.n_cov, n_frames=r.n_frames, cps=[r.cov_cp])
    if cur is not None: segs.append(cur)
S = pd.DataFrame(segs).drop(columns=["key", "cps"])
S["length"] = S.end - S.start
S["all_cov_iso_clean"] = S.cov_isoforms.apply(lambda s: all(i in iso_ok for i in s.split(",")))
S["all_excl_iso_clean"] = S.spanning_excluders.apply(lambda s: s == "" or all(i in iso_ok for i in s.split(",")))
S.to_csv(OUT("segments.tsv"), sep="\t", index=False)
print(A["class"].value_counts().to_string()); print(S.cls.value_counts().to_string())
