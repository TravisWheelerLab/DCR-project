#!/usr/bin/env python3
"""prep_gene_models.py -- build the CSV inputs for the GNAS and GABBR1 gene-model
figures from the executed DOOVER pipeline outputs (protein-based DCR calls).

Frames are derived ONLY from the SwissProt isoform models mapped by MIRAGE2
(DOOVER BED12, whole block span = coding), never from GENCODE CDS annotation.

Frame label convention (reproduces Figure 2E, MADD: canonical '+2', alternative '0'):
  c = codon position (0,1,2) of the first DCR base in transcript direction,
      counted from the first coded base of the isoform model;
  frame offset = (3 - c) % 3 = offset (nt) of the first complete codon inside the
      DCR interval, in transcript direction.

Usage: python prep_gene_models.py [--doover DIR] [--hg38 hg38.fa] [--datadir DIR]
samtools must be on PATH. Runtime < 1 min.
Manuscript: Supplementary Figures S6 (GABBR1) and S7 (GNAS).

Outputs (--datadir, default $WORK/figure_data):
  gm_<gene>_exons.csv   isoform blocks (0-based half-open) + DCR exon frame offset/role
  gm_<gene>_dcr.csv     DCR interval
  gm_<gene>_iso.csv     isoform metadata from New_TableS1 (accession, evidence, ECO)
  gnas_mouse_frames.csv stop-codon positions per forward frame of the mouse GNAS ortholog
  gnas_mouse_summary.csv ORF lengths + 4-mer peptide homology of each mouse frame to the
                        human canonical (ALEX) / alternative (XLas) frame
"""
import os, subprocess, sys
import pandas as pd
import argparse as _argparse
_doover = os.environ.get("DOOVER", ".")
_ap = _argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
_ap.add_argument("--doover", default=_doover, help="DOOVER executed-pipeline directory (default $DOOVER)")
_ap.add_argument("--hg38", default=os.environ.get("HG38", os.path.join(os.environ.get("REFDIR", os.path.join(_doover, "data")), "hg38.fa")))
_ap.add_argument("--datadir", default=os.path.join(os.environ.get("WORK", "work"), "figure_data"),
                 help="output directory for the derived tables (the former data/)")
_args = _ap.parse_args()

DOOVER = _args.doover
HG38 = _args.hg38
OUT = _args.datadir
os.makedirs(OUT, exist_ok=True)

BED = f"{DOOVER}/output/01_Description/FINAL_Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguousGenes_noAmbiguousEnds.extBed.bed"
DCRBED = f"{DOOVER}/output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed"
TS1 = f"{DOOVER}/output/New_TableS1.tsv"
MM_BED = f"{DOOVER}/output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed"
MM_FA = f"{DOOVER}/output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.fa"

iso = pd.read_csv(BED, sep="\t", header=None, usecols=range(12))
iso.columns = "chrom start end iname score strand ts te rgb nb sizes starts".split()
dcr = pd.read_csv(DCRBED, sep="\t", header=None, names="chrom start end id strand gene nf".split())
ts1 = pd.read_csv(TS1, sep="\t", header=None).iloc[:, :10]
ts1.columns = ["iname", "dcr_id", "frame_group", "accession", "evidence", "uniprot_status",
               "canonical_group", "frame_role", "eco", "nmd"]


# isoforms absent from New_TableS1 (no DCR): SwissProt accession by entry name
# (GNAS3_HUMAN = O95467, NESP55; checked against rest.uniprot.org)
EXTRA_ACC = {"homo_sapiens|gnas/gnas3_human|1": "O95467-1"}


def blocks(r):
    s = [int(x) for x in str(r.sizes).strip(",").split(",")]
    o = [int(x) for x in str(r.starts).strip(",").split(",")]
    return [(r.start + a, r.start + a + b) for a, b in zip(o, s)]


def codon_pos(r, g):
    bl = blocks(r)
    i = 0
    for a, b in (bl if r.strand == "+" else list(reversed(bl))):
        if a <= g < b:
            return (i + (g - a if r.strand == "+" else b - 1 - g)) % 3
        i += b - a
    return None


def build(gene_key, dcr_id, label):
    d = dcr[dcr.id == dcr_id].iloc[0]
    g0 = d.start if d.strand == "+" else d.end - 1
    sub = iso[iso.iname.str.contains(f"\\|{gene_key}/", regex=True)].copy()
    meta = ts1[ts1.dcr_id == dcr_id].set_index("iname")
    rows = []
    for _, r in sub.iterrows():
        acc = meta.accession.get(r.iname, None)
        if acc is None:   # isoform not part of this DCR: accession from any DCR row of Table S1
            hit = ts1.loc[ts1.iname == r.iname, "accession"]
            acc = hit.iloc[0] if len(hit) else EXTRA_ACC.get(r.iname,
                  r.iname.split("|")[1].split("/")[1].upper() + "-" + r.iname.split("|")[2])
        role = meta.frame_role.get(r.iname, None)
        c = codon_pos(r, g0)
        off = None if c is None else (3 - c) % 3
        for a, b in blocks(r):
            hit = (a < d.end and b > d.start) and role is not None
            rows.append(dict(isoform=r.iname, accession=acc, chrom=r.chrom, strand=r.strand,
                             start=a, end=b, dcr_exon=hit,
                             frame_offset=off if hit else None, frame_role=role if hit else None,
                             iso_start=r.start, iso_end=r.end, in_dcr=role is not None))
    ex = pd.DataFrame(rows)
    ex.to_csv(f"{OUT}/gm_{label}_exons.csv", index=False)
    pd.DataFrame([d]).to_csv(f"{OUT}/gm_{label}_dcr.csv", index=False)
    meta.reset_index().to_csv(f"{OUT}/gm_{label}_iso.csv", index=False)
    # our frame partition must equal the pipeline's frame groups
    chk = ex[ex.dcr_exon].drop_duplicates("isoform").merge(meta.reset_index(), left_on="isoform", right_on="iname")
    po = sorted(sorted(s) for s in chk.groupby("frame_offset").isoform.apply(list))
    pp = sorted(sorted(s) for s in chk.groupby("frame_group").isoform.apply(list))
    assert po == pp, (label, po, pp)
    print(label, chk[["accession_x", "frame_offset", "frame_role_x"]].to_dict("records"))
    return d


build("gabbr1", 1812, "gabbr1")
d201 = build("gnas", 201, "gnas")

# ---------------- mouse GNAS ortholog of DCR 201 ----------------
_B = "TCAG"
CODE = dict(zip([a + b + c for a in _B for b in _B for c in _B],
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def translate(s):
    return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


mmb = pd.read_csv(MM_BED, sep="\t", header=None, names="chrom start end name score strand".split())
row = mmb[mmb.name.str.startswith("DCRid=201;")].iloc[0]
key = f"{row.chrom}:{row.start}-{row.end}({row.strand})"
fa, name = {}, None
for line in open(MM_FA):
    if line.startswith(">"):
        name = line[1:].strip(); fa[name] = []
    else:
        fa[name].append(line.strip().upper())
mseq = "".join(fa[key])
hseq = subprocess.run(["samtools", "faidx", HG38, f"{d201.chrom}:{d201.start + 1}-{d201.end}"],
                      capture_output=True, text=True, check=True).stdout.split("\n", 1)[1].replace("\n", "").upper()
assert d201.strand == "+"
ex = pd.read_csv(f"{OUT}/gm_gnas_exons.csv")
hum = ex[ex.dcr_exon].drop_duplicates("isoform")
hframes = {}
for role, grp in hum.groupby("frame_role"):
    off = int(grp.frame_offset.iloc[0])
    hframes[role] = (off, translate(hseq[off:]), ",".join(sorted(grp.accession.unique())))


def kmers(p, k=4):
    return {p[i:i + k] for i in range(len(p) - k + 1) if "*" not in p[i:i + k]}


frows, srows = [], []
for f in (1, 2, 3):
    off = f - 1
    pep = translate(mseq[off:])
    stops = [off + 3 * i for i, a in enumerate(pep) if a == "*"]
    frows += [dict(frame=f, pos=s, type="stop") for s in stops]
    km = kmers(pep)
    share = {role: len(km & kmers(hp)) / max(1, len(kmers(hp))) for role, (_, hp, _) in hframes.items()}
    srows.append(dict(frame=f, offset=off, n_codons=len(pep), n_stops=len(stops),
                      full_length_orf_aa=(len(pep) if not stops else None),
                      longest_orf_aa=max(len(x) for x in pep.split("*")),
                      kmer_share_canonical=round(share.get("canonical", 0), 3),
                      kmer_share_alternative=round(share.get("alternative", 0), 3),
                      mouse_chrom=row.chrom, mouse_start=row.start, mouse_end=row.end,
                      mouse_strand=row.strand, mouse_len=len(mseq), human_len=len(hseq)))
pd.DataFrame(frows, columns=["frame", "pos", "type"]).to_csv(f"{OUT}/gnas_mouse_frames.csv", index=False)
summ = pd.DataFrame(srows)
for role, (off, hp, acc) in hframes.items():
    summ[f"human_{role}_offset"] = off
    summ[f"human_{role}_accessions"] = acc
    summ[f"human_{role}_stops"] = hp.count("*")
summ.to_csv(f"{OUT}/gnas_mouse_summary.csv", index=False)
print(summ.T.to_string())

# ---------------- GABBR1 DCR sequence + both translations (transcript orientation) --------
def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


d = dcr[dcr.id == 1812].iloc[0]
raw = subprocess.run(["samtools", "faidx", HG38, f"{d.chrom}:{d.start + 1}-{d.end}"],
                     capture_output=True, text=True, check=True).stdout.split("\n", 1)[1].replace("\n", "").upper()
tx = revcomp(raw) if d.strand == "-" else raw
ex = pd.read_csv(f"{OUT}/gm_gabbr1_exons.csv")
fr = ex[ex.dcr_exon].drop_duplicates("isoform").groupby("frame_role").agg(
    offset=("frame_offset", "first"), accessions=("accession", lambda x: ",".join(sorted(set(x))))).reset_index()
rows = []
for i, b in enumerate(tx):
    gpos = (d.end - i) if d.strand == "-" else (d.start + 1 + i)   # 1-based genomic
    rows.append(dict(track="nt", i=i, char=b, gpos=gpos, role=None, codon=None))
for _, f in fr.iterrows():
    off = int(f.offset)
    for k, j in enumerate(range(off, len(tx) - 2, 3)):
        rows.append(dict(track="aa", i=j, char=CODE[tx[j:j + 3]], gpos=None, role=f.frame_role, codon=tx[j:j + 3]))
seq = pd.DataFrame(rows)
assert "*" not in "".join(seq[seq.track == "aa"].char), "stop codon inside DCR"
seq.to_csv(f"{OUT}/gm_gabbr1_dcr_seq.csv", index=False)
fr.to_csv(f"{OUT}/gm_gabbr1_frames.csv", index=False)
print(tx, {r.frame_role: "".join(seq[(seq.track == 'aa') & (seq.role == r.frame_role)].char) for _, r in fr.iterrows()})
