#!/usr/bin/env python3
"""splice_site_phylop.py -- phyloP100way conservation at DCR vs control splice sites (reviewer R2-M4b).

Purpose  Mean phyloP100way over each splice-site interval written by the pipeline
         (01_Description/extract_DCE_splice_sites.sh and extract_control_exons_splice_site.py),
         then DCR vs control per junction (upstream/downstream) x site (donor/acceptor) with a
         two-sided Mann-Whitney U test, rank-biserial effect size and Bonferroni adjustment.
Inputs   --dcr-sites   $DOOVER/output/01_Description/SPLICE_SITES_filtered_DCE_splice_sites.bed
         --ctrl-sites  $DOOVER/output/01_Description/SPLICE_SITES_control_splice_sites.bed
         --phylop      hg38.phyloP100way.bw (default $PHYLOP_BW, else $REFDIR/hg38.phyloP100way.bw)
                       Copy the 9.9 GB bigWig to a local disk first: random reads over a network /
                       cloud-synced drive run at ~1.6 queries/s (~11 h for this step) versus
                       ~10 s from a local SSD.
Outputs  --outdir/splice_site_phylop.tsv        one row per site: chrom start end junction site cls phyloP
         --outdir/splice_site_conservation.csv  per junction x site tests (+ ALL)
Manuscript  Supplementary Figure S9 panel A (figures/supp_splice_expression_tissue.R); Results,
            splice-site conservation paragraph.
"""
import argparse, os, re, time
import numpy as np, pandas as pd, pyBigWig
from scipy import stats


def load(path, cls):
    rows = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4: continue
        j = re.search(r'junction=(\w+)', f[3]); s = re.search(r'site=(\w+)', f[3])
        rows.append((f[0], int(f[1]), int(f[2]), j.group(1) if j else None, s.group(1) if s else None, cls))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "junction", "site", "cls"])


def rb(a, b):
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return 1 - 2*u/(len(a)*len(b))


def main():
    doover = os.environ.get("DOOVER", ".")
    refdir = os.environ.get("REFDIR", os.path.join(doover, "data"))
    d1 = os.path.join(doover, "output", "01_Description")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dcr-sites", default=os.path.join(d1, "SPLICE_SITES_filtered_DCE_splice_sites.bed"))
    ap.add_argument("--ctrl-sites", default=os.path.join(d1, "SPLICE_SITES_control_splice_sites.bed"))
    ap.add_argument("--phylop", default=os.environ.get("PHYLOP_BW", os.path.join(refdir, "hg38.phyloP100way.bw")))
    ap.add_argument("--outdir", default=os.path.join(os.environ.get("WORK", "work"), "controls"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    ss = pd.concat([load(a.dcr_sites, "DCR"), load(a.ctrl_sites, "control")], ignore_index=True)
    ss = ss.sort_values(["chrom", "start"]).reset_index(drop=True)

    bw = pyBigWig.open(a.phylop); chroms = set(bw.chroms())
    t0 = time.time(); vals = []
    for i, r in enumerate(ss.itertuples()):
        if r.chrom not in chroms: vals.append(np.nan); continue
        try:
            v = bw.stats(r.chrom, int(r.start), int(r.end), type="mean")[0]
            vals.append(np.nan if v is None else v)
        except Exception:
            vals.append(np.nan)
        if i and i % 20000 == 0: print(f"{i}/{len(ss)} {time.time()-t0:.0f}s", flush=True)
    bw.close()
    ss["phyloP"] = vals
    ss.to_csv(os.path.join(a.outdir, "splice_site_phylop.tsv"), sep="\t", index=False)
    print(f"done {len(ss)} in {time.time()-t0:.0f}s | scored {ss.phyloP.notna().sum()}")
    print(ss.groupby("cls").phyloP.agg(["count", "mean", "median"]).round(3).to_string())

    d = ss.dropna(subset=["phyloP"])
    rows = []
    for key, g in list(d.groupby(["junction", "site"])) + [(("ALL", "ALL"), d)]:
        A = g.loc[g.cls.eq("DCR"), "phyloP"].values
        B = g.loc[g.cls.eq("control"), "phyloP"].values
        if len(A) < 10 or len(B) < 10:
            continue
        p = stats.mannwhitneyu(A, B, alternative="two-sided").pvalue
        rows.append(dict(junction=key[0], site=key[1], n_DCR=len(A), n_ctrl=len(B),
                         med_DCR=np.median(A), med_ctrl=np.median(B),
                         mean_DCR=A.mean(), mean_ctrl=B.mean(),
                         diff_med=np.median(A)-np.median(B), rank_biserial=rb(A, B), p=p))
    res = pd.DataFrame(rows)
    res["p_adj"] = np.minimum(res.p * len(res), 1.0)
    res.to_csv(os.path.join(a.outdir, "splice_site_conservation.csv"), index=False)
    print(res.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
