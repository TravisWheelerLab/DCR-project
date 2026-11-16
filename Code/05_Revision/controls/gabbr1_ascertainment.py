#!/usr/bin/env python3
"""gabbr1_ascertainment.py -- why showcase DCRs (GABBR1, PTBP2, ...) are or are not in the selection test.

Purpose  Joins the per-DCR-interval exclusion fates of the alternate-frame selection test (v1/v2 table
         dcr_exclusion_fates.csv: kept / too_short / lift_failed / indel / ...) to the interval
         coordinates of dcr_regions.fa, flags non-primary contigs (_alt/_random/_fix/chrUn), and
         reports the candidate showcase genes with their DCR intervals, frame count, fate, length and GC.
         Used to choose a showcase example other than GABBR1 (26-nt DCR, below the 30-nt length
         filter) and to document the ascertainment in docs/A3_results.md.
Inputs   --dcr-regions   dcr_regions.fa           (default $WORK/selection_inputs/dcr_regions.fa)
         --fates         dcr_exclusion_fates.csv  (default $WORK/selection_inputs/dcr_exclusion_fates.csv)
         --dcr-bed       $DOOVER/output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed
         --genes         comma-separated lower-case gene symbols (default: the 10 candidates examined)
Outputs  --outdir/selection_interval_fates_with_coords.csv  (idx, chrom, start, end, strand, fate, len_h, gc, nonprimary)
         --outdir/A3_showcase_candidates.csv
Manuscript  Supplementary Figure S6 (GABBR1) legend and Results text on the GABBR1 ascertainment;
            docs/A3_results.md.
"""
import argparse, os, re
import pandas as pd

CAND = "gabbr1,ptbp2,ccser2,plcb4,mark4,magi1,pbx1,atp2b1,atp2b4,madd"


def main():
    doover = os.environ.get("DOOVER", ".")
    work = os.environ.get("WORK", "work")
    si = os.path.join(work, "selection_inputs")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dcr-regions", default=os.path.join(si, "dcr_regions.fa"))
    ap.add_argument("--fates", default=os.path.join(si, "dcr_exclusion_fates.csv"))
    ap.add_argument("--dcr-bed", default=os.path.join(doover, "output", "01_Description",
                                                      "filtered_DCR_overlapping_regions_with_gene_and_dceID.bed"))
    ap.add_argument("--genes", default=CAND)
    ap.add_argument("--outdir", default=os.path.join(work, "controls"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    hdrs = [l[1:].strip() for l in open(a.dcr_regions) if l.startswith(">")]
    rec = []
    for h in hdrs:
        pre, coord = h.split("::", 1)
        idx = int(pre.split("|")[-1])
        m = re.match(r'(\S+?):(\d+)-(\d+)\(([+-])\)', coord)
        rec.append(dict(idx=idx, chrom=m.group(1), start=int(m.group(2)), end=int(m.group(3)), strand=m.group(4)))
    reg = pd.DataFrame(rec)

    fates = pd.read_csv(a.fates)
    fates["idx"] = fates.name.str.split("|").str[-1].astype(int)
    j = reg.merge(fates[["idx", "fate", "len_h", "gc"]], on="idx", how="left")
    j["nonprimary"] = j.chrom.str.contains("_alt|_random|_fix|^chrUn", regex=True)
    j.to_csv(os.path.join(a.outdir, "selection_interval_fates_with_coords.csv"), index=False)

    bed = pd.read_csv(a.dcr_bed, sep="\t", header=None,
                      names=["chrom", "start", "end", "dcrid", "strand", "gene", "nfr"])
    bed["len"] = bed.end - bed.start
    cand = a.genes.split(",")
    c = bed[bed.gene.isin(cand)].merge(j[["chrom", "start", "end", "strand", "fate", "len_h", "gc", "nonprimary"]],
                                       on=["chrom", "start", "end"], how="left", suffixes=("", "_r"))
    c = c[["gene", "chrom", "start", "end", "len", "nfr", "fate", "len_h", "gc", "nonprimary"]].sort_values(
        ["fate", "len_h"], ascending=[True, False])
    c.to_csv(os.path.join(a.outdir, "A3_showcase_candidates.csv"), index=False)
    print(c.to_string(index=False))


if __name__ == "__main__":
    main()
