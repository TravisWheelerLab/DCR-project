#!/usr/bin/env python3
"""gtex_fetch_medians.py -- per-gene GTEx v8 median expression summaries for DCR and background genes.

Purpose  Standalone replacement for the in-platform 'expression' connector used during the revision.
         (1) gene universe: DCR genes = gene column of the final DCR catalogue; background = every gene
             symbol of the SwissProt isoform models mapped by MIRAGE2 (UniprotKB_isoforms_mapping.bed);
         (2) symbol -> GENCODE v26 id via the public GTEx Portal API v2 (/reference/gene), keeping one
             'protein coding' id per symbol (--tie-break, see below);
         (3) median TPM per tissue (/expression/medianGeneExpression, datasetId=gtex_v8, 54 tissue
             sites), summarised per gene: n_tissues, mean, median (median of tissue medians), max,
             tau = sum(1 - x_i/max)/(n-1), top_tissue.
Inputs   --dcr-bed  $DOOVER/output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed
         --iso-bed  $DOOVER/output/01_Description/UniprotKB_isoforms_mapping.bed
         network access to https://gtexportal.org/api/v2 (no key required)
Outputs  --out  gtex_v8_gene_medians.json  {"genes": {gencodeId: {...}}, "n_failed": int,
                                            "n_symbols": int, "n_resolved": int, ...}
         --genes-out  a4_genes.json         {"dcr": [...], "background": [...]}
Runtime  ~3 min (about 90 + 170 HTTP requests); re-runs reuse --cache-dir responses.
Manuscript  feeds controls/expression_tau.py -> Supplementary Figure S9 B-C, Results (expression class
            and tissue specificity; brain dominance by mouse-conservation stratum).
Note     The values used in the manuscript were fetched through the connector (same GTEx v8 /
         GENCODE v26 endpoints). This script reproduces that procedure; compare with --check-against.
"""
import argparse, hashlib, json, os, statistics as st, sys, time, urllib.error, urllib.parse, urllib.request

API = "https://gtexportal.org/api/v2"


def get_json(path, params, cache_dir=None, retries=5):
    q = urllib.parse.urlencode(params, doseq=True)
    url = f"{API}{path}?{q}"
    cf = None
    if cache_dir:
        cf = os.path.join(cache_dir, hashlib.sha1(url.encode()).hexdigest() + ".json")
        if os.path.exists(cf):
            return json.load(open(cf))
    for k in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}),
                                        timeout=120) as r:
                d = json.load(r)
            if cf:
                json.dump(d, open(cf, "w"))
            return d
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            if k == retries - 1:
                raise
            time.sleep(2 ** k)


def get_all_pages(path, params, cache_dir=None, per_page=100000):
    out, page = [], 0
    while True:
        d = get_json(path, dict(params, page=page, itemsPerPage=per_page), cache_dir)
        out.extend(d.get("data", []))
        pi = d.get("paging_info", {})
        page += 1
        if page >= pi.get("numberOfPages", 1):
            return out


def gene_universe(dcr_bed, iso_bed):
    dcr = set()
    for line in open(dcr_bed):
        f = line.rstrip("\n").split("\t")
        if len(f) > 5:
            dcr.add(f[5].lower())
    bg = set()
    with open(iso_bed) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 3 and "|" in f[3]:
                part = f[3].split("|")
                if len(part) > 1 and "/" in part[1]:
                    bg.add(part[1].split("/")[0].strip().lower())
    return sorted(dcr), sorted(x for x in bg if x)


def main():
    doover = os.environ.get("DOOVER", ".")
    d1 = os.path.join(doover, "output", "01_Description")
    work = os.environ.get("WORK", "work")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dcr-bed", default=os.path.join(d1, "filtered_DCR_overlapping_regions_with_gene_and_dceID.bed"))
    ap.add_argument("--iso-bed", default=os.path.join(d1, "UniprotKB_isoforms_mapping.bed"))
    ap.add_argument("--out", default=os.path.join(work, "controls", "gtex_v8_gene_medians.json"))
    ap.add_argument("--genes-out", default=os.path.join(work, "controls", "a4_genes.json"))
    ap.add_argument("--cache-dir", default=os.path.join(work, "controls", "gtex_cache"))
    ap.add_argument("--resolve-batch", type=int, default=200)
    ap.add_argument("--median-batch", type=int, default=100)
    ap.add_argument("--dataset", default="gtex_v8")
    ap.add_argument("--gencode-version", default="v26")
    ap.add_argument("--tie-break", choices=["lowest", "first"], default="lowest",
                    help="symbol with several protein-coding ids: lowest ENSG accession (default) or first API hit")
    ap.add_argument("--limit", type=int, default=0, help="first N background symbols only (testing)")
    ap.add_argument("--check-against", default=None, help="previous JSON to compare per-gene values with")
    a = ap.parse_args()
    for p in (os.path.dirname(a.out), os.path.dirname(a.genes_out), a.cache_dir):
        if p:
            os.makedirs(p, exist_ok=True)

    dcr_genes, bg = gene_universe(a.dcr_bed, a.iso_bed)
    json.dump({"dcr": dcr_genes, "background": bg}, open(a.genes_out, "w"))
    dcr = set(dcr_genes)
    syms = bg[: a.limit] if a.limit else bg
    print(f"background {len(bg)} | DCR {len(dcr)} | DCR within background {len(dcr & set(bg))}", flush=True)

    # 1) symbols -> GENCODE ids (first protein-coding hit per symbol, as in the connector run)
    sym2id, hits, t0, rfail = {}, {}, time.time(), 0
    for i in range(0, len(syms), a.resolve_batch):
        chunk = [s.upper() for s in syms[i:i + a.resolve_batch]]
        try:
            rows = get_all_pages("/reference/gene", dict(geneId=chunk, gencodeVersion=a.gencode_version,
                                                         genomeBuild="GRCh38/hg38"), a.cache_dir, per_page=1000)
        except Exception as e:
            rfail += len(chunk); print("resolve fail", i, type(e).__name__, str(e)[:100], flush=True); continue
        for x in rows:
            s = (x.get("geneSymbol") or "").lower()
            if s and x.get("geneType") == "protein coding":
                hits.setdefault(s, []).append(x["gencodeId"])
    for s, ids_ in hits.items():
        # ~0.2% of symbols (28 of 17,472 here) have >1 protein-coding GENCODE v26 id; 'first' keeps API order (what the
        # connector run did, with the connector's order), 'lowest' keeps the lowest ENSG accession.
        sym2id[s] = ids_[0] if a.tie_break == "first" else min(ids_, key=lambda g: g.split(".")[0])
    n_multi = sum(len(set(v)) > 1 for v in hits.values())
    print(f"resolved {len(sym2id)}/{len(syms)} in {time.time()-t0:.0f}s ({n_multi} symbols with >1 "
          f"protein-coding id, tie-break={a.tie_break})", flush=True)

    # 2) median TPM -> per-gene summaries
    id2sym = {v: k for k, v in sym2id.items()}
    ids = list(id2sym)
    out, fails, t0 = {}, [], time.time()
    for bi, i in enumerate(range(0, len(ids), a.median_batch)):
        batch = ids[i:i + a.median_batch]
        try:
            rows = get_all_pages("/expression/medianGeneExpression", dict(gencodeId=batch, datasetId=a.dataset),
                                 a.cache_dir)
        except Exception as e:
            fails.extend(batch); print("median fail", i, type(e).__name__, str(e)[:100], flush=True); continue
        per = {}
        for row in rows:
            per.setdefault(row["gencodeId"], []).append((row["tissueSiteDetailId"], row["median"]))
        for gid, vals in per.items():
            v = [x[1] for x in vals if x[1] is not None]
            if not v: continue
            mx = max(v)
            tau = sum(1 - x/mx for x in v)/(len(v)-1) if mx > 0 and len(v) > 1 else None
            top = max(vals, key=lambda x: (x[1] is not None, x[1]))[0]
            out[gid] = dict(symbol=id2sym[gid], n_tissues=len(v), mean=st.mean(v),
                            median=st.median(v), max=mx, tau=tau, top_tissue=top,
                            is_dcr=id2sym[gid] in dcr)
        if bi % 20 == 0:
            print(f"batch {bi} genes={len(out)} {time.time()-t0:.0f}s", flush=True)

    json.dump({"genes": out, "n_failed": len(fails), "n_resolve_failed": rfail, "n_symbols": len(syms),
               "n_resolved": len(sym2id), "dataset": a.dataset, "gencode_version": a.gencode_version,
               "source": API}, open(a.out, "w"))
    nd = sum(1 for v in out.values() if v["is_dcr"])
    print(f"DONE genes={len(out)} dcr={nd} nondcr={len(out)-nd} failed={len(fails)}")

    if a.check_against:
        old = json.load(open(a.check_against))["genes"]
        common = set(old) & set(out)
        diff = [g for g in common if any(
            (old[g][k] is None) != (out[g][k] is None) or
            (old[g][k] is not None and abs(float(old[g][k]) - float(out[g][k])) > 1e-6)
            for k in ("median", "max", "tau"))]
        top = sum(old[g]["top_tissue"] != out[g]["top_tissue"] for g in common)
        print(f"check: old {len(old)} new {len(out)} common {len(common)} | only_old {len(set(old)-set(out))} "
              f"only_new {len(set(out)-set(old))} | numeric mismatches {len(diff)} | top_tissue mismatches {top}")


if __name__ == "__main__":
    sys.exit(main())
