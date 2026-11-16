#!/usr/bin/env python3
"""
A2 step 3 - NMD (pipeline 50-nt rule, identical code for all arms), conservation (phyloP100way, length/GC
matched) and descriptive length / frame-preservation comparisons between DCRs and regular AS.
Inputs: out/A2_events.tsv (A2_event_tables.py), out/NMDtable_allgenes.tsv (nmd_allgenes.py = verbatim
        analyze_isoforms() from DOOVER isoform_NMD.py run on all genes), in/isoforms_to_DCR_metadata.tsv,
        in/DCR_NMD_analysis.out (official pipeline DCR-level NMD result), out/isoform_frames.tsv
Outputs: out/A2_events_annotated.csv, out/A2_table_*.csv
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
_args = _ap.parse_args()
IN = lambda f: _os.path.join(_args.indir, f)
OUT = lambda f: _os.path.join(_args.outdir, f)
_os.makedirs(_args.outdir, exist_ok=True)
import pandas as pd, numpy as np, statsmodels.api as sm, statsmodels.formula.api as smf
from scipy import stats
from statsmodels.discrete.conditional_models import ConditionalLogit
rng = np.random.default_rng(1)

E = pd.read_csv(OUT("A2_events.tsv"), sep="\t", low_memory=False)
dcr = pd.read_csv(IN("filtered_DCR_overlapping_regions_with_gene_and_dceID.bed"), sep="\t", header=None,
                  names=["chr", "start", "end", "DCR", "strand", "gene", "frames"])
N = pd.read_csv(OUT("NMDtable_allgenes.tsv"), sep="\t", header=None, names=["gene", "isoform", "strand", "dist", "pred"])
O = pd.read_csv(IN("NMDtable.tsv"), sep="\t", header=None, names=["gene", "isoform", "strand", "dist", "pred"])
qc = pd.read_csv(OUT("isoform_frames.tsv"), sep="\t")
md = pd.read_csv(IN("isoforms_to_DCR_metadata.tsv"), sep="\t", header=None, names=["isoform", "frame", "termini", "DCR"])
md = md[md.DCR.isin(set(dcr.DCR))]
off = pd.read_csv(IN("DCR_NMD_analysis.out"), sep="\t"); off.columns = ["DCR", "outcome", "reason", "frames"]
off = off[off.DCR.isin(set(dcr.DCR))]

# ---------- validation of the NMD rerun against the pipeline's own table
v = O.merge(N, on="isoform", suffixes=("_o", "_n"))
VAL = dict(orig_isoforms=len(O), rerun_found=len(v), pred_identical=int((v.pred_o == v.pred_n).sum()),
           dist_identical=int((v.dist_o.astype(str) == v.dist_n.astype(str)).sum()),
           rerun_isoforms_total=len(N), rerun_NMD=int((N.pred == "NMD").sum()),
           official_DCR_in_NMD_zone=int((off.reason == "in_NMD_zone").sum()), official_DCR_n=off.DCR.nunique())

# genes whose merged model is a single exon: the verbatim code re-uses EJCindex of the previous gene -> order
# dependent call; excluded from NMD comparisons
iso_nmd = pd.read_csv(IN("UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed"), sep="\t", header=None,
                      names=["chr", "start", "end", "isoform", "strand", "gene"])
def nmerged(t):
    n, e = 0, -1
    for a, b in zip(t.start, t.end):
        if a > e: n += 1
        e = max(e, b)
    return n
nm = iso_nmd.sort_values("start").groupby("gene").apply(nmerged)
single = set(nm[nm == 1].index)
pred = dict(zip(N.isoform, N.pred))
niso = qc.groupby("gene").size()

def arm_rule(arms):
    """pipeline core rule (analyze_DCR): one arm 100% NMD while another arm is not."""
    r = []
    for a in arms:
        p = [pred[i] for i in a if i in pred]
        if not p: return (False, False, None)
        r.append(np.mean([x == "NMD" for x in p]))
    allN = [k for k, x in enumerate(r) if x == 1]
    return (True, bool(allN) and any(x < 1 for x in r), allN)

rec = {}
for D_, t in md.groupby("DCR"):
    arms = [list(s.isoform) for _, s in t.groupby("frame")]
    ev, y, _ = arm_rule(arms)
    rec["DCR" + str(D_)] = dict(nmd_evaluable=ev, nmd_arm_specific=y, nmd_arm=None,
                                n_carrier_iso=t.isoform.nunique(), min_arm=min(len(a) for a in arms))
for r in E[E.arm.str.startswith("regAS")].itertuples():
    inc = r.cov_isoforms.split(","); exc = r.spanning_excluders.split(",")
    ev, y, allN = arm_rule([inc, exc])
    rec[r.event_id] = dict(nmd_evaluable=ev, nmd_arm_specific=y,
                           nmd_arm=(None if not y else ("inclusion" if 0 in allN else "exclusion")),
                           n_carrier_iso=len(inc), min_arm=min(len(inc), len(exc)))
E = E.merge(pd.DataFrame(rec).T.rename_axis("event_id").reset_index(), on="event_id", how="left")
for c_ in ["n_carrier_iso", "min_arm"]: E[c_] = pd.to_numeric(E[c_])
for c_ in ["nmd_evaluable", "nmd_arm_specific"]: E[c_] = E[c_].astype("boolean")
E["n_iso_gene"] = E.gene.map(niso)
E["incl_frac"] = E.n_carrier_iso / E.n_iso_gene
E["single_model_exon_gene"] = E.gene.isin(single)
E = E.merge(off[["DCR", "reason"]].rename(columns={"reason": "official_NMD_reason"}), on="DCR", how="left")
E["arm2"] = np.where(E.arm.str.startswith("regAS"),
                     np.where(E.mult3, E.arm.str.replace("regAS", "regASfp"), E.arm.str.replace("regAS", "fsAS")), E.arm)
E["loglen"] = np.log10(E.length); E["logiso"] = np.log2(E.n_iso_gene); E["logminarm"] = np.log2(E.min_arm.astype(float))
E.to_csv(OUT("A2_events_annotated.csv"), index=False)

ORDER = ["DCR", "regASfp_DCRgenes", "regASfp_nonDCRgenes", "fsAS_DCRgenes", "fsAS_nonDCRgenes",
         "const_DCRgenes", "const_nonDCRgenes"]
def cliffs(x, y):
    u = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 2 * u.statistic / (len(x) * len(y)) - 1, u.pvalue

# ---------- (c) descriptive
desc = E.groupby("arm2").agg(n_events=("event_id", "size"), n_genes=("gene", "nunique"),
                             length_median=("length", "median"), length_q1=("length", lambda x: x.quantile(.25)),
                             length_q3=("length", lambda x: x.quantile(.75)), frac_mult3=("mult3", "mean"),
                             GC_median=("GC", "median"), isoforms_per_gene_median=("n_iso_gene", "median"),
                             incl_frac_median=("incl_frac", "median")).reindex(ORDER)
lenrows = []
for ref in ORDER[1:5]:
    a = E.loc[E.arm2 == "DCR", "length"]; b = E.loc[E.arm2 == ref, "length"]
    dl, pl = cliffs(a, b)
    m3a = E.loc[E.arm2 == "DCR", "mult3"]; m3 = E.loc[E.arm.isin([ref.replace("regASfp", "regAS").replace("fsAS", "regAS")]), "mult3"]
    lenrows.append(dict(reference=ref, DCR_len_median=a.median(), ref_len_median=b.median(), cliffs_delta_len=dl, p_len=pl))
desc.to_csv(OUT("A2_table_descriptive.csv"))
pd.DataFrame(lenrows).to_csv(OUT("A2_table_length_tests.csv"), index=False)
# frame preservation: DCR vs all internal regular AS (before the frame-preservation split)
fp = []
for ref in ["regAS_DCRgenes", "regAS_nonDCRgenes"]:
    a = E.loc[E.arm == "DCR", "mult3"]; b = E.loc[E.arm == ref, "mult3"]
    ct = np.array([[a.sum(), (~a).sum()], [b.sum(), (~b).sum()]])
    orr, p = stats.fisher_exact(ct)
    fp.append(dict(reference=ref, DCR_frac_mult3=a.mean(), ref_frac_mult3=b.mean(), n_DCR=len(a), n_ref=len(b), OR=orr, p=p))
pd.DataFrame(fp).to_csv(OUT("A2_table_frame_preservation.csv"), index=False)

# ---------- (a) NMD
T = E[E.arm2.isin(ORDER[:5]) & ~E.single_model_exon_gene & E.nmd_evaluable.fillna(False).astype(bool)].copy()
T["y"] = T.nmd_arm_specific.astype(bool).astype(int)
nmdrows = []
for a2 in ORDER[:5]:
    t = T[T.arm2 == a2]
    nmdrows.append(dict(arm=a2, n=len(t), genes=t.gene.nunique(), nmd=int(t.y.sum()), frac=t.y.mean(),
                        inclusion_arm=int((t.nmd_arm == "inclusion").sum()), exclusion_arm=int((t.nmd_arm == "exclusion").sum())))
nmd_desc = pd.DataFrame(nmdrows)
tests = []
for ref in ORDER[1:5]:
    d = T[T.arm2.isin(["DCR", ref])].copy(); d["isDCR"] = (d.arm2 == "DCR").astype(int)
    p1, p0 = d[d.isDCR == 1].y.mean(), d[d.isDCR == 0].y.mean()
    ct = np.array([[d[d.isDCR == 1].y.sum(), (1 - d[d.isDCR == 1].y).sum()], [d[d.isDCR == 0].y.sum(), (1 - d[d.isDCR == 0].y).sum()]])
    orr, pf = stats.fisher_exact(ct)
    g = smf.gee("y ~ isDCR + loglen + logiso + logminarm", groups="gene", data=d, family=sm.families.Binomial(),
                cov_struct=sm.cov_struct.Exchangeable()).fit()
    ci = g.conf_int().loc["isDCR"]
    row = dict(reference=ref, DCR_frac=p1, ref_frac=p0, risk_diff=p1 - p0, cohens_h=2 * np.arcsin(np.sqrt(p1)) - 2 * np.arcsin(np.sqrt(p0)),
               OR_raw=orr, p_fisher=pf, OR_adj_GEE=np.exp(g.params.isDCR), OR_adj_lo=np.exp(ci[0]), OR_adj_hi=np.exp(ci[1]), p_adj_GEE=g.pvalues.isDCR)
    if ref.endswith("_DCRgenes"):
        cl = ConditionalLogit(d.y, d[["isDCR", "loglen", "logminarm"]], groups=d.gene).fit(disp=0)
        cci = cl.conf_int().loc["isDCR"]
        row.update(OR_within_gene=np.exp(cl.params.isDCR), OR_wg_lo=np.exp(cci[0]), OR_wg_hi=np.exp(cci[1]), p_within_gene=cl.pvalues.isDCR,
                   n_informative_genes=int(d.groupby("gene").y.nunique().gt(1).sum()))
    tests.append(row)
nmd_desc.to_csv(OUT("A2_table_NMD_counts.csv"), index=False)
pd.DataFrame(tests).to_csv(OUT("A2_table_NMD_tests.csv"), index=False)

# ---------- (b) conservation, length + GC matched
C = E[E.phyloP_bp_frac >= 0.9].copy()
def nn_match(case, pool, k=5, cal_len=0.1, cal_gc=0.03):
    X = pool[["loglen", "GC"]].to_numpy(); out = []
    for r in case.itertuples():
        dl = np.abs(X[:, 0] - r.loglen); dg = np.abs(X[:, 1] - r.GC)
        ok = np.where((dl <= cal_len) & (dg <= cal_gc))[0]
        if len(ok) == 0: continue
        dist = (dl[ok] / cal_len) ** 2 + (dg[ok] / cal_gc) ** 2
        sel = ok[np.argsort(dist + rng.random(len(ok)) * 1e-9)[:k]]
        out.append((r.event_id, r.phyloP_mean, pool.phyloP_mean.to_numpy()[sel].mean(), len(sel),
                    pool.loglen.to_numpy()[sel].mean(), pool.GC.to_numpy()[sel].mean()))
    return pd.DataFrame(out, columns=["event_id", "pp_DCR", "pp_match", "k", "loglen_match", "GC_match"])
cons = []; matched_all = []
case = C[C.arm2 == "DCR"]
for ref in ORDER[1:]:
    pool = C[C.arm2 == ref]
    if len(pool) > 20000: pool = pool.sample(20000, random_state=1)   # constitutive context arm only
    dr, pr = cliffs(case.phyloP_mean, C[C.arm2 == ref].phyloP_mean)
    M = nn_match(case, pool); M["reference"] = ref; matched_all.append(M)
    dd = M.pp_DCR - M.pp_match
    w = stats.wilcoxon(dd)
    rpos = stats.rankdata(np.abs(dd)); rb = (rpos[dd > 0].sum() - rpos[dd < 0].sum()) / rpos.sum()
    cc = case.set_index("event_id").loc[M.event_id]
    smd_len = (cc.loglen.mean() - M.loglen_match.mean()) / case.loglen.std()
    smd_gc = (cc.GC.mean() - M.GC_match.mean()) / case.GC.std()
    d = C[C.arm2.isin(["DCR", ref])].copy(); d["isDCR"] = (d.arm2 == "DCR").astype(int)
    ols = smf.ols("phyloP_mean ~ isDCR + loglen + GC + I(GC**2) + logiso", data=d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.gene)[0]})
    row = dict(reference=ref, n_DCR=len(case), n_ref=int((C.arm2 == ref).sum()),
               DCR_median=case.phyloP_mean.median(), ref_median=C[C.arm2 == ref].phyloP_mean.median(),
               cliffs_delta_raw=dr, p_raw=pr, n_DCR_matched=len(M), median_k=M.k.median(),
               SMD_loglen_after=smd_len, SMD_GC_after=smd_gc,
               matched_median_diff=dd.median(), matched_mean_diff=dd.mean(), matched_rank_biserial=rb, p_matched_wilcoxon=w.pvalue,
               OLS_beta_DCR=ols.params.isDCR, OLS_lo=ols.conf_int().loc["isDCR", 0], OLS_hi=ols.conf_int().loc["isDCR", 1], OLS_p=ols.pvalues.isDCR)
    if ref.endswith("_DCRgenes") and not ref.startswith("const"):
        both = d.groupby(["gene", "isDCR"]).phyloP_mean.mean().unstack().dropna()
        gw = stats.wilcoxon(both[1] - both[0])
        row.update(within_gene_n=len(both), within_gene_median_diff=(both[1] - both[0]).median(), within_gene_p=gw.pvalue)
        # sensitivity: additionally adjust for inclusion fraction
        o2 = smf.ols("phyloP_mean ~ isDCR + loglen + GC + I(GC**2) + logiso + incl_frac", data=d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.gene)[0]})
        row.update(OLS_inclfrac_beta_DCR=o2.params.isDCR, OLS_inclfrac_p=o2.pvalues.isDCR)
    cons.append(row)
pd.DataFrame(cons).to_csv(OUT("A2_table_conservation.csv"), index=False)

# ---------- (b') inclusion-level control. 66% of DCRs are carried by every isoform of their gene (incl_frac = 1),
# a state that regular AS cannot have by definition. Compare like with like:
#   DCRs with incl_frac < 1  vs regular AS matched on length, GC AND inclusion fraction (caliper 0.1)
#   DCRs with incl_frac = 1  vs constitutive segments (inclusion 1) matched on length and GC
def nn_match3(case, pool, k=5, cal=(0.1, 0.03, 0.1)):
    X = pool[["loglen", "GC", "incl_frac"]].to_numpy(); out = []
    for r in case.itertuples():
        dv = np.abs(X - np.array([r.loglen, r.GC, r.incl_frac]))
        ok = np.where((dv <= np.array(cal)).all(1))[0]
        if len(ok) == 0: continue
        dist = ((dv[ok] / np.array(cal)) ** 2).sum(1)
        sel = ok[np.argsort(dist)[:k]]
        out.append((r.event_id, r.phyloP_mean, pool.phyloP_mean.to_numpy()[sel].mean(), len(sel)))
    return pd.DataFrame(out, columns=["event_id", "pp_DCR", "pp_match", "k"])
inc = []
strata = C.copy(); strata["incl_stratum"] = pd.cut(strata.incl_frac, [0, 0.5, 0.75, 0.999, 1.01], right=False,
                                                   labels=["<0.5", "0.5-0.75", "0.75-<1", "1"])
strata[strata.arm2.isin(ORDER[:5])].groupby(["incl_stratum", "arm2"], observed=True).phyloP_mean.agg(["size", "median", "mean"]) \
    .reset_index().to_csv(OUT("A2_table_conservation_by_inclusion.csv"), index=False)
case_lt1 = C[(C.arm2 == "DCR") & (C.incl_frac < 0.999)]
case_eq1 = C[(C.arm2 == "DCR") & (C.incl_frac >= 0.999)]
for ref, cs, fn in [("regASfp_DCRgenes", case_lt1, nn_match3), ("regASfp_nonDCRgenes", case_lt1, nn_match3),
                    ("fsAS_DCRgenes", case_lt1, nn_match3), ("fsAS_nonDCRgenes", case_lt1, nn_match3),
                    ("const_DCRgenes", case_eq1, nn_match), ("const_nonDCRgenes", case_eq1, nn_match)]:
    pool = C[C.arm2 == ref]
    if len(pool) > 20000: pool = pool.sample(20000, random_state=1)
    M = fn(cs, pool); dd = M.pp_DCR - M.pp_match; w = stats.wilcoxon(dd)
    rpos = stats.rankdata(np.abs(dd)); rb = (rpos[dd > 0].sum() - rpos[dd < 0].sum()) / rpos.sum()
    inc.append(dict(DCR_subset="incl_frac<1" if cs is case_lt1 else "incl_frac=1", reference=ref, n_DCR_subset=len(cs),
                    n_matched=len(M), matched_median_diff=dd.median(), matched_mean_diff=dd.mean(),
                    matched_rank_biserial=rb, p_matched_wilcoxon=w.pvalue))
    M["reference"] = ref; M["DCR_subset"] = inc[-1]["DCR_subset"]; matched_all.append(M)
pd.DataFrame(inc).to_csv(OUT("A2_table_conservation_inclusion_matched.csv"), index=False)
# NMD with inclusion fraction as covariate (DCRs with incl_frac < 1 only, so the covariate overlaps)
nmd_inc = []
for ref in ORDER[1:5]:
    d = T[T.arm2.isin(["DCR", ref]) & (T.incl_frac < 0.999)].copy(); d["isDCR"] = (d.arm2 == "DCR").astype(int)
    g = smf.gee("y ~ isDCR + loglen + logiso + logminarm + incl_frac", groups="gene", data=d, family=sm.families.Binomial(),
                cov_struct=sm.cov_struct.Exchangeable()).fit()
    ci = g.conf_int().loc["isDCR"]
    nmd_inc.append(dict(reference=ref, n_DCR=int(d.isDCR.sum()), DCR_frac=d[d.isDCR == 1].y.mean(), ref_frac=d[d.isDCR == 0].y.mean(),
                        OR_adj=np.exp(g.params.isDCR), lo=np.exp(ci[0]), hi=np.exp(ci[1]), p=g.pvalues.isDCR))
pd.DataFrame(nmd_inc).to_csv(OUT("A2_table_NMD_inclusion_adjusted.csv"), index=False)
T.groupby([pd.cut(T.incl_frac, [0, 0.5, 0.75, 0.999, 1.01], right=False, labels=["<0.5", "0.5-0.75", "0.75-<1", "1"]), "arm2"],
          observed=True).y.agg(["size", "sum", "mean"]).reset_index().to_csv(OUT("A2_table_NMD_by_inclusion.csv"), index=False)
print(pd.DataFrame(inc).round(4).to_string()); print(pd.DataFrame(nmd_inc).round(4).to_string())
pd.concat(matched_all).to_csv(OUT("A2_matched_pairs.csv"), index=False)
pd.Series(VAL).to_csv(OUT("A2_table_NMD_validation.csv"), header=["value"])
print(pd.Series(VAL).to_string()); print(desc.round(3).to_string()); print(pd.DataFrame(fp).round(4).to_string())
print(nmd_desc.round(3).to_string()); print(pd.DataFrame(tests).round(4).T.to_string()); print(pd.DataFrame(cons).round(4).T.to_string())
