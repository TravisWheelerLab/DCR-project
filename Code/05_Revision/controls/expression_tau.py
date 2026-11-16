#!/usr/bin/env python3
"""expression_tau.py -- DCR incidence by expression class, tissue specificity (tau) and brain dominance
by mouse-conservation stratum (reviewer items R2-M5 / R2-M6).

Purpose  From per-gene GTEx v8 summaries (controls/gtex_fetch_medians.py):
         (a) % DCR genes per quintile of median TPM (+ binomial 95% CI) and a logistic regression
             is_dcr ~ log10(median TPM + 1);
         (b) tau of DCR vs non-DCR genes within expression quintiles (Mann-Whitney U, rank-biserial);
         (c) per-gene table (symbol, median TPM, tau, top tissue, brain-top flag, quintile);
         (d) % genes whose top tissue is a brain region: non-DCR genes vs DCR genes split by mouse
             conservation of their DCRs (DCRmiceStatus.tsv; Fisher exact tests printed).
Inputs   --expr       gtex_v8_gene_medians.json (default $WORK/controls/gtex_v8_gene_medians.json)
         --mice-status $DOOVER/output/02_Conservation/DCRmiceStatus.tsv
         --dce2gene    $DOOVER/output/02_Conservation/DCEtoGene.txt
Outputs  --outdir/A4_dcr_by_expression.csv, A4_tau_by_expression.csv, A4_gene_expression_tau.csv,
         A4_brain_by_conservation.csv, A4_tests.txt (logistic + Fisher results printed by the session)
Manuscript  Supplementary Figure S9 B-C (figures/supp_splice_expression_tissue.R); Results (expression
            class / tissue specificity); docs/A4_results.md.
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm

QLAB = ["Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"]


def main():
    doover = os.environ.get("DOOVER", ".")
    oc = os.path.join(doover, "output", "02_Conservation")
    work = os.environ.get("WORK", "work")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--expr", default=os.path.join(work, "controls", "gtex_v8_gene_medians.json"))
    ap.add_argument("--mice-status", default=os.path.join(oc, "DCRmiceStatus.tsv"))
    ap.add_argument("--dce2gene", default=os.path.join(oc, "DCEtoGene.txt"))
    ap.add_argument("--outdir", default=os.path.join(work, "controls"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    out = lambda f: os.path.join(a.outdir, f)
    log = open(out("A4_tests.txt"), "w")
    def say(s):
        print(s); log.write(s + "\n")

    E = json.load(open(a.expr))
    ex = pd.DataFrame(E["genes"]).T.reset_index(names="gencodeId")
    for c in ["n_tissues", "mean", "median", "max", "tau"]: ex[c] = pd.to_numeric(ex[c], errors="coerce")
    ex["is_dcr"] = ex.is_dcr.astype(bool)
    ex["logexp"] = np.log10(ex["median"] + 1)
    say(f"genes={len(ex)} dcr={ex.is_dcr.sum()} | tau non-null={ex.tau.notna().sum()} | failed batches={E['n_failed']}")

    # (a) R2-M5: DCR incidence across expression classes (quintiles of median TPM)
    ex["expr_bin"] = pd.qcut(ex["median"], 5, labels=QLAB)
    g = ex.groupby("expr_bin", observed=True).agg(n=("is_dcr", "size"), n_dcr=("is_dcr", "sum"),
                                                  med_TPM=("median", "median"))
    g["pct_dcr"] = 100*g.n_dcr/g.n
    se = np.sqrt(g.pct_dcr/100*(1-g.pct_dcr/100)/g.n)*100
    g["lo"], g["hi"] = g.pct_dcr-1.96*se, g.pct_dcr+1.96*se
    say("=== R2-M5: DCR incidence by expression quintile ===\n" + g.round(3).to_string())
    lr = sm.Logit(ex.is_dcr.astype(int), sm.add_constant(ex[["logexp"]])).fit(disp=0)
    say(f"logistic is_dcr ~ log10(medTPM+1): beta={lr.params.logexp:.3f} OR={np.exp(lr.params.logexp):.3f} "
        f"p={lr.pvalues.logexp:.3g}")
    g.to_csv(out("A4_dcr_by_expression.csv"))

    # (b) R2-M6: tau within expression quintiles
    t = ex.dropna(subset=["tau"]).copy()
    rows = []
    for b, gg in t.groupby("expr_bin", observed=True):
        a_, b_ = gg.loc[gg.is_dcr, "tau"], gg.loc[~gg.is_dcr, "tau"]
        if len(a_) < 10: continue
        uu = stats.mannwhitneyu(a_, b_, alternative="two-sided")
        rows.append(dict(expr_bin=b, n_dcr=len(a_), n_non=len(b_), tau_dcr=a_.median(),
                         tau_non=b_.median(), diff=a_.median()-b_.median(),
                         rank_biserial=1-2*uu.statistic/(len(a_)*len(b_)), p=uu.pvalue))
    pd.DataFrame(rows).to_csv(out("A4_tau_by_expression.csv"), index=False)

    # (c) per-gene table
    t["brain_top"] = t.top_tissue.str.startswith("Brain")
    t[["gencodeId", "symbol", "is_dcr", "median", "tau", "top_tissue", "brain_top", "expr_bin"]].to_csv(
        out("A4_gene_expression_tau.csv"), index=False)

    # (d) brain dominance by mouse-conservation stratum
    rows = []
    with open(a.mice_status) as fh:
        next(fh)
        for line in fh:
            p_ = line.split()
            if len(p_) < 3: continue
            rows.append((int(p_[0]), " ".join(p_[1:-1]), p_[-1]))
    ms = pd.DataFrame(rows, columns=["DCE_ID", "mice_status", "group"])
    d2g = pd.read_csv(a.dce2gene, sep="\t", header=None, names=["DCE_ID", "gene"])
    mg = ms.merge(d2g, on="DCE_ID", how="left"); mg["gene"] = mg.gene.str.lower()
    cons = set(mg.loc[mg.mice_status.eq("DCR"), "gene"].dropna())
    noncons = set(mg.loc[mg.mice_status.eq("not DCR"), "gene"].dropna()) - cons
    tt = t.copy()
    tt["grp"] = np.where(~tt.is_dcr, "non-DCR gene",
                 np.where(tt.symbol.isin(cons), "DCR, mouse-conserved",
                 np.where(tt.symbol.isin(noncons), "DCR, not conserved", "DCR, unscored")))
    sm_ = tt.groupby("grp").agg(n=("brain_top", "size"), brain_pct=("brain_top", lambda s: 100*s.mean()),
                                tau_med=("tau", "median"))
    B_ = tt.loc[tt.grp.eq("non-DCR gene"), "brain_top"]
    for a_lab in ["DCR, mouse-conserved", "DCR, not conserved"]:
        A_ = tt.loc[tt.grp.eq(a_lab), "brain_top"]
        orr, pv = stats.fisher_exact([[A_.sum(), len(A_)-A_.sum()], [B_.sum(), len(B_)-B_.sum()]])
        say(f"{a_lab} (n={len(A_)}) vs non-DCR: OR={orr:.3f} p={pv:.4g}")
    sm_.to_csv(out("A4_brain_by_conservation.csv"))
    log.close()


if __name__ == "__main__":
    main()
