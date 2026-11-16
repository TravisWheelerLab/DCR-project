"""selection_v3.py -- OLGenie-style alternate-frame selection test for DCRs (v3).

Reimplementation of the lost v2 code (Goubert et al. DCR revision).
Manuscript: Figure 2B(iii) (alternate-frame selection), Methods "Alternate-frame selection test",
Results/Supplementary tables of the selection test. Methods detail: docs/selection_v3_methods.md and
docs/OLGenie_methods_explained.md.

Usage
  python selection_v3.py --inputs-dir $WORK/selection_inputs --outdir $WORK/selection [--boot 2000]
  (each input can be overridden individually, see --help). Runtime ~1 min on 1 core.

Inputs (prior-session intermediates; stage them with stage_selection_inputs.py, md5 in
selection_inputs_manifest.tsv):
  dcr_regions.fa, D_ucsc.fa, dcr_protein_frames.csv, dcr_exclusion_fates.csv, D_cds.ov, D_cds_any.ov,
  N_human.fa, N_ucsc.fa, N_cds.ov, cds_exons.tsv, nondcr_pool.fa
Outputs (--outdir):
  selection_omega_v3.csv              headline + sensitivity comparisons, both contig scopes  -> Fig 2B(iii)
  selection_by_framerel_v3.csv        ss13 (+1) / ss12 (+2) frame relations
  selection_length_matched_v3.csv     exact-length reweighting + k:1 matched bootstrap
  selection_by_lenbin_v3.csv          per length bin (not saved as an artifact in the original run)
  selection_v2_reproduction_check.csv reimplementation vs v2 published values
  selection_site_counts_v3.csv        per-interval site/difference counts and fates (~11 MB)

Design (Nelson, Ardern & Wei 2020 MBE, OLGenie):
  every nucleotide position whose codon is complete in BOTH frames contributes 1/3 site per
  alternative base to one of NN/NS/SN/SS (letter 1 = reference frame, letter 2 = alternate frame);
  observed human-mouse differences are classified the same way, single-site in the HUMAN codon
  context. Counts are pooled over intervals, p = diffs/sites, Jukes-Cantor corrected, then
    omega_alt = dNN/dNS ; omega_alt = dSN/dSS ; omega_ref (positive control) = dNN/dSN.

Arms
  D  (DCRs): human interval (dcr_regions.fa) vs UCSC-liftOver mouse (D_ucsc.fa); BOTH frames from
     the SwissProt isoform models (dcr_protein_frames.csv, offsets are 0-based in-interval offsets of
     the first full codon). Reference = lowest offset, alternate = second offset (v2 convention;
     arbitrary for DCRs - a swapped-role sensitivity is reported).
  N  (latent null): latent dual-coding intervals in non-DCR coding exons (N_human.fa) vs UCSC
     liftOver mouse (N_ucsc.fa). Reference frame = GENCODE v48 CDS frame: cumulative preceding CDS
     length -> codon position cpos at the interval's 5' base -> offset = (3 - cpos) % 3. Only
     intervals whose fully-containing CDS blocks give exactly ONE frame. Alternate frame = an open
     (stop-free) non-reference frame of the human interval.
     v2-reproducing rule (reverse-engineered, see methods): reference frame must itself be open;
     one open alternative -> use it; two open alternatives -> lowest offset when ref offset is 0/1,
     interval DROPPED when ref offset is 2.  Principled sensitivity rule ("orfpair"): the alternate
     frame is the frame of the ORF pair that defines the interval in its source exon.
Filters: length < 30 nt; failed lift; human/mouse length differ (indel).  v3 adds: non-primary
contigs (_alt, _random, _fix, chrUn*) removed from BOTH arms.
"""
import argparse, os, sys
import numpy as np, pandas as pd, ast, bisect, re, json
from collections import defaultdict, Counter

# ----------------------------------------------------------------------------- config
# CFG is filled from the command line in main(); the functions below read it as a module global.
CFG = {}
B_BOOT, SEED = 2000, 20261006
MIN_LEN = 30
NONPRIMARY = re.compile(r"(_alt|_random|_fix)$|^chrUn")
BINS = [0, 30, 45, 60, 90, 120, 180, 300, 100000]

# ----------------------------------------------------------------------------- genetic code
_B = "TCAG"; _AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a+b+c: _AA[16*i+4*j+k] for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}
STOPS = {"TAA", "TAG", "TGA"}
def tr(s, o): return "".join(CODON.get(s[i:i+3], "X") for i in range(o, len(s)-2, 3))
def is_open(s, o): return tr(s, o)[:-1].count("*") == 0     # no stop among complete codons, last codon excluded
def fix(cpos): return (3 - cpos) % 3                          # codon position of 5' base -> offset of first full codon

def read_fa(p):
    d = {}; k = None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): k = l[1:]; d[k] = []
        elif k: d[k].append(l.upper())
    return {k: "".join(v) for k, v in d.items()}
def parse_coord(h):
    c = h.split("::")[1]; ch, rest = c.split(":"); st = rest[-2]; s, e = map(int, rest[:-3].split("-"))
    return ch, s, e, st

def esl_orfs(seq):
    """pure-python replica of `esl-translate -l 1 --watson`: (start1, end1, frame1..3)."""
    s = seq.upper(); n = len(s); out = []
    for frame in (1, 2, 3):
        off = frame-1; st = off; L = 0; i = off
        while i+3 <= n:
            if s[i:i+3] in STOPS:
                if L: out.append((st+1, i, frame))
                L = 0; st = i+3
            else: L += 1
            i += 3
        if L: out.append((st+1, st+3*L, frame))
    return out

# ----------------------------------------------------------------------------- OLGenie counts
CLS = {"NN": 0, "NS": 1, "SN": 2, "SS": 3}
def _eff(c0, c1): return "N" if CODON.get(c0, "X") != CODON.get(c1, "X") else "S"
def olg_counts(h, m, r, t):
    L = len(h); sites = np.zeros(4); diffs = np.zeros(4)
    for i in range(max(r, t), L):
        cr = r+3*((i-r)//3); ca = t+3*((i-t)//3)
        if cr+3 > L or ca+3 > L: continue
        pr, pa = i-cr, i-ca; R = h[cr:cr+3]; A = h[ca:ca+3]
        for b in "ACGT":
            if b == h[i]: continue
            sites[CLS[_eff(R, R[:pr]+b+R[pr+1:]) + _eff(A, A[:pa]+b+A[pa+1:])]] += 1/3
        if h[i] != m[i]:
            b = m[i]
            diffs[CLS[_eff(R, R[:pr]+b+R[pr+1:]) + _eff(A, A[:pa]+b+A[pa+1:])]] += 1
    return sites, diffs

def jc(p):
    p = np.asarray(p, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(p < 0.75, -0.75*np.log(1-4*p/3), np.nan)
def omegas(S, Dd):
    """S, Dd: pooled (...,4) arrays -> dict of ratios (vectorised over leading axes)."""
    d = jc(Dd/S)
    return dict(dNN=d[..., 0], dNS=d[..., 1], dSN=d[..., 2], dSS=d[..., 3],
                w_dNNdNS=d[..., 0]/d[..., 1], w_dSNdSS=d[..., 2]/d[..., 3], w_ref=d[..., 0]/d[..., 2])
CNT = ["s_NN", "s_NS", "s_SN", "s_SS", "d_NN", "d_NS", "d_SN", "d_SS"]

# ----------------------------------------------------------------------------- D arm
def build_D():
    Dh = {h.split("::")[0].split("|")[1]: (h, s) for h, s in read_fa(CFG["D_human"]).items()}
    Dm = {h.split("::")[0].split("|")[1]: (h, s) for h, s in read_fa(CFG["D_mouse"]).items()}
    fr = pd.read_csv(CFG["D_frames"]); fr["idx"] = fr.name.str.split("|").str[1]
    fr["offs"] = fr.prot_offs.apply(ast.literal_eval)
    rows = []
    for r in fr.itertuples():
        h, s = Dh[r.idx]; ch, st, en, strand = parse_coord(h); L = len(s)
        lifted = r.idx in Dm; mlen = len(Dm[r.idx][1]) if lifted else np.nan
        ref, alt = r.offs[0], r.offs[1]
        rows.append(dict(arm="D", id=r.idx, name=r.name, chrom=ch, start=st, end=en, strand=strand, L=L,
                         n_prot_frames=len(r.offs), prot_offsets=str(r.offs), ref=ref, alt=alt, rel=(alt-ref) % 3,
                         nonprimary=bool(NONPRIMARY.search(ch)), too_short=L < MIN_LEN, lift_failed=not lifted,
                         mouse_len=mlen, indel=lifted and mlen != L,
                         mouse_coord=Dm[r.idx][0].split("::")[1] if lifted else ""))
    D = pd.DataFrame(rows)
    D["v2_scored"] = ~D.too_short & ~D.lift_failed & ~D.indel
    cnt = np.full((len(D), 8), np.nan); cnt_sw = np.full((len(D), 8), np.nan)
    for i, r in enumerate(D.itertuples()):
        if r.v2_scored:
            h = Dh[r.id][1]; m = Dm[r.id][1]
            s, d = olg_counts(h, m, r.ref, r.alt); cnt[i] = np.r_[s, d]
            s, d = olg_counts(h, m, r.alt, r.ref); cnt_sw[i] = np.r_[s, d]    # roles swapped
    D[CNT] = cnt; D[[c+"_swap" for c in CNT]] = cnt_sw
    return D

# ----------------------------------------------------------------------------- N arm
def build_N():
    Nh = {h.split("::")[0]: (h, s) for h, s in read_fa(CFG["N_human"]).items()}
    Nm = {h.split("::")[0]: (h, s) for h, s in read_fa(CFG["N_mouse"]).items()}
    ce = pd.read_csv(CFG["cds_exons"], sep="\t", header=None, names=["tx", "i", "chrom", "s", "e", "len", "strand", "name"])
    ce = ce.sort_values(["tx", "i"]); ce["cum"] = ce.groupby("tx")["len"].cumsum() - ce["len"]
    cum = dict(zip(ce.name, ce.cum)); del ce
    ov = pd.read_csv(CFG["N_cds"], sep="\t", header=None,
                     names=["c", "s", "e", "nm", "x", "st", "c2", "bs", "be", "bname", "x2", "bst", "olen"])
    ov = ov[(ov.bs <= ov.s) & (ov.be >= ov.e)].copy()                       # fully containing CDS blocks
    ov["prec"] = np.where(ov.st == "+", ov.bname.map(cum) + (ov.s - ov.bs), ov.bname.map(cum) + (ov.be - ov.e))
    ov["off"] = fix((ov.prec % 3).astype(int))
    ann = ov.groupby("nm").off.agg(lambda x: tuple(sorted(set(x))))
    # ORF pair that defines each interval in its source exon (for the 'orfpair' sensitivity rule)
    byc = defaultdict(list)
    for ph, ps in read_fa(CFG["pool"]).items():
        ch, s, e, st = parse_coord(ph); byc[(ch, st)].append((s, e, ps))
    for k in byc: byc[k].sort()
    starts = {k: [x[0] for x in v] for k, v in byc.items()}
    def orfpair(nh):
        ch, s2, e2, st = parse_coord(nh); Lst = byc.get((ch, st), []); i = bisect.bisect_right(starts.get((ch, st), []), s2)
        for j in range(i-1, max(-1, i-60), -1):
            s, e, ps = Lst[j]
            if s <= s2 and e2 <= e:
                a1, b1 = (s2-s+1, e2-s) if st == "+" else (e-e2+1, e-s2)
                o = esl_orfs(ps)
                for x in range(len(o)):
                    for y in range(x+1, len(o)):
                        if o[x][2] != o[y][2] and max(o[x][0], o[y][0]) == a1 and min(o[x][1], o[y][1]) == b1:
                            return tuple(sorted(((f-1)-(a1-1)) % 3 for f in (o[x][2], o[y][2])))
        return None
    rows = []
    for n, (h, s) in Nh.items():
        ch, st, en, strand = parse_coord(h); L = len(s)
        a = ann.get(n, ()); lifted = n in Nm; mlen = len(Nm[n][1]) if lifted else np.nan
        opens = tuple(f for f in range(3) if is_open(s, f))
        rows.append(dict(arm="N", id=n, chrom=ch, start=st, end=en, strand=strand, L=L, ann_frames=str(a),
                         n_ann=len(a), ref=a[0] if len(a) == 1 else -1, open_frames=str(opens), orfpair=str(orfpair(h)),
                         nonprimary=bool(NONPRIMARY.search(ch)), too_short=L < MIN_LEN, lift_failed=not lifted,
                         mouse_len=mlen, indel=lifted and mlen != L,
                         mouse_coord=Nm[n][0].split("::")[1] if lifted else ""))
    N = pd.DataFrame(rows)
    def alt_v2(r):
        if r.n_ann != 1: return -1
        op = ast.literal_eval(r.open_frames); a = r.ref
        if a not in op: return -1
        c = [f for f in op if f != a]
        if len(c) == 1: return c[0]
        if len(c) == 2: return -1 if a == 2 else min(c)
        return -1
    def alt_orfpair(r):
        p = ast.literal_eval(r.orfpair) if r.orfpair != "None" else None
        if r.n_ann != 1 or p is None or r.ref not in p: return -1
        return [f for f in p if f != r.ref][0]
    N["alt"] = N.apply(alt_v2, axis=1); N["alt_orfpair"] = N.apply(alt_orfpair, axis=1)
    N["rel"] = np.where(N.alt >= 0, (N.alt - N.ref) % 3, -1)
    N["rel_orfpair"] = np.where(N.alt_orfpair >= 0, (N.alt_orfpair - N.ref) % 3, -1)
    tech = ~N.too_short & ~N.lift_failed & ~N.indel
    N["v2_scored"] = tech & (N.alt >= 0); N["orfpair_scored"] = tech & (N.alt_orfpair >= 0)
    cnt = np.full((len(N), 8), np.nan); cnt2 = np.full((len(N), 8), np.nan)
    for i, r in enumerate(N.itertuples()):
        if r.v2_scored:
            s, d = olg_counts(Nh[r.id][1], Nm[r.id][1], r.ref, r.alt); cnt[i] = np.r_[s, d]
        if r.orfpair_scored:
            if r.v2_scored and r.alt_orfpair == r.alt: cnt2[i] = cnt[i]
            else:
                s, d = olg_counts(Nh[r.id][1], Nm[r.id][1], r.ref, r.alt_orfpair); cnt2[i] = np.r_[s, d]
    N[CNT] = cnt; N[[c+"_orfpair" for c in CNT]] = cnt2
    return N

# ----------------------------------------------------------------------------- statistics
def pooled(X): return omegas(X[:, :4].sum(0), X[:, 4:].sum(0))
def boot(X, B=B_BOOT, seed=SEED, chunk=250):
    """bootstrap over intervals -> dict of (B,) arrays"""
    rng = np.random.default_rng(seed); n = len(X); out = []
    for b0 in range(0, B, chunk):
        b = min(chunk, B-b0)
        W = rng.multinomial(n, np.full(n, 1/n), size=b).astype(float)   # resampling weights
        out.append(W @ X)
    T = np.vstack(out); return omegas(T[:, :4], T[:, 4:])
def ci(a): a = a[np.isfinite(a)]; return np.percentile(a, [2.5, 97.5]) if len(a) else (np.nan, np.nan)
def pboot(diff):
    diff = diff[np.isfinite(diff)]; B = len(diff)
    p = 2*min((diff <= 0).mean(), (diff >= 0).mean()); return max(p, 1/B)

def compare(XD, XN, label, seed=SEED):
    eD, eN = pooled(XD), pooled(XN); bD, bN = boot(XD, seed=seed), boot(XN, seed=seed+1)
    row = dict(comparison=label, n_D=len(XD), n_N=len(XN))
    for k in ["w_dNNdNS", "w_dSNdSS", "w_ref", "dNN", "dNS", "dSN", "dSS"]:
        for arm, e, b in (("D", eD, bD), ("N", eN, bN)):
            row[f"{arm}_{k}"] = float(e[k])
            if k.startswith("w_"):
                lo, hi = ci(b[k]); row[f"{arm}_{k}_lo"], row[f"{arm}_{k}_hi"] = lo, hi
        if k in ("w_dNNdNS", "w_dSNdSS"):
            dd = bN[k] - bD[k]; lo, hi = ci(dd)
            row[f"diff_N_minus_D_{k}"] = float(eN[k] - eD[k]); row[f"diff_{k}_lo"], row[f"diff_{k}_hi"] = lo, hi
            row[f"p_{k}"] = pboot(dd)
    return row

def length_matched(D, N, cols, k=10, B=1000, seed=SEED):
    """(a) deterministic exact-length reweighting of N to the D length distribution (D restricted to lengths
    present in N); (b) bootstrap of a k:1 exact-length matched N sample redrawn every replicate, with D resampled."""
    LD, LN = D.L.values, N.L.values; XD, XN = D[cols].values, N[cols].values
    cN = Counter(LN); mD = np.array([L in cN for L in LD]); cD = Counter(LD[mD])
    w = np.array([cD.get(L, 0)/cN[L] for L in LN])
    eD = pooled(XD[mD]); eN = omegas((XN*w[:, None])[:, :4].sum(0), (XN*w[:, None])[:, 4:].sum(0))
    byL = defaultdict(list)
    for i, L in enumerate(LN): byL[L].append(i)
    rng = np.random.default_rng(seed); XDm = XD[mD]; LDm = LD[mD]; res = defaultdict(list)
    for _ in range(B):
        bi = rng.integers(0, len(XDm), len(XDm))
        nidx = np.concatenate([rng.choice(byL[L], k) for L in LDm[bi]])
        oD = pooled(XDm[bi]); oN = pooled(XN[nidx])
        for key in ("w_dNNdNS", "w_dSNdSS"):
            res["D_"+key].append(oD[key]); res["N_"+key].append(oN[key]); res["diff_"+key].append(oN[key]-oD[key])
    row = dict(n_D_matched=int(mD.sum()), n_N_pool=len(XN), k_per_D=k)
    for key in ("w_dNNdNS", "w_dSNdSS"):
        row["D_"+key] = float(eD[key]); row["N_"+key+"_weighted"] = float(eN[key])
        for pre in ("D_", "N_", "diff_"):
            lo, hi = ci(np.array(res[pre+key])); row[pre+key+"_lo"], row[pre+key+"_hi"] = lo, hi
        row["p_"+key] = pboot(np.array(res["diff_"+key]))
    return row

def by_bin(D, N, cols, seed=SEED):
    rows = []
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        d = D[(D.L > lo) & (D.L <= hi)]; n = N[(N.L > lo) & (N.L <= hi)]
        r = dict(len_bin=f"({lo}, {hi if hi < 100000 else 'inf'}]", n_D=len(d), n_N=len(n))
        for arm, x, sd in (("D", d, seed), ("N", n, seed+1)):
            X = x[cols].values
            if len(X) == 0: continue
            e = pooled(X); b = boot(X, B=1000, seed=sd)
            for k in ("w_dNNdNS", "w_dSNdSS"):
                r[f"{arm}_{k}"] = float(e[k]); r[f"{arm}_{k}_lo"], r[f"{arm}_{k}_hi"] = ci(b[k])
        rows.append(r)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- CLI driver
# (was run_selection_v3.py + the output cells of the revision session)
INPUT_FILES = dict(D_human="dcr_regions.fa", D_mouse="D_ucsc.fa", D_frames="dcr_protein_frames.csv",
                   fates_old="dcr_exclusion_fates.csv",
                   D_cds="D_cds.ov", D_cds_any="D_cds_any.ov", N_human="N_human.fa", N_mouse="N_ucsc.fa",
                   N_cds="N_cds.ov", cds_exons="cds_exons.tsv", pool="nondcr_pool.fa")

V2 = dict(n_D=1228, D_w_dNNdNS=0.2437, D_w_dSNdSS=0.4465, D_w_ref=0.2604,
          n_N=21340, N_w_dNNdNS=0.7981, N_w_dSNdSS=0.9839, N_w_ref=0.1097,
          ss13_n_D=661, ss13_D=0.4427, ss13_n_N=11559, ss13_N=0.8269,
          ss12_n_D=567, ss12_D=0.4511, ss12_n_N=9781, ss12_N=1.1955)


def parse_args(argv=None):
    work = os.environ.get("WORK", "work")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs-dir", default=os.path.join(work, "selection_inputs"),
                    help="directory holding the 11 staged inputs (default $WORK/selection_inputs)")
    for k, fn in INPUT_FILES.items():
        ap.add_argument("--" + k.replace("_", "-"), dest=k, default=None, help=f"override path to {fn}")
    ap.add_argument("--outdir", default=os.path.join(work, "selection"))
    ap.add_argument("--boot", type=int, default=2000, help="bootstrap replicates (default 2000)")
    ap.add_argument("--seed", type=int, default=20261006)
    ap.add_argument("--skip-site-counts", action="store_true", help="do not write the 11 MB per-interval table")
    return ap.parse_args(argv)


def main(argv=None):
    global B_BOOT, SEED
    a = parse_args(argv)
    for k, fn in INPUT_FILES.items():
        CFG[k] = getattr(a, k) or os.path.join(a.inputs_dir, fn)
    missing = [CFG[k] for k in INPUT_FILES if not os.path.exists(CFG[k])]
    if missing:
        sys.exit("missing inputs:\n  " + "\n  ".join(missing))
    B_BOOT, SEED = a.boot, a.seed
    # boot()/compare() bind B and seed defaults at definition time; rebind them here
    boot.__defaults__ = (B_BOOT, SEED, 250)
    compare.__defaults__ = (SEED,)
    length_matched.__defaults__ = (10, 1000, SEED)
    by_bin.__defaults__ = (SEED,)
    os.makedirs(a.outdir, exist_ok=True)
    out = lambda f: os.path.join(a.outdir, f)

    D = build_D(); N = build_N()
    print(f"D arm: {len(D)} intervals, {int(D.v2_scored.sum())} scored (all contigs); "
          f"N arm: {len(N)} intervals, {int(N.v2_scored.sum())} scored", flush=True)

    def arm_sets(primary_only):
        d = D[D.v2_scored & (~D.nonprimary if primary_only else True)]
        n = N[N.v2_scored & (~N.nonprimary if primary_only else True)]
        return d, n

    SCOPES = (("all_contigs (v2 scope)", False), ("primary_contigs_only (v3)", True))
    # ---------------- headline table, both scopes side by side + sensitivities
    rows = []
    for scope, po in SCOPES:
        d, n = arm_sets(po)
        r = compare(d[CNT].values, n[CNT].values, "main"); r["scope"] = scope; rows.append(r)
        # sensitivity 1: principled null (ORF-pair latent frame, no ref-offset-2 drop)
        n2 = N[N.orfpair_scored & (~N.nonprimary if po else True)]
        r = compare(d[CNT].values, n2[[c+"_orfpair" for c in CNT]].values, "sens_null_orfpair_latent"); r["scope"] = scope; rows.append(r)
        # sensitivity 2: DCR frame roles swapped (reference = higher offset)
        r = compare(d[[c+"_swap" for c in CNT]].values, n[CNT].values, "sens_DCR_roles_swapped"); r["scope"] = scope; rows.append(r)
    omega = pd.DataFrame(rows); omega.insert(0, "scope", omega.pop("scope"))
    omega.to_csv(out("selection_omega_v3.csv"), index=False)

    # ---------------- by frame relation
    fr_rows = []
    for scope, po in SCOPES:
        d, n = arm_sets(po)
        for rel, lab in ((1, "ss13 (+1)"), (2, "ss12 (+2)")):
            r = compare(d[d.rel == rel][CNT].values, n[n.rel == rel][CNT].values, lab); r["scope"] = scope; fr_rows.append(r)
    framerel = pd.DataFrame(fr_rows); framerel.insert(0, "scope", framerel.pop("scope"))
    framerel.to_csv(out("selection_by_framerel_v3.csv"), index=False)

    # ---------------- length matched + bins
    lm_rows = []
    for scope, po in SCOPES:
        d, n = arm_sets(po); r = length_matched(d, n, CNT); r["scope"] = scope; lm_rows.append(r)
    lenmatch = pd.DataFrame(lm_rows); lenmatch.insert(0, "scope", lenmatch.pop("scope"))
    lenmatch.to_csv(out("selection_length_matched_v3.csv"), index=False)
    bins_tab = []
    for scope, po in SCOPES:
        d, n = arm_sets(po); b = by_bin(d, n, CNT); b.insert(0, "scope", scope); bins_tab.append(b)
    pd.concat(bins_tab).to_csv(out("selection_by_lenbin_v3.csv"), index=False)

    # ---------------- reproduction check against the published v2 numbers
    rep = []
    ma = omega[(omega.comparison == "main") & (omega.scope.str.startswith("all"))].iloc[0]
    fa_ = framerel[framerel.scope.str.startswith("all")].set_index("comparison")
    pairs = [("DCR n", ma.n_D, 1228), ("DCR dNN/dNS", ma.D_w_dNNdNS, 0.24365235977561348),
             ("DCR dNN/dNS CI lo", ma.D_w_dNNdNS_lo, 0.22331127604097556), ("DCR dNN/dNS CI hi", ma.D_w_dNNdNS_hi, 0.2640573934596554),
             ("DCR dSN/dSS", ma.D_w_dSNdSS, 0.4464720083337628), ("DCR dSN/dSS CI lo", ma.D_w_dSNdSS_lo, 0.37573987726465374),
             ("DCR dSN/dSS CI hi", ma.D_w_dSNdSS_hi, 0.529163517854424), ("DCR omega_ref", ma.D_w_ref, 0.26036179559993633),
             ("null n", ma.n_N, 21340), ("null dNN/dNS", ma.N_w_dNNdNS, 0.7980769978895574),
             ("null dSN/dSS", ma.N_w_dSNdSS, 0.9839237195317997), ("null dSN/dSS CI lo (methods doc)", ma.N_w_dSNdSS_lo, 0.951),
             ("null dSN/dSS CI hi (methods doc)", ma.N_w_dSNdSS_hi, 1.017), ("null omega_ref", ma.N_w_ref, 0.10969223218379394),
             ("ss13 DCR n", fa_.loc["ss13 (+1)", "n_D"], 661), ("ss13 DCR dSN/dSS", fa_.loc["ss13 (+1)", "D_w_dSNdSS"], 0.4426740309922797),
             ("ss13 null n", fa_.loc["ss13 (+1)", "n_N"], 11559), ("ss13 null dSN/dSS", fa_.loc["ss13 (+1)", "N_w_dSNdSS"], 0.8268923729860755),
             ("ss12 DCR n", fa_.loc["ss12 (+2)", "n_D"], 567), ("ss12 DCR dSN/dSS", fa_.loc["ss12 (+2)", "D_w_dSNdSS"], 0.4511129765316218),
             ("ss12 null n", fa_.loc["ss12 (+2)", "n_N"], 9781), ("ss12 null dSN/dSS", fa_.loc["ss12 (+2)", "N_w_dSNdSS"], 1.1955066677082553),
             ("length-matched DCR dSN/dSS", lenmatch.iloc[0].D_w_dSNdSS, 0.439),
             ("length-matched null dSN/dSS", lenmatch.iloc[0].N_w_dSNdSS_weighted, 1.034)]
    for lab, new, old in pairs:
        isn = lab.endswith(" n")
        d = new - old; rel = d / old
        ok = (abs(rel) <= 0.01) if isn else (abs(d) <= 0.02)
        rep.append(dict(quantity=lab, v2=old, reimplementation=new, diff=d, rel_diff=rel,
                        tolerance=("n within 1%" if isn else "|diff|<=0.02"), within_tolerance=ok))
    pd.DataFrame(rep).to_csv(out("selection_v2_reproduction_check.csv"), index=False)

    # ---------------- per-interval site counts + fates
    if not a.skip_site_counts:
        fo = pd.read_csv(CFG["fates_old"]); fo["id"] = fo.name.str.split("|").str[1]
        Dm_ = D.merge(fo[["id", "fate", "nf", "len_h"]].rename(columns={"fate": "fate_old_table", "nf": "nf_gencode_old"}),
                      on="id", how="left")

        def cas(r, np_first):
            if np_first and r.nonprimary: return "nonprimary_contig"
            if r.too_short: return "too_short"
            if r.lift_failed: return "lift_failed"
            if r.indel: return "indel"
            return "kept"

        def ncas(r, npf):
            if npf and r.nonprimary: return "nonprimary_contig"
            if r.too_short: return "too_short"
            if r.n_ann != 1: return "not_one_annotated_frame"
            if r.alt < 0: return "no_usable_latent_frame"
            if r.lift_failed: return "lift_failed"
            if r.indel: return "indel"
            return "kept"

        Dm_["fate_v2_repro"] = Dm_.apply(lambda r: cas(r, False), axis=1)
        Dm_["fate_v3"] = Dm_.apply(lambda r: cas(r, True), axis=1)
        N["fate_v2_repro"] = N.apply(lambda r: ncas(r, False), axis=1)
        N["fate_v3"] = N.apply(lambda r: ncas(r, True), axis=1)
        F = Dm_.set_index("id")
        Dsc = D.assign(scored_v2=D.v2_scored, scored_v3=D.v2_scored & ~D.nonprimary,
                       fate_v3=F.loc[D.id, "fate_v3"].values, fate_v2_repro=F.loc[D.id, "fate_v2_repro"].values)
        Nsc = N.assign(scored_v2=N.v2_scored, scored_v3=N.v2_scored & ~N.nonprimary)
        common = ["arm", "id", "chrom", "start", "end", "strand", "L", "nonprimary", "ref", "alt", "rel",
                  "too_short", "lift_failed", "indel", "mouse_coord", "fate_v2_repro", "fate_v3", "scored_v2", "scored_v3"] + CNT
        Dout = Dsc[common + [c+"_swap" for c in CNT] + ["prot_offsets"]]
        Nout = Nsc[common + ["ann_frames", "open_frames", "orfpair", "alt_orfpair", "rel_orfpair", "orfpair_scored"]
                   + [c+"_orfpair" for c in CNT]]
        SC = pd.concat([Dout, Nout], ignore_index=True)
        SC["rel"] = SC.rel.map({1: "ss13 (+1)", 2: "ss12 (+2)"})
        SC["len_bin"] = pd.cut(SC.L, BINS).astype(str)
        SC.to_csv(out("selection_site_counts_v3.csv"), index=False)

    m = omega[(omega.comparison == "main")].set_index("scope")
    for scope in m.index:
        r = m.loc[scope]
        print(f"{scope}: DCR n={int(r.n_D)} dSN/dSS={r.D_w_dSNdSS:.3f} [{r.D_w_dSNdSS_lo:.3f},{r.D_w_dSNdSS_hi:.3f}] | "
              f"null n={int(r.n_N)} dSN/dSS={r.N_w_dSNdSS:.3f} | p={r.p_w_dSNdSS:.3g}")


if __name__ == "__main__":
    main()
