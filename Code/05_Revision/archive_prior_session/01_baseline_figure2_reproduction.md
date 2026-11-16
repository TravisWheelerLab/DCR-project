# 01_baseline_figure2_reproduction — as executed, not refactored

**Reproduction of the published Figure 2A-2C (baseline before any control); defines helper functions reused by later steps.**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[145..205]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 25 (msg 145 … msg 199) |
| Language / environment | bash (env: python); python (env: python) |
| Saved outputs | `baseline_fig2B.png`, `baseline_fig2C.png`, `baseline_counts.csv`, `baseline_fig2b_dcr.csv`, `baseline_fig2b_random.csv`, `baseline_fig2b_shuffle.csv`, `baseline_fig2b_expected.csv` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 15 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

---

### [msg 145] python · env: python · Validating DCR bed and mouse status

~~~~python
import pandas as pd, numpy as np
INP=f"{WS}/tierA/inputs"

# --- DCR interval bed ---
bed = pd.read_csv(f"{INP}/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed", sep="\t", header=None)
print("bed shape:", bed.shape)
print("bed sample:"); print(bed.head(2).to_string())
# columns per script: chr start end DCRid score strand(gene?) ... col7 = nFrames per awk in 02
# DCRid is col4 (index 3)
dcr_ids = bed[3].astype(str).unique()
print("unique DCR ids in bed:", len(dcr_ids))
print("interval rows:", len(bed))

# --- mouse status (drives Fig 2A/2B) ---
mice = pd.read_csv(f"{INP}/DCRmiceStatus.tsv", sep=r"\s+", header=0)
print("\nDCRmiceStatus cols:", list(mice.columns), "rows:", len(mice))
print(mice['mice_status'].value_counts().to_string())
~~~~

### [msg 147] bash · env: python · Inspecting mice status file format

~~~~bash
sed -n '1,8p' "<OLD_WORKSPACE>/tierA/inputs/DCRmiceStatus.tsv" | cat -A | head -8
echo "--- unique status values (col2) ---"
awk 'NR>1{print $2}' "<OLD_WORKSPACE>/tierA/inputs/DCRmiceStatus.tsv" | sort | uniq -c
~~~~

### [msg 149] python · env: python · Reconciling mouse conservation breakdown

~~~~python
mice = pd.read_csv(f"{INP}/DCRmiceStatus.tsv", sep=r"\s+", header=0,
                   names=["DCE_ID","mice_status","group","extra"], engine="python", skiprows=1)
# some rows have "too short" (2 words) -> shift; rebuild robustly from raw
rows=[]
with open(f"{INP}/DCRmiceStatus.tsv") as f:
    header=f.readline()
    for line in f:
        p=line.split()
        did=p[0]; grp=p[-1]
        status=" ".join(p[1:-1])
        rows.append((did,status,grp))
mice=pd.DataFrame(rows,columns=["DCE_ID","mice_status","group"])
print("rows:", len(mice))
print(mice['mice_status'].value_counts().to_string())
n=len(mice); dcr=(mice.mice_status=='DCR').sum()
print(f"\nconserved-dual-coding: {dcr}/{n} = {100*dcr/n:.1f}%")
# paper: 981/1294 (80%) — 1294 = DCRs that lifted over & are testable
# here 'too short' excluded from denominator typically
testable = mice[mice.mice_status!='too short']
print(f"excluding 'too short': {dcr}/{len(testable)} = {100*dcr/len(testable):.1f}%")
~~~~

### [msg 152] python · env: python · Parsing mouse status breakdown

~~~~python
rows=[]
with open(f"{INP}/DCRmiceStatus.tsv") as f:
    header=f.readline()
    for line in f:
        p=line.split()
        did=p[0]; grp=p[-1]; status=" ".join(p[1:-1])
        rows.append((did,status,grp))
mice=pd.DataFrame(rows,columns=["DCE_ID","mice_status","group"])
print("header was:", header.strip())
print("rows:", len(mice))
print(mice['mice_status'].value_counts().to_string())
n=len(mice); dcr=(mice.mice_status=='DCR').sum()
testable=mice[mice.mice_status!='too short']
print(f"\nconserved-dual-coding all: {dcr}/{n} = {100*dcr/n:.1f}%")
print(f"excluding too-short: {dcr}/{len(testable)} = {100*dcr/len(testable):.1f}%")
~~~~

### [msg 154] bash · env: python · Reading Fig 2 figure script logic

~~~~bash
S="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/scripts/02_Conservation/Mice_conservation.r"
echo "=== how status categories & Fig 2A/2B are built ==="
grep -n "DCR\|not DCR\|not mapped\|too short\|mice_status\|Figure2\|F2\|ggsave\|geom_bar\|nrow\|percent\|/\s*sum\|table(" "$S" 2>/dev/null | head -40
~~~~

### [msg 156] bash · env: python · Reading Fig 2A/2B plotting logic

~~~~bash
S="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/scripts/02_Conservation/Mice_conservation.r"
sed -n '222,300p' "$S"
~~~~

### [msg 158] bash · env: python · Reading DND_DCE proportion function

~~~~bash
S="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/scripts/02_Conservation/Mice_conservation.r"
echo "=== lines 80-190 (the DND_DCE / >1 ORF proportion logic) ==="
sed -n '80,190p' "$S"
~~~~

### [msg 160] bash · env: python · Reading isDualCoding definition

~~~~bash
S="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/scripts/02_Conservation/Mice_conservation.r"
sed -n '32,79p' "$S"
~~~~

### [msg 164] bash · env: python · Copying control translated TSVs

~~~~bash
WS="<OLD_WORKSPACE>"
C="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation"
DEST="$WS/tierA/inputs"
# random-region control tsv
cp "$C/fake_DCE_intervals.fa_translated.tsv" "$DEST/" && echo "OK fake_DCE_intervals tsv ($(wc -l < "$DEST/fake_DCE_intervals.fa_translated.tsv") lines)"
# codon-shuffle control tsv (the renamedWithCoord one Mice_conservation.r consumes as -e)
F="$C/fake_mm39_exons.fasta_translated.fa_linearized_renamedWithCoord.tsv"
ls -lh "$F" | awk '{print "shuffle tsv size:", $5}'
cp "$F" "$DEST/" && echo "OK fake_mm39 renamedWithCoord tsv ($(wc -l < "$DEST/fake_mm39_exons.fasta_translated.fa_linearized_renamedWithCoord.tsv") lines)"
~~~~

### [msg 166] python · env: python · Loading translated ORF tables

~~~~python
import pandas as pd, numpy as np
INP=f"{WS}/tierA/inputs"
cols=["ORF","chr","start","end","strand","ExLength","ORFstart","ORFend","ORFlen_nt","ORF2exon_ratio","ORF_len_AA","frame","AA_seq"]

def load_translated(path):
    df=pd.read_csv(path,sep="\t",header=None,names=cols,dtype={"chr":str})
    return df

DCE = load_translated(f"{INP}/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv")
ND  = load_translated(f"{INP}/fake_DCE_intervals.fa_translated.tsv")
print("DCE ORFs:", len(DCE), "| random-ctrl ORFs:", len(ND))
print("DCE unique regions:", DCE.groupby(['chr','start','end']).ngroups)
print("random unique regions:", ND.groupby(['chr','start','end']).ngroups)
print(DCE.head(2).to_string())
~~~~

### [msg 168] python · env: python · Loading codon-shuffle control

~~~~python
# Load codon-shuffle control (renamedWithCoord format). Per Conservation.sh awk it has same 13 cols
E = load_translated(f"{INP}/fake_mm39_exons.fasta_translated.fa_linearized_renamedWithCoord.tsv")
print("shuffle ORFs:", len(E), "| unique regions:", E.groupby(['chr','start','end']).ngroups)
print(E.head(3).to_string())
print("\nExLength dtype:", E.ExLength.dtype, "| any NA in ORFstart/ORFend/ExLength:",
      E[['ORFstart','ORFend','ExLength']].isna().any().to_dict())
~~~~

### [msg 170] python · env: python · Computing baseline dual-coding proportions

~~~~python
def compute_dualcoding(df, region_keys):
    """Apply the exact Mice_conservation.r rule.
    isDualCoding: ORFstart-1-3<=0 AND ExLength-ORFend<=3  (R sets ORFstart-=1 first).
    A region is dual-coding if it has >1 such ORFs in >1 distinct frame? 
    NOTE: For Fig2B (OpE_df) the rule is ORFs>1 (count of isDualCoding ORFs), 
    but the mouse-status (Fig2A) uses frames>1. Fig2B uses ORFs>1 count of dualcoding ORFs."""
    d = df.copy()
    # R: Exons$ORFstart <- Exons$ORFstart - 1 ; then isDualCoding uses that
    orfstart0 = d['ORFstart'] - 1
    d['isDualCoding'] = ((orfstart0 - 3) <= 0) & ((d['ExLength'] - d['ORFend']) <= 3)
    dc = d[d['isDualCoding']].copy()
    # region id
    dc['RID'] = dc[region_keys].astype(str).agg(':'.join, axis=1)
    grp = dc.groupby('RID')
    per_region = grp.agg(n_dc_orfs=('ORF','size'),
                         n_frames=('frame','nunique'),
                         ExLength=('ExLength','first')).reset_index()
    per_region['is_dualcoding_region'] = per_region['n_dc_orfs'] > 1  # Fig2B rule: ORFs>1
    return per_region

# DCR mouse orthologs: filter to DCE_ID regions with human sub-region >=21 handled upstream; 
# here group by chr:start:end
dce_reg = compute_dualcoding(DCE, ['chr','start','end'])
nd_reg  = compute_dualcoding(ND,  ['chr','start','end'])
e_reg   = compute_dualcoding(E,   ['chr'])   # shuffle: region id = chr col

for nm,r in [("DCR-ortholog",dce_reg),("random-ctrl",nd_reg),("shuffle-ctrl",e_reg)]:
    tot=len(r); dc=r['is_dualcoding_region'].sum()
    print(f"{nm:14s} regions={tot:7d}  dual-coding(>1 full ORF)={dc:7d} ({100*dc/tot:.1f}%)")
~~~~

### [msg 172] python · env: python · Reproducing binned Fig 2B curves

~~~~python
def make_breaks(maxlen, binlen=3):
    breaks=[]; i=1
    while i < maxlen:
        if i<=100:
            bw=binlen
        else:
            bw=int(np.ceil(np.sqrt(i)/2))
        breaks.append(i+bw); i+=bw
    return np.array(breaks)

maxlen = int(max(dce_reg.ExLength.max(), nd_reg.ExLength.max(), e_reg.ExLength.max()))
breaks = make_breaks(maxlen, binlen=3)
# midpoints of consecutive breaks (pd.cut style)
def binned_prop(reg):
    r=reg.copy()
    cats=pd.cut(r.ExLength, bins=breaks)
    r['mid']=cats.apply(lambda x: x.mid if pd.notna(x) else np.nan).astype(float)
    g=r.dropna(subset=['mid']).groupby('mid')
    out=g['is_dualcoding_region'].agg(prop='mean', n='size').reset_index()
    return out

dce_b=binned_prop(dce_reg); nd_b=binned_prop(nd_reg); e_b=binned_prop(e_reg)
# theoretical expectation curve: 1-(1-(61/64)^(L/3))^2
Lgrid=np.arange(1,5001)
pexp=1-(1-(61/64)**(Lgrid/3))**2
exp_df=pd.DataFrame({'ExLength':Lgrid,'p':pexp})
exp_df['mid']=pd.cut(exp_df.ExLength,bins=breaks).apply(lambda x:x.mid if pd.notna(x) else np.nan).astype(float)
exp_b=exp_df.dropna(subset=['mid']).groupby('mid')['p'].mean().reset_index()

print("DCR curve bins (21-500):")
print(dce_b[(dce_b.mid>=21)&(dce_b.mid<=120)].head(10).to_string())
# save curves
for nm,df in [("dcr",dce_b),("random",nd_b),("shuffle",e_b),("expected",exp_b.rename(columns={'p':'prop'}))]:
    df.to_csv(f"{WS}/tierA/baseline_fig2b_{nm}.csv",index=False)
print("\nsaved 4 baseline curve CSVs")
~~~~

### [msg 176] python · env: python · Rendering baseline Fig 2B

~~~~python
import matplotlib as mpl, matplotlib.pyplot as plt
apply_figure_style()

# Baseline Fig 2B: proportion dual-coding vs region length, three groups
fig, ax = plt.subplots(figsize=(5.2,3.8))
groups = [("dcr", dce_b, "#FFA800", "human DCR (mouse ortholog)"),
          ("random", nd_b, "#9a9a9a", "control: random regions"),
          ("shuffle", e_b, "#4a4a4a", "control: shuffled codons")]
for key,df,col,lab in groups:
    m=(df.mid>=21)&(df.mid<=500)
    ax.scatter(df.mid[m], df.prop[m], s=14, color=col, label=lab, zorder=3, edgecolor='none')
# expected theoretical curve
m=(exp_b.mid>=21)&(exp_b.mid<=500)
ax.plot(exp_b.mid[m], exp_b.p[m], color="#7fb3d5", lw=1.4, ls="--", label="null: random codons (analytic)", zorder=2)
ax.axvline(21, color="#FFA800", lw=1, ls=":", alpha=0.7)
ax.set_xlim(21,500); ax.set_ylim(0,1.02)
ax.set_xlabel("Exon region length (bp)")
ax.set_ylabel("Proportion of regions with \u22652 full-length ORFs")
ax.set_title("Baseline Fig 2B: dual-coding potential vs region length", loc="left", fontsize=8)
ax.legend(frameon=False, fontsize=6.2, loc="lower right")
ax.margins(0.02)
fig.tight_layout()
fig.savefig(f"{WS}/tierA/baseline_fig2B.png", dpi=300)
print("saved baseline_fig2B.png")
# quick bbox check
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
ov=[(a,b) for i,(a,ba) in enumerate(texts) for b,bb in texts[i+1:] if ba.overlaps(bb)]
print("text overlaps:", len(ov))
plt.close(fig)
~~~~

### [msg 178] python · env: python · Repositioning legend on baseline Fig 2B

~~~~python
fig, ax = plt.subplots(figsize=(5.2,3.8))
for key,df,col,lab in groups:
    m=(df.mid>=21)&(df.mid<=500)
    ax.scatter(df.mid[m], df.prop[m], s=14, color=col, label=lab, zorder=3, edgecolor='none')
m=(exp_b.mid>=21)&(exp_b.mid<=500)
ax.plot(exp_b.mid[m], exp_b.p[m], color="#7fb3d5", lw=1.4, ls="--", label="null: random codons (analytic)", zorder=2)
ax.axvline(21, color="#FFA800", lw=1, ls=":", alpha=0.7)
ax.set_xlim(21,500); ax.set_ylim(0,1.02)
ax.set_xlabel("Exon region length (bp)")
ax.set_ylabel("Proportion with \u22652 full-length ORFs")
ax.set_title("Baseline (length-only): DCR orthologs retain dual-coding\npotential above both composition controls", loc="left", fontsize=8)
ax.legend(frameon=False, fontsize=6.2, loc="upper right", bbox_to_anchor=(1.0,0.98))
ax.margins(0.02)
fig.tight_layout()
fig.savefig(f"{WS}/tierA/baseline_fig2B.png", dpi=300)
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
tls={ax: set(ax.get_xticklabels()+ax.get_yticklabels()) for ax in fig.axes}
ov=[(a,b) for i,(a,ba) in enumerate(texts) for b,bb in texts[i+1:] if ba.overlaps(bb)]
print("text-text overlaps:", len(ov))
plt.close(fig)
~~~~

### [msg 180] python · env: python · Diagnosing flagged overlaps

~~~~python
fig, ax = plt.subplots(figsize=(5.2,3.8))
for key,df,col,lab in groups:
    m=(df.mid>=21)&(df.mid<=500)
    ax.scatter(df.mid[m], df.prop[m], s=14, color=col, label=lab, zorder=3, edgecolor='none')
m=(exp_b.mid>=21)&(exp_b.mid<=500)
ax.plot(exp_b.mid[m], exp_b.p[m], color="#7fb3d5", lw=1.4, ls="--", label="null: random codons (analytic)", zorder=2)
ax.axvline(21, color="#FFA800", lw=1, ls=":", alpha=0.7)
ax.set_xlim(21,500); ax.set_ylim(0,1.02)
ax.set_xlabel("Exon region length (bp)"); ax.set_ylabel("Proportion with \u22652 full-length ORFs")
ax.set_title("Baseline (length-only): DCR orthologs retain dual-coding\npotential above both composition controls", loc="left", fontsize=8)
ax.legend(frameon=False, fontsize=6.2, loc="upper right", bbox_to_anchor=(1.0,0.98))
ax.margins(0.02); fig.tight_layout()
fig.savefig(f"{WS}/tierA/baseline_fig2B.png", dpi=300)
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
ov=[(a.get_text(),b.get_text()) for i,(a,ba) in enumerate(texts) for b,bb in texts[i+1:] if ba.overlaps(bb)]
for o in ov: print(o)
plt.close(fig)
~~~~

### [msg 182] python · env: python · Loading per-DCE phyloP scores

~~~~python
from scipy import stats
# Fig 2C: mean phyloP per DCE vs non-DCR control classes. Read meanPP_per_DCE + the big out file
pp = pd.read_csv(f"{INP}/meanPP_per_DCE.tsv", sep="\t", header=None)
print("meanPP_per_DCE shape:", pp.shape)
print(pp.head(3).to_string())
~~~~

### [msg 184] bash · env: python · Inspecting phyloP output structure

~~~~bash
C="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation"
echo "=== DCEvsNonDCE.meanPhyloP100.out head ==="
head -5 "$C/DCEvsNonDCE.meanPhyloP100.out"
echo "cols: $(head -1 "$C/DCEvsNonDCE.meanPhyloP100.out" | awk '{print NF}')"
echo "lines: $(wc -l < "$C/DCEvsNonDCE.meanPhyloP100.out")"
echo ""
echo "=== phyloP section of 02_Conservation.sh ==="
grep -n "phyloP\|phylo\|bigWig\|bwtool\|DCEvsNonDCE\|meanPhyloP\|Rscript\|\.r " "<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/scripts/02_Conservation.sh" | tail -30
~~~~

### [msg 186] bash · env: python · Reading phyloP boxplot script and classes

~~~~bash
C="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation"
echo "=== classes present in the phyloP out file ==="
cut -f3 "$C/DCEvsNonDCE.meanPhyloP100.out" | sort | uniq -c
echo ""
echo "=== PhyloP_boxplots.r ==="
cat "<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/scripts/02_Conservation/PhyloP_boxplots.r"
~~~~

### [msg 188] python · env: python · Loading phyloP data and size restriction

~~~~python
# Load the full phyloP out file (300918 exons) — copy needed since it's in DOOVER
import subprocess
src=f"<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation/DCEvsNonDCE.meanPhyloP100.out"
subprocess.run(["cp",src,f"{INP}/DCEvsNonDCE.meanPhyloP100.out"],check=True)
PP=pd.read_csv(f"{INP}/DCEvsNonDCE.meanPhyloP100.out",sep="\t",header=None,
               names=["exon","meanPP","type","len"])
print("phyloP rows:", len(PP))
print(PP.type.value_counts().to_string())
# size restriction as in script
lowDCE = PP[PP.type=="DCE_exon"].len.min()+10
maxDCE = PP[PP.type=="DCE_exon"].len.max()
print(f"\nDCE size range for restriction: [{lowDCE}, {maxDCE}]")
PPs = PP[(PP.len>=lowDCE)&(PP.len<=maxDCE)].copy()
print("after size restriction:", len(PPs))
print(PPs.groupby('type').meanPP.agg(['count','mean','median']).to_string())
~~~~

### [msg 190] bash · env: python · Copying phyloP out file

~~~~bash
WS="<OLD_WORKSPACE>"
C="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation"
cp "$C/DCEvsNonDCE.meanPhyloP100.out" "$WS/tierA/inputs/" && echo "copied ($(wc -l < "$WS/tierA/inputs/DCEvsNonDCE.meanPhyloP100.out") lines)"
~~~~

### [msg 192] python · env: python · Reproducing Fig 2C phyloP tests

~~~~python
PP=pd.read_csv(f"{INP}/DCEvsNonDCE.meanPhyloP100.out",sep="\t",header=None,
               names=["exon","meanPP","type","len"])
print("phyloP rows:", len(PP))
print(PP.type.value_counts().to_string())
lowDCE = PP[PP.type=="DCE_exon"].len.min()+10
maxDCE = PP[PP.type=="DCE_exon"].len.max()
PPs = PP[(PP.len>=lowDCE)&(PP.len<=maxDCE)].copy()
print(f"\nsize range [{lowDCE},{maxDCE}] -> {len(PPs)} exons")
print(PPs.groupby('type').meanPP.agg(['count','mean','median']).to_string())

# Reproduce the paper's pairwise t-tests (overall, matching stat_compare_means t.test)
from scipy import stats
D  = PPs[PPs.type=="DCE_exon"].meanPP
Ds = PPs[PPs.type=="nDCE_exon_DCE_gene"].meanPP
ND = PPs[PPs.type=="nDCE_exon_nDCE_gene"].meanPP
for a,b,na,nb in [(D,Ds,"D","D*"),(Ds,ND,"D*","ND"),(D,ND,"D","ND")]:
    t,p=stats.ttest_ind(a,b,equal_var=False)
    print(f"t-test {na} vs {nb}: t={t:.3f} p={p:.3e}  (median {a.median():.3f} vs {b.median():.3f})")
~~~~

### [msg 194] python · env: python · Cleaning phyloP and retesting

~~~~python
PP['meanPP']=pd.to_numeric(PP['meanPP'],errors='coerce')
print("NA phyloP:", PP.meanPP.isna().sum(), "of", len(PP))
PPc=PP.dropna(subset=['meanPP'])
lowDCE = PPc[PPc.type=="DCE_exon"].len.min()+10
maxDCE = PPc[PPc.type=="DCE_exon"].len.max()
PPs = PPc[(PPc.len>=lowDCE)&(PPc.len<=maxDCE)].copy()
D  = PPs[PPs.type=="DCE_exon"].meanPP
Ds = PPs[PPs.type=="nDCE_exon_DCE_gene"].meanPP
ND = PPs[PPs.type=="nDCE_exon_nDCE_gene"].meanPP
print(f"\nn: D={len(D)} D*={len(Ds)} ND={len(ND)}")
print(f"medians: D={D.median():.3f} D*={Ds.median():.3f} ND={ND.median():.3f}")
print(f"means:   D={D.mean():.3f} D*={Ds.mean():.3f} ND={ND.mean():.3f}")
import numpy as np
for a,b,na,nb in [(D,Ds,"D","D*"),(Ds,ND,"D*","ND"),(D,ND,"D","ND")]:
    t,p=stats.ttest_ind(a.values,b.values,equal_var=False)
    u,pu=stats.mannwhitneyu(a.values,b.values,alternative='two-sided')
    print(f"{na} vs {nb}: t-test p={p:.3e} | MWU p={pu:.3e}")
~~~~

### [msg 196] python · env: python · Rendering baseline Fig 2C and counts

~~~~python
# Fig 2C baseline boxplot: phyloP by exon class (D / D* / ND), matching paper palette
apply_figure_style()
order=["DCE_exon","nDCE_exon_DCE_gene","nDCE_exon_nDCE_gene"]
labs={"DCE_exon":"D\nDCR exon","nDCE_exon_DCE_gene":"D*\nnon-DCR exon,\nDCR gene","nDCE_exon_nDCE_gene":"ND\nnon-DCR exon,\nnon-DCR gene"}
cols={"DCE_exon":"#FFA800","nDCE_exon_DCE_gene":"#FFF0BD","nDCE_exon_nDCE_gene":"#cccccc"}
fig, ax = plt.subplots(figsize=(4.2,3.6))
data=[PPs[PPs.type==t].meanPP.values for t in order]
bp=ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.6,
              medianprops=dict(color="black",lw=1.2))
for patch,t in zip(bp['boxes'],order): patch.set_facecolor(cols[t]); patch.set_alpha(0.7)
# jitter DCR points
import numpy as np
xj=np.random.normal(1,0.06,len(data[0]))
ax.scatter(xj,data[0],s=3,color="#FFA800",alpha=0.35,zorder=3,edgecolor='none')
ax.set_xticks([1,2,3]); ax.set_xticklabels([labs[t] for t in order],fontsize=6.2)
ax.set_ylabel("mean phyloP (100-way)")
ax.set_ylim(-3,7)
for i,t in enumerate(order):
    ax.text(i+1,-2.6,f"n={len(data[i])}",ha='center',fontsize=6)
ax.set_title("Baseline Fig 2C: DCR exons less conserved\nthan non-DCR controls (length-restricted)",loc="left",fontsize=8)
fig.tight_layout(); fig.savefig(f"{WS}/tierA/baseline_fig2C.png",dpi=300)
print("saved baseline_fig2C.png")
plt.close(fig)

# baseline counts CSV
baseline=pd.DataFrame({
 "metric":["intervals","DCRs","mouse_conserved_DCR","mouse_not_DCR","mouse_not_mapped","mouse_too_short",
           "fig2b_DCR_pct","fig2b_random_pct","fig2b_shuffle_pct",
           "phyloP_D_median","phyloP_Dstar_median","phyloP_ND_median",
           "phyloP_D_mean","phyloP_Dstar_mean","phyloP_ND_mean"],
 "value":[1741,1407,981,250,112,63,
          round(100*dce_reg.is_dualcoding_region.mean(),1),
          round(100*nd_reg.is_dualcoding_region.mean(),1),
          round(100*e_reg.is_dualcoding_region.mean(),1),
          round(D.median(),3),round(Ds.median(),3),round(ND.median(),3),
          round(D.mean(),3),round(Ds.mean(),3),round(ND.mean(),3)]})
baseline.to_csv(f"{WS}/tierA/baseline_counts.csv",index=False)
print(baseline.to_string(index=False))
~~~~

### [msg 199] python · env: python · Fixing n-label placement Fig 2C

~~~~python
apply_figure_style()
fig, ax = plt.subplots(figsize=(4.2,3.6))
bp=ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.6,
              medianprops=dict(color="black",lw=1.2))
for patch,t in zip(bp['boxes'],order): patch.set_facecolor(cols[t]); patch.set_alpha(0.7)
xj=np.random.normal(1,0.06,len(data[0]))
ax.scatter(xj,data[0],s=3,color="#FFA800",alpha=0.35,zorder=3,edgecolor='none')
ax.set_xticks([1,2,3]); ax.set_xticklabels([labs[t] for t in order],fontsize=6.2)
ax.set_ylabel("mean phyloP (100-way)")
ax.set_ylim(-3.2,7)
for i,t in enumerate(order):
    ax.text(i+1,6.6,f"n={len(data[i])}",ha='center',fontsize=6)
ax.set_title("Baseline Fig 2C: DCR exons less conserved\nthan non-DCR controls (length-restricted)",loc="left",fontsize=8)
fig.tight_layout(); fig.savefig(f"{WS}/tierA/baseline_fig2C.png",dpi=300)
plt.close(fig)
print("re-saved")
~~~~
