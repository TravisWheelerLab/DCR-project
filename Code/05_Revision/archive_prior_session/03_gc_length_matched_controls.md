# 03_gc_length_matched_controls — as executed, not refactored

**GC+length-matched random-region controls and GC-matched Figure 2B curves.**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[231..247]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 7 (msg 231 … msg 243) |
| Language / environment | bash (env: python); python (env: dcr-tierA) |
| Saved outputs | `gc_matched_fig2B.png`, `gc_matched_fig2b.csv` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 1 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

---

### [msg 231] python · env: dcr-tierA · Building dual-coding calls and per-region GC

~~~~python
cols=["ORF","chr","start","end","strand","ExLength","ORFstart","ORFend","ORFlen_nt","ORF2exon_ratio","ORF_len_AA","frame","AA_seq"]
INP=f"{WS}/tierA/inputs"
def load_tr(p): return pd.read_csv(p,sep="\t",header=None,names=cols,dtype={"chr":str})
DCE=load_tr(f"{INP}/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv")
ND =load_tr(f"{INP}/fake_DCE_intervals.fa_translated.tsv")

def dc_regions(df, keys):
    d=df.copy()
    o0=d.ORFstart-1
    d['isDC']=((o0-3)<=0)&((d.ExLength-d.ORFend)<=3)
    dc=d[d.isDC].copy()
    dc['RID']=dc[keys].astype(str).agg(':'.join,axis=1)
    g=dc.groupby('RID').agg(n=('ORF','size'),ExLength=('ExLength','first')).reset_index()
    g['dcflag']=g.n>1
    return g

dce_reg=dc_regions(DCE,['chr','start','end'])
nd_reg =dc_regions(ND,['chr','start','end'])

# per-region GC keyed by chr:start:end matching RID. FASTA header ">chr16:26529230-26529366(+)"
def fasta_gc_keyed(path, strip_strand=True):
    out={}; hdr=None; cur=[]
    def flush():
        if hdr is not None:
            s=''.join(cur).upper()
            if len(s)>0:
                g=s.count('G')+s.count('C'); n=g+s.count('A')+s.count('T')
                if n>0: out[hdr]=(len(s),100.0*g/n)
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                flush(); h=line[1:].strip()
                # ">chr16:26529230-26529366(+)" -> key chr16:26529230:26529366
                h2=h.split('(')[0]  # drop strand
                c,rng=h2.split(':'); a,b=rng.split('-')
                hdr=f"{c}:{a}:{b}"; cur=[]
        flush()
    return out

gc_mouse=fasta_gc_keyed(f"{GC}/mouse_DCR_ortholog.fa")
gc_rand =fasta_gc_keyed(f"{GC}/random_ctrl.fa")
print("mouse GC keys:",len(gc_mouse)," random GC keys:",len(gc_rand))
# check RID format matches
print("sample dce RID:",dce_reg.RID.iloc[0], "| in gc_mouse:",dce_reg.RID.iloc[0] in gc_mouse)
print("sample nd RID:",nd_reg.RID.iloc[0], "| in gc_rand:",nd_reg.RID.iloc[0] in gc_rand)
~~~~

### [msg 233] bash · env: python · Comparing FASTA headers vs TSV coords

~~~~bash
WS="<OLD_WORKSPACE>"
echo "=== mouse ortholog fasta headers ==="; grep '>' "$WS/tierA/gc/mouse_DCR_ortholog.fa" | head -3
echo "=== random ctrl fasta headers ==="; grep '>' "$WS/tierA/gc/random_ctrl.fa" | head -3
echo "=== DCE translated tsv first rows (chr/start/end) ==="; head -3 "$WS/tierA/inputs/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv" | cut -f2-5
echo "=== random translated tsv ==="; head -3 "$WS/tierA/inputs/fake_DCE_intervals.fa_translated.tsv" | cut -f2-5
~~~~

### [msg 235] python · env: dcr-tierA · Fixing FASTA parser and joining GC

~~~~python
def fasta_gc_keyed(path):
    out={}; hdr=None; cur=[]
    def flush():
        if hdr is not None:
            s=''.join(cur).upper()
            if len(s)>0:
                g=s.count('G')+s.count('C'); n=g+s.count('A')+s.count('T')
                if n>0: out[hdr]=(len(s),100.0*g/n)
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                flush(); h=line[1:].strip().split('(')[0]
                c,rng=h.split(':'); a,b=rng.split('-')
                hdr=f"{c}:{a}:{b}"; cur=[]
            else:
                cur.append(line.strip())
        flush()
    return out

gc_mouse=fasta_gc_keyed(f"{GC}/mouse_DCR_ortholog.fa")
gc_rand =fasta_gc_keyed(f"{GC}/random_ctrl.fa")
print("mouse GC keys:",len(gc_mouse)," random GC keys:",len(gc_rand))
print("dce RID in gc_mouse:", dce_reg.RID.iloc[0] in gc_mouse)
print("nd RID in gc_rand:", nd_reg.RID.iloc[0] in gc_rand)
# attach GC to region tables
dce_reg['gc']=dce_reg.RID.map(lambda k: gc_mouse.get(k,(None,None))[1])
nd_reg['gc'] =nd_reg.RID.map(lambda k: gc_rand.get(k,(None,None))[1])
print("dce_reg GC coverage:", dce_reg.gc.notna().mean(), " nd_reg:", nd_reg.gc.notna().mean())
~~~~

### [msg 237] python · env: dcr-tierA · GC+length matching random control

~~~~python
# GC+length-matched resampling of the random-region control pool to the DCR orthologs' joint distribution.
# Bin on length (bp) and GC (%), sample control regions to match DCR per-bin counts.
rng=np.random.default_rng(42)
dce=dce_reg.dropna(subset=['gc']).copy()
nd =nd_reg.dropna(subset=['gc']).copy()

# joint bins: length in exon-appropriate bins, GC in 2.5% bins
len_edges=np.array([21,30,45,60,90,120,180,300,500,10000])
gc_edges=np.arange(20,92.5,5.0)
def binidx(df):
    li=np.digitize(df.ExLength, len_edges)
    gi=np.digitize(df.gc, gc_edges)
    return li*1000+gi
dce['bin']=binidx(dce); nd['bin']=binidx(nd)

# For each DCR bin, sample same number of control regions (with replacement if pool thin)
matched=[]
thin=0
for b,grp in dce.groupby('bin'):
    need=len(grp)
    pool=nd[nd.bin==b]
    if len(pool)==0:
        thin+=need; continue
    repl=len(pool)<need
    if repl: thin+=0
    idx=rng.choice(pool.index, size=need, replace=(len(pool)<need))
    matched.append(nd.loc[idx])
matched=pd.concat(matched)
print(f"DCR regions with GC: {len(dce)}; matched control regions drawn: {len(matched)}; DCR in empty control bins: {thin}")
print(f"\nDCR dual-coding: {dce.dcflag.mean()*100:.1f}%")
print(f"GC+len-matched random control dual-coding: {matched.dcflag.mean()*100:.1f}%")
print(f"(length-only baseline random control was 60.5%)")
# GC/length balance check
print(f"\nGC medians  DCR={dce.gc.median():.2f}  matched-ctrl={matched.gc.median():.2f}")
print(f"Len medians DCR={dce.ExLength.median():.0f}  matched-ctrl={matched.ExLength.median():.0f}")
~~~~

### [msg 239] python · env: dcr-tierA · NN GC-matching with bootstrap CI

~~~~python
# Robust version: nearest-neighbor matching in standardized (length, GC) space, with bootstrap CI.
from scipy.spatial import cKDTree
# standardize features on the DCR distribution
feat_cols=['ExLength','gc']
mu=dce[feat_cols].mean().values; sd=dce[feat_cols].std().values
def std(df): return (df[feat_cols].values-mu)/sd
nd_std=std(nd); dce_std=std(dce)
tree=cKDTree(nd_std)

boot=[]
rng=np.random.default_rng(7)
for it in range(200):
    # for each DCR region pick one of its k nearest control neighbors at random (k=15) w/o double-heavy reuse
    k=15
    dd,ii=tree.query(dce_std,k=k)
    pick=ii[np.arange(len(ii)),rng.integers(0,k,len(ii))]
    ctrl_dc=nd.iloc[pick].dcflag.values.mean()
    boot.append(ctrl_dc)
boot=np.array(boot)*100
print(f"DCR observed dual-coding: {dce.dcflag.mean()*100:.1f}%")
print(f"GC+length NN-matched control: {boot.mean():.1f}% (95% CI {np.percentile(boot,2.5):.1f}-{np.percentile(boot,97.5):.1f})")
print(f"Length-only baseline control: 60.5%")
# check matching quality on the k=1 NN
dd1,ii1=tree.query(dce_std,k=1)
m1=nd.iloc[ii1]
print(f"\nNN(k=1) matched control GC median={m1.gc.median():.2f} (DCR {dce.gc.median():.2f}); len median={m1.ExLength.median():.0f} (DCR {dce.ExLength.median():.0f})")
gap=dce.dcflag.mean()*100 - boot.mean()
print(f"\nGap surviving GC+length control: {gap:.1f} points (baseline gap was {dce.dcflag.mean()*100-60.5:.1f})")
~~~~

### [msg 241] python · env: dcr-tierA · Building GC-matched Fig 2B curves

~~~~python
# Build GC+length-matched control per length bin for the Fig 2B-style curve.
# For each DCR region, its k-NN matched control set; aggregate matched control dual-coding by length bin.
def make_breaks(maxlen, binlen=3):
    breaks=[]; i=1
    while i<maxlen:
        bw=binlen if i<=100 else int(np.ceil(np.sqrt(i)/2))
        breaks.append(i+bw); i+=bw
    return np.array(breaks)
maxlen=int(max(dce.ExLength.max(), nd.ExLength.max()))
breaks=make_breaks(maxlen,3)
def binprop(df,flagcol='dcflag'):
    r=df.copy(); cats=pd.cut(r.ExLength,bins=breaks)
    r['mid']=cats.apply(lambda x:x.mid if pd.notna(x) else np.nan).astype(float)
    g=r.dropna(subset=['mid']).groupby('mid')
    return g[flagcol].agg(prop='mean',n='size').reset_index()

# matched control regions: aggregate the multiset drawn by NN (k=15 random pick, one realization for the curve)
k=15; dd,ii=tree.query(dce_std,k=k)
rng=np.random.default_rng(11)
pick=ii[np.arange(len(ii)),rng.integers(0,k,len(ii))]
matched_curve=nd.iloc[pick].copy()
# assign each matched control the length of its DCR partner? No—use matched control's own length (it's GC+len matched so ~same)
dce_b=binprop(dce); nd_base_b=binprop(nd_reg.dropna(subset=['gc'])); matched_b=binprop(matched_curve)

# save table
gcmatch=pd.DataFrame({
 'comparison':['DCR_observed','random_ctrl_length_only','random_ctrl_GClength_matched'],
 'dualcoding_pct':[round(dce.dcflag.mean()*100,1),60.5,round(boot.mean(),1)],
 'ci_low':[np.nan,np.nan,round(np.percentile(boot,2.5),1)],
 'ci_high':[np.nan,np.nan,round(np.percentile(boot,97.5),1)]})
gcmatch.to_csv(f"{WS}/tierA/gc_matched_fig2b.csv",index=False)
print(gcmatch.to_string(index=False))
print(f"\nGap: baseline {round(dce.dcflag.mean()*100-60.5,1)} pts -> GC-matched {round(dce.dcflag.mean()*100-boot.mean(),1)} pts")
~~~~

### [msg 243] python · env: dcr-tierA · Rendering GC-matched Fig 2B

~~~~python
fig,axes=plt.subplots(1,2,figsize=(8.0,3.6))
# Panel A: length-binned curves (DCR vs length-only ctrl vs GC+length matched ctrl)
ax=axes[0]
for df,col,lab in [(dce_b,"#FFA800","human DCR (mouse ortholog)"),
                   (nd_base_b,"#bdbdbd","control: length-only"),
                   (matched_b,"#2c7fb8","control: GC + length matched")]:
    m=(df.mid>=21)&(df.mid<=300)
    ax.scatter(df.mid[m],df.prop[m],s=13,color=col,label=lab,edgecolor='none',zorder=3)
ax.axvline(21,color="#FFA800",ls=':',lw=1,alpha=0.6)
ax.set_xlim(21,300);ax.set_ylim(0,1.02)
ax.set_xlabel("Exon region length (bp)");ax.set_ylabel("Proportion with \u22652 full-length ORFs")
ax.set_title("Fig 2B under GC matching",loc="left")
ax.legend(frameon=False,fontsize=6.2,loc="upper right")
# Panel B: aggregate bars with CI
ax=axes[1]
labs=["DCR\nobserved","control\nlength-only","control\nGC+length"]
vals=[dce.dcflag.mean()*100,60.5,boot.mean()]
cols_=["#FFA800","#bdbdbd","#2c7fb8"]
bars=ax.bar([0,1,2],vals,color=cols_,width=0.62)
ax.errorbar(2,boot.mean(),yerr=[[boot.mean()-np.percentile(boot,2.5)],[np.percentile(boot,97.5)-boot.mean()]],
            fmt='none',ecolor='black',capsize=3,lw=1)
for i,v in enumerate(vals): ax.text(i,v+1.5,f"{v:.1f}%",ha='center',fontsize=7)
ax.set_xticks([0,1,2]);ax.set_xticklabels(labs,fontsize=7)
ax.set_ylim(0,95);ax.set_ylabel("% regions dual-coding")
ax.set_title("GC explains ~5 of 22 pts;\n~18-pt DCR excess remains",loc="left",fontsize=8)
# bracket
ax.annotate("",xy=(0,88),xytext=(2,88),arrowprops=dict(arrowstyle='-',lw=0.8))
ax.text(1,89.5,f"+{dce.dcflag.mean()*100-boot.mean():.0f} pts",ha='center',fontsize=7)
fig.suptitle("R2-M1: DCR dual-coding excess survives GC+length control",fontsize=9,x=0.02,ha='left')
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(f"{WS}/tierA/gc_matched_fig2B.png",dpi=300)
print("saved")
plt.close(fig)
~~~~
