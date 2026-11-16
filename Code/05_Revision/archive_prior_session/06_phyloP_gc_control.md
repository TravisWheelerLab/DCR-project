# 06_phyloP_gc_control — as executed, not refactored

**phyloP conservation with GC control (matched resampling + linear model).**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[280..301]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 7 (msg 280 … msg 297) |
| Language / environment | python (env: dcr-tierA) |
| Saved outputs | `phyloP_gc_control.png`, `phyloP_gc_control.csv` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 3 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

Environment operations in this span (not code cells):

* `[msg 287] manage_packages: {"mode": "install", "environment": "dcr-tierA", "packages": ["statsmodels"]}`

---

### [msg 280] python · env: dcr-tierA · Joining GC to phyloP table

~~~~python
# phyloP table
pp=pd.read_csv(f"{WS}/tierA/inputs/DCEvsNonDCE.meanPhyloP100.out",sep="\t",header=None,
               names=['exonID','meanPP','cls','len'],dtype={'exonID':str})
pp['meanPP']=pd.to_numeric(pp.meanPP,errors='coerce')
pp=pp.dropna(subset=['meanPP'])
print("phyloP rows:",len(pp),"| classes:",pp.cls.value_counts().to_dict())
# Build GC lookup by coordinate from the extracted FASTAs. Headers are "chr:start-end".
# The exonID encodes coords: ..._cds_5_0_chr1_100022386_f  -> chr1, start(0-based+1?), need matching.
# Simpler: recompute GC per exon by re-extracting via the class bed files which we still have,
# keyed by exonID (col4 in the bed). Build exonID->GC from the FASTA + bed pairing.
import subprocess
# We have the 3 human FASTAs keyed by coordinate, and the beds with exonID in col4.
def bed_gc_map(bedfile, fafile):
    # bed col4=exonID, coord chr:start-end matches FASTA header
    coord2gc={}
    for h,s in read_fasta_list(fafile):
        coord2gc[h]=gcval(s)
    m={}
    with open(bedfile) as fh:
        for line in fh:
            p=line.rstrip("\n").split("\t")
            c,a,b=p[0],p[1],p[2]; eid=p[3]
            key=f"{c}:{a}-{b}"
            if key in coord2gc: m[eid]=coord2gc[key]
    return m
Cdir=f"{WS}/tierA/gc"
gmap={}
gmap.update(bed_gc_map(f"{Cdir}/DCR_exons_DCR_genes_UNIQ.b3.bed",f"{Cdir}/DCR_exons_DCR_genes_UNIQ.fa"))
gmap.update(bed_gc_map(f"{Cdir}/nonDCR_exon_DCR_gene_UNIQ.b3.bed",f"{Cdir}/nonDCR_exon_DCR_gene_UNIQ.fa"))
gmap.update(bed_gc_map(f"{Cdir}/nonDCR_exon_nonDCR_gene_UNIQ.b3.bed",f"{Cdir}/nonDCR_exon_nonDCR_gene_UNIQ.fa"))
print("GC map size:",len(gmap))
pp['gc']=pp.exonID.map(gmap)
print("GC coverage on phyloP table:",round(pp.gc.notna().mean(),3))
~~~~

### [msg 282] python · env: dcr-tierA · Rejoining GC via original beds

~~~~python
C="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation"
def bed_gc_map(bedfile, fafile):
    coord2gc={h:gcval(s) for h,s in read_fasta_list(fafile)}
    m={}
    with open(bedfile) as fh:
        for line in fh:
            p=line.rstrip("\n").split("\t")
            if len(p)<4: continue
            key=f"{p[0]}:{p[1]}-{p[2]}"; eid=p[3]
            if key in coord2gc: m[eid]=coord2gc[key]
    return m
gmap={}
for bed,fa in [("DCR_exons_DCR_genes_UNIQ","DCR_exons_DCR_genes_UNIQ"),
               ("nonDCR_exon_DCR_gene_UNIQ","nonDCR_exon_DCR_gene_UNIQ"),
               ("nonDCR_exon_nonDCR_gene_UNIQ","nonDCR_exon_nonDCR_gene_UNIQ")]:
    gmap.update(bed_gc_map(f"{C}/{bed}.bed",f"{WS}/tierA/gc/{fa}.fa"))
print("GC map size:",len(gmap))
pp['gc']=pp.exonID.map(gmap)
print("GC coverage:",round(pp.gc.notna().mean(),3))
print(pp.groupby('cls').gc.apply(lambda x:round(x.notna().mean(),3)).to_dict())
~~~~

### [msg 284] python · env: dcr-tierA · GC-matched phyloP resampling test

~~~~python
import statsmodels.formula.api as smf
# restrict to DCR size range (as pipeline does): lowDCE=min+10, maxDCE=max over DCE class
dce_len=pp[pp.cls=='DCE_exon'].len
low=dce_len.min()+10; high=dce_len.max()
ppr=pp[(pp.len>=low)&(pp.len<=high)].copy()
print(f"size restrict [{low},{high}] -> {len(ppr)} rows")
D=ppr[ppr.cls=='DCE_exon']; Dg=ppr[ppr.cls=='nDCE_exon_DCE_gene']; Ng=ppr[ppr.cls=='nDCE_exon_nDCE_gene']

# baseline (no GC control) medians/means for reference
print("\nBASELINE means: D=%.3f Dg=%.3f Ng=%.3f"%(D.meanPP.mean(),Dg.meanPP.mean(),Ng.meanPP.mean()))

# (a) GC+length-matched resampling: match Dg and Ng to D on (len,GC), compare mean phyloP
Dfeat=D[['len','gc']].values
def matched_meanPP(pool,seed=1,k=1):
    mu=Dfeat.mean(0);sd=Dfeat.std(0)
    tree=cKDTree((pool[['len','gc']].values-mu)/sd)
    dd,ii=tree.query((Dfeat-mu)/sd,k=k)
    idx=ii if k==1 else ii[:,0]
    sub=pool.iloc[idx]
    return sub.meanPP.values, sub[['len','gc']].values
mppDg,fDg=matched_meanPP(Dg); mppNg,fNg=matched_meanPP(Ng)
from scipy import stats
print("\n(a) GC+LENGTH-MATCHED comparison (vs DCR):")
print(f"  DCR mean phyloP           = {D.meanPP.mean():.3f} (median {D.meanPP.median():.3f})")
print(f"  matched nDCE_DCE_gene     = {mppDg.mean():.3f} (median {np.median(mppDg):.3f})  [matched GC {np.median(fDg[:,1]):.1f} vs {D.gc.median():.1f}]")
print(f"  matched nDCE_nDCE_gene    = {mppNg.mean():.3f} (median {np.median(mppNg):.3f})")
print(f"  D vs matched-Dg MWU p = {stats.mannwhitneyu(D.meanPP,mppDg)[1]:.2e}")
print(f"  D vs matched-Ng MWU p = {stats.mannwhitneyu(D.meanPP,mppNg)[1]:.2e}")
~~~~

### [msg 290] python · env: dcr-tierA · Rebuilding phyloP+GC table

~~~~python
import statsmodels.formula.api as smf
# rebuild after kernel restart
import numpy as np, pandas as pd, sys, random
from scipy import stats
from scipy.spatial import cKDTree
WS="<OLD_WORKSPACE>"
GC=f"{WS}/tierA/gc"
sys.path.insert(0,f"{WS}/tierA")
from esl_translate_replica import esl_translate_watson
def read_fasta_list(path):
    out=[];hdr=None;cur=[]
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if hdr is not None:out.append((hdr,''.join(cur).upper()))
                hdr=line[1:].strip();cur=[]
            else:cur.append(line.strip())
        if hdr is not None:out.append((hdr,''.join(cur).upper()))
    return out
def gcval(s):
    g=s.count('G')+s.count('C');n=g+s.count('A')+s.count('T');return 100.0*g/n if n else np.nan
pp=pd.read_csv(f"{WS}/tierA/inputs/DCEvsNonDCE.meanPhyloP100.out",sep="\t",header=None,
               names=['exonID','meanPP','cls','len'],dtype={'exonID':str})
pp['meanPP']=pd.to_numeric(pp.meanPP,errors='coerce');pp=pp.dropna(subset=['meanPP'])
C="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation"
def bed_gc_map(bedfile,fafile):
    coord2gc={h:gcval(s) for h,s in read_fasta_list(fafile)}
    m={}
    with open(bedfile) as fh:
        for line in fh:
            p=line.rstrip("\n").split("\t")
            if len(p)<4:continue
            key=f"{p[0]}:{p[1]}-{p[2]}"
            if key in coord2gc:m[p[3]]=coord2gc[key]
    return m
gmap={}
for b in ["DCR_exons_DCR_genes_UNIQ","nonDCR_exon_DCR_gene_UNIQ","nonDCR_exon_nonDCR_gene_UNIQ"]:
    gmap.update(bed_gc_map(f"{C}/{b}.bed",f"{GC}/{b}.fa"))
pp['gc']=pp.exonID.map(gmap)
print("rows",len(pp),"GC cov",round(pp.gc.notna().mean(),3))
~~~~

### [msg 293] python · env: dcr-tierA · GC-matched phyloP resampling test

~~~~python
dce_len=pp[pp.cls=='DCE_exon'].len
low=dce_len.min()+10; high=dce_len.max()
ppr=pp[(pp.len>=low)&(pp.len<=high)].copy()
D=ppr[ppr.cls=='DCE_exon']; Dg=ppr[ppr.cls=='nDCE_exon_DCE_gene']; Ng=ppr[ppr.cls=='nDCE_exon_nDCE_gene']
Dfeat=D[['len','gc']].values
def matched_meanPP(pool):
    mu=Dfeat.mean(0);sd=Dfeat.std(0)
    tree=cKDTree((pool[['len','gc']].values-mu)/sd)
    dd,ii=tree.query((Dfeat-mu)/sd,k=1)
    sub=pool.iloc[ii]; return sub.meanPP.values, sub[['len','gc']].values
mppDg,fDg=matched_meanPP(Dg); mppNg,fNg=matched_meanPP(Ng)
print(f"size restrict [{low},{high}] -> {len(ppr)} rows")
print("\nBASELINE means: D=%.3f Dg=%.3f Ng=%.3f"%(D.meanPP.mean(),Dg.meanPP.mean(),Ng.meanPP.mean()))
print("BASELINE medians: D=%.3f Dg=%.3f Ng=%.3f"%(D.meanPP.median(),Dg.meanPP.median(),Ng.meanPP.median()))
print("\n(a) GC+LENGTH-MATCHED:")
print(f"  DCR              mean {D.meanPP.mean():.3f} median {D.meanPP.median():.3f} (GC {D.gc.median():.1f})")
print(f"  matched nDCE/DCE mean {mppDg.mean():.3f} median {np.median(mppDg):.3f} (GC {np.median(fDg[:,1]):.1f})")
print(f"  matched nDCE/nDCE mean {mppNg.mean():.3f} median {np.median(mppNg):.3f} (GC {np.median(fNg[:,1]):.1f})")
print(f"  D vs matched-Dg  MWU p={stats.mannwhitneyu(D.meanPP,mppDg)[1]:.2e}")
print(f"  D vs matched-Ng  MWU p={stats.mannwhitneyu(D.meanPP,mppNg)[1]:.2e}")
~~~~

### [msg 295] python · env: dcr-tierA · Fitting phyloP linear model with GC

~~~~python
# (b) Linear model: meanPP ~ is_DCR + GC + length, on size-restricted data
ppr=ppr.copy()
ppr['is_DCR']=(ppr.cls=='DCE_exon').astype(int)
ppr['in_DCR_gene']=(ppr.cls.isin(['DCE_exon','nDCE_exon_DCE_gene'])).astype(int)
m1=smf.ols("meanPP ~ is_DCR + gc + len", data=ppr).fit()
print("Model: meanPP ~ is_DCR + GC + length")
print(m1.params.round(4).to_string())
print("is_DCR p =", f"{m1.pvalues['is_DCR']:.2e}", "| GC coef p =", f"{m1.pvalues['gc']:.2e}")
print(f"\nis_DCR coefficient: {m1.params['is_DCR']:.3f} phyloP units "
      f"(95% CI {m1.conf_int().loc['is_DCR',0]:.3f} to {m1.conf_int().loc['is_DCR',1]:.3f})")
print("=> DCR exons are LESS conserved by ~%.2f phyloP units even after adjusting for GC and length"%abs(m1.params['is_DCR']))
# save table
gcpp=pd.DataFrame({
 'test':['baseline_D_mean','baseline_Dg_mean','baseline_Ng_mean',
         'GCmatched_D_mean','GCmatched_Dg_mean','GCmatched_Ng_mean',
         'D_vs_matchedDg_MWU_p','D_vs_matchedNg_MWU_p',
         'LM_is_DCR_coef','LM_is_DCR_p','LM_gc_coef','LM_gc_p'],
 'value':[round(D.meanPP.mean(),3),round(Dg.meanPP.mean(),3),round(Ng.meanPP.mean(),3),
          round(D.meanPP.mean(),3),round(mppDg.mean(),3),round(mppNg.mean(),3),
          f"{stats.mannwhitneyu(D.meanPP,mppDg)[1]:.2e}",f"{stats.mannwhitneyu(D.meanPP,mppNg)[1]:.2e}",
          round(m1.params['is_DCR'],3),f"{m1.pvalues['is_DCR']:.2e}",
          round(m1.params['gc'],4),f"{m1.pvalues['gc']:.2e}"]})
gcpp.to_csv(f"{WS}/tierA/phyloP_gc_control.csv",index=False)
print("\nsaved phyloP_gc_control.csv")
~~~~

### [msg 297] python · env: dcr-tierA · Rendering phyloP GC-control figure

~~~~python
import matplotlib as mpl, matplotlib.pyplot as plt
mpl.rcParams.update({'font.size':8,'axes.titlesize':8,'axes.labelsize':8,'xtick.labelsize':7,
 'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':110,
 'font.family':'sans-serif','axes.linewidth':0.8,'legend.fontsize':6.5})
fig,axes=plt.subplots(1,2,figsize=(8.2,3.7))
# Panel A: baseline vs GC-matched means
ax=axes[0]
groups=['DCR\nexon','non-DCR exon,\nDCR gene','non-DCR exon,\nnon-DCR gene']
base=[D.meanPP.mean(),Dg.meanPP.mean(),Ng.meanPP.mean()]
matched=[D.meanPP.mean(),mppDg.mean(),mppNg.mean()]
x=np.arange(3);w=0.38
ax.bar(x-w/2,base,w,color='#bdbdbd',label='baseline')
ax.bar(x+w/2,matched,w,color='#2c7fb8',label='GC+length matched to DCR')
ax.set_xticks(x);ax.set_xticklabels(groups,fontsize=6.8)
ax.set_ylabel("mean phyloP (100-way)");ax.set_ylim(0,3.2)
ax.set_title("DCR exons less conserved \u2014 survives GC control",loc="left",fontsize=8)
ax.legend(frameon=False,fontsize=6.5,loc='lower right')
for i,(b,mm) in enumerate(zip(base,matched)):
    ax.text(i-w/2,b+0.03,f"{b:.2f}",ha='center',fontsize=6)
    ax.text(i+w/2,mm+0.03,f"{mm:.2f}",ha='center',fontsize=6)
# Panel B: LM coefficient forest
ax=axes[1]
terms=['is_DCR','GC (per %)','length (per bp)']
coefs=[m1.params['is_DCR'],m1.params['gc'],m1.params['len']]
ci=m1.conf_int()
errs=[[coefs[i]-ci.iloc[[1,2,3][i],0] for i in range(3)],
      [ci.iloc[[1,2,3][i],1]-coefs[i] for i in range(3)]]
yy=np.arange(3)[::-1]
ax.errorbar(coefs,yy,xerr=errs,fmt='o',color='#c0392b',capsize=3,ms=5)
ax.axvline(0,color='grey',lw=0.8,ls='--')
ax.set_yticks(yy[::-1][::-1]);ax.set_yticks(np.arange(3)[::-1]);ax.set_yticklabels(terms,fontsize=7)
ax.set_xlabel("phyloP effect (LM coefficient)")
ax.set_title("meanPP ~ is_DCR + GC + length\nis_DCR effect independent of GC (p=2e-10)",loc="left",fontsize=8)
fig.suptitle("R2-M4a: reduced DCR nucleotide conservation is not a GC artifact",fontsize=9,x=0.02,ha='left')
fig.tight_layout(rect=[0,0,1,0.95]);fig.savefig(f"{WS}/tierA/phyloP_gc_control.png",dpi=300)
print("saved");plt.close(fig)
~~~~
