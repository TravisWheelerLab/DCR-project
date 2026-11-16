# 09_definition_reconciliation — as executed, not refactored

**Reconciliation of the 82% vs 48% dual-coding definitions (strict vs sub-interval rule).**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[344..392]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 12 (msg 348 … msg 387) |
| Language / environment | bash (env: dcr-tierA); bash (env: python); python (env: dcr-tierA) |
| Saved outputs | `plan_revised.json`, `definition_reconciliation.png`, `definition_reconciliation.csv` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 6 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

---

### [msg 348] bash · env: python · Checking NMD/TrEMBL/isoform data availability

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
echo "=== New_TableS1 header ==="
head -1 "$D/output/New_TableS1.tsv" | tr '\t' '\n' | head -40
echo ""
echo "=== output dir subfolders ==="
ls "$D/output/"
echo ""
echo "=== grep for NMD / canonical / trembl / tallyman / bath in outputs & scripts ==="
ls "$D/output/"*/ 2>/dev/null | grep -i -E "nmd|canon|trembl|tally|bath|junction|isoform" | head -30
echo "--- scripts ---"
find "$D/scripts" -type f 2>/dev/null | xargs grep -l -i -E "nmd|trembl|tallyman|bath" 2>/dev/null | head
~~~~

### [msg 367] bash · env: dcr-tierA · Extracting human DCR region sequences

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
echo "=== existing gc fastas ==="
ls -la "$WS/tierA/gc/" | head -20
echo ""
echo "=== DCR interval bed (1741) head ==="
head -3 "$D/output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed"
echo ""
echo "=== extract HUMAN DCR REGIONS from hg38 ==="
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"dce"$7"_"NR,0,$5}' "$D/output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed" > "$WS/tierA/gc/human_DCR_regions.bed"
wc -l < "$WS/tierA/gc/human_DCR_regions.bed"
head -2 "$WS/tierA/gc/human_DCR_regions.bed"
bedtools getfasta -fi "$REF/hg38.fa" -bed "$WS/tierA/gc/human_DCR_regions.bed" -s -fo "$WS/tierA/gc/human_DCR_regions.fa" 2>&1 | head -3
grep -c '^>' "$WS/tierA/gc/human_DCR_regions.fa"
head -2 "$WS/tierA/gc/human_DCR_regions.fa" | cut -c1-80
~~~~

### [msg 369] python · env: dcr-tierA · Running strict rule across three input sets

~~~~python
import sys, numpy as np, pandas as pd
WS="<OLD_WORKSPACE>"
GC=f"{WS}/tierA/gc"
sys.path.insert(0,f"{WS}/tierA")
from esl_translate_replica import esl_translate_watson

def read_fa(path):
    out={};hdr=None;cur=[]
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if hdr is not None: out[hdr]=''.join(cur).upper()
                hdr=line[1:].strip(); cur=[]
            else: cur.append(line.strip())
        if hdr is not None: out[hdr]=''.join(cur).upper()
    return out

def strict_rule(seq):
    """Pipeline rule: >1 ORF spanning the WHOLE input sequence (+/-3bp)."""
    L=len(seq)
    orfs=esl_translate_watson(seq,min_aa=1)
    full=[o for o in orfs if (o[0]-1-3)<=0 and (L-o[1])<=3]
    return len(full)>1, len(full), L

def subinterval_rule(seq, min_overlap=21):
    """Do >=2 ORFs in DIFFERENT frames overlap by >= min_overlap nt?
       Returns (bool, max_overlap). This is 'is any sub-interval dual-coding'."""
    orfs=esl_translate_watson(seq,min_aa=1)
    best=0
    for i in range(len(orfs)):
        s1,e1,f1,_=orfs[i]
        for j in range(i+1,len(orfs)):
            s2,e2,f2,_=orfs[j]
            if f1==f2: continue
            ov=min(e1,e2)-max(s1,s2)+1
            if ov>best: best=ov
    return best>=min_overlap, best

sets={
 "1_mouse_DCR_region":  f"{GC}/mouse_DCR_ortholog.fa",
 "2_human_DCR_region":  f"{GC}/human_DCR_regions.fa",
 "3_human_DCR_exon":    f"{GC}/DCR_exons_DCR_genes_UNIQ.fa",
}
res={}
for name,path in sets.items():
    seqs=read_fa(path)
    rows=[]
    for h,s in seqs.items():
        if len(s)<6 or set(s)-set("ACGTN"): 
            if set(s)-set("ACGTN"): pass
        strict,nfull,L=strict_rule(s)
        rows.append((h,L,strict,nfull))
    df=pd.DataFrame(rows,columns=['id','len','strict','nfull'])
    res[name]=df
    print(f"{name}: n={len(df)}, median_len={df.len.median():.0f}, strict_dualcoding={df.strict.mean()*100:.1f}%")
~~~~

### [msg 371] bash · env: dcr-tierA · Computing per-exon DCR overlap fraction

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
cd "$WS/tierA/gc"
# per-exon: total DCR bp overlapping it
sort -k1,1 -k2,2n DCR_exons_DCR_genes_UNIQ.b3.bed > ex.sorted.bed
sort -k1,1 -k2,2n human_DCR_regions.bed > dcr.sorted.bed
bedtools intersect -a ex.sorted.bed -b dcr.sorted.bed -wao | \
  awk 'BEGIN{OFS="\t"}{key=$1":"$2"-"$3; ov[key]+=$NF; len[key]=$3-$2} END{for(k in ov) print k,len[k],ov[k]}' > exon_dcr_overlap.tsv
wc -l < exon_dcr_overlap.tsv
head -3 exon_dcr_overlap.tsv
~~~~

### [msg 373] python · env: dcr-tierA · Testing DCR-coverage dependence of strict rule

~~~~python
ov=pd.read_csv(f"{GC}/exon_dcr_overlap.tsv",sep="\t",header=None,names=['key','exlen','dcrbp'])
ov['ratio']=ov.dcrbp/ov.exlen
ex=res['3_human_DCR_exon'].copy()
ex['key']=ex.id.str.replace(r'\(.\)$','',regex=True)
m=ex.merge(ov,on='key',how='inner')
print("joined:",len(m),"of",len(ex))
# ratio dependence: strict-rule pass rate by DCR/exon coverage decile
m['dec']=pd.qcut(m.ratio,10,labels=False,duplicates='drop')
tab=m.groupby('dec').agg(n=('strict','size'),cov=('ratio','median'),pass_rate=('strict','mean')).reset_index()
print("\nStrict-rule pass rate by DCR-coverage decile:")
print(tab.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
from scipy import stats
r,p=stats.pearsonr(m.ratio,m.strict.astype(int))
print(f"\ncorr(DCR coverage fraction, strict pass) r={r:.3f} p={p:.3e}")
print("pass rate when exon is ~entirely DCR (ratio>=0.95): %.1f%% (n=%d)"%(m[m.ratio>=0.95].strict.mean()*100,(m.ratio>=0.95).sum()))
print("pass rate when DCR is small part of exon (ratio<0.5): %.1f%% (n=%d)"%(m[m.ratio<0.5].strict.mean()*100,(m.ratio<0.5).sum()))
~~~~

### [msg 375] python · env: dcr-tierA · Running sub-interval scan on DCR exons

~~~~python
import time
# Sub-interval scan on human DCR exons: do >=2 ORFs in different frames overlap by >=21nt?
seqs=read_fa(f"{GC}/DCR_exons_DCR_genes_UNIQ.fa")
t0=time.time(); rows=[]
for h,s in seqs.items():
    ok,best=subinterval_rule(s,min_overlap=21)
    rows.append((h,len(s),ok,best))
sub=pd.DataFrame(rows,columns=['id','len','sub','maxov'])
print(f"human DCR exon, sub-interval rule: {sub.sub.mean()*100:.1f}%  (n={len(sub)}, {time.time()-t0:.0f}s)")
sub['key']=sub.id.str.replace(r'\(.\)$','',regex=True)
m2=sub.merge(ov,on='key',how='inner')
r2,p2=stats.pearsonr(m2.ratio,m2['sub'].astype(int))
print(f"corr(DCR coverage, sub-interval pass) r={r2:.3f} p={p2:.3e}")
print("  ratio<0.5: %.1f%%   ratio>=0.95: %.1f%%"%(m2[m2.ratio<0.5]['sub'].mean()*100, m2[m2.ratio>=0.95]['sub'].mean()*100))
~~~~

### [msg 377] python · env: dcr-tierA · Comparing sub-interval vs strict rule

~~~~python
sub=sub.rename(columns={'sub':'dual_sub'})
print(f"human DCR exon, sub-interval rule: {sub.dual_sub.mean()*100:.1f}%  (n={len(sub)})")
m2=sub.merge(ov,on='key',how='inner')
r2,p2=stats.pearsonr(m2.ratio,m2.dual_sub.astype(int))
print(f"corr(DCR coverage, sub-interval pass) r={r2:.3f} p={p2:.3e}")
print("  ratio<0.5: %.1f%%   ratio>=0.95: %.1f%%"%(m2[m2.ratio<0.5].dual_sub.mean()*100, m2[m2.ratio>=0.95].dual_sub.mean()*100))
print("\nSTRICT rule for comparison:  overall %.1f%% | ratio<0.5 %.1f%% | ratio>=0.95 %.1f%%"%(
    m.strict.mean()*100, m[m.ratio<0.5].strict.mean()*100, m[m.ratio>=0.95].strict.mean()*100))
~~~~

### [msg 379] python · env: dcr-tierA · Completing four-way rule comparison

~~~~python
sub['key']=sub.id.str.replace(r'\(.\)$','',regex=True)
m2=sub.merge(ov,on='key',how='inner')
r2,p2=stats.pearsonr(m2.ratio,m2.dual_sub.astype(int))
print(f"SUB-INTERVAL rule: overall {sub.dual_sub.mean()*100:.1f}% | ratio<0.5 {m2[m2.ratio<0.5].dual_sub.mean()*100:.1f}% | ratio>=0.95 {m2[m2.ratio>=0.95].dual_sub.mean()*100:.1f}%")
print(f"   corr with DCR coverage r={r2:.3f} p={p2:.3e}   <- should be ~0 if geometry-independent")
print(f"STRICT rule:       overall {m.strict.mean()*100:.1f}% | ratio<0.5 {m[m.ratio<0.5].strict.mean()*100:.1f}% | ratio>=0.95 {m[m.ratio>=0.95].strict.mean()*100:.1f}%")
print(f"   corr with DCR coverage r={r:.3f} p={p:.3e}")
# also run sub-interval on mouse + human regions for the 4-way table
for nm,path in [("1_mouse_DCR_region",f"{GC}/mouse_DCR_ortholog.fa"),("2_human_DCR_region",f"{GC}/human_DCR_regions.fa")]:
    sq=read_fa(path); v=[subinterval_rule(s,21)[0] for s in sq.values()]
    res[nm]['dual_sub']=v
    print(f"{nm}: sub-interval {np.mean(v)*100:.1f}%")
~~~~

### [msg 381] python · env: dcr-tierA · Calibrating dual-coding overlap threshold

~~~~python
hr=res['2_human_DCR_region']
print("Human DCR regions: strict %.1f%%, sub %.1f%%"%(hr.strict.mean()*100,hr.dual_sub.mean()*100))
print("regions <21nt: n=%d"%(hr.len<21).sum())
short=hr[hr.len<24]
print("of regions <24nt: strict %.1f%%, sub %.1f%%  (n=%d)"%(short.strict.mean()*100,short.dual_sub.mean()*100,len(short)))
long_=hr[hr.len>=24]
print("of regions >=24nt: strict %.1f%%, sub %.1f%%  (n=%d)"%(long_.strict.mean()*100,long_.dual_sub.mean()*100,len(long_)))

# Threshold sweep: DCR exons vs matched non-DCR exons. Sample non-DCR for speed.
import random
rng=random.Random(7)
nd=read_fa(f"{GC}/nonDCR_exon_DCR_gene_UNIQ.fa")
nd_items=list(nd.items()); rng.shuffle(nd_items); nd_s=dict(nd_items[:4000])
dcr_ex=read_fa(f"{GC}/DCR_exons_DCR_genes_UNIQ.fa")
def maxov(seq):
    orfs=esl_translate_watson(seq,min_aa=1); best=0
    for i in range(len(orfs)):
        s1,e1,f1,_=orfs[i]
        for j in range(i+1,len(orfs)):
            s2,e2,f2,_=orfs[j]
            if f1==f2: continue
            o=min(e1,e2)-max(s1,s2)+1
            if o>best: best=o
    return best
ov_dcr=np.array([maxov(s) for s in dcr_ex.values()])
ov_nd =np.array([maxov(s) for s in nd_s.values()])
len_dcr=np.array([len(s) for s in dcr_ex.values()]); len_nd=np.array([len(s) for s in nd_s.values()])
print("\nmedian max-overlap: DCR exon %d nt | non-DCR exon %d nt"%(np.median(ov_dcr),np.median(ov_nd)))
print("\nthreshold sweep (absolute nt overlap):")
print(" thr   DCRexon%  nonDCR%   gap")
for thr in [21,30,45,60,90,120,150]:
    a=(ov_dcr>=thr).mean()*100; b=(ov_nd>=thr).mean()*100
    print(f" {thr:4d}   {a:6.1f}   {b:6.1f}   {a-b:5.1f}")
print("\nthreshold sweep (overlap as FRACTION of exon length):")
fd=ov_dcr/len_dcr; fn=ov_nd/len_nd
print(" frac  DCRexon%  nonDCR%   gap")
for thr in [0.2,0.3,0.5,0.7,0.9]:
    a=(fd>=thr).mean()*100; b=(fn>=thr).mean()*100
    print(f" {thr:4.1f}   {a:6.1f}   {b:6.1f}   {a-b:5.1f}")
~~~~

### [msg 383] python · env: dcr-tierA · Rendering reconciliation figure and table

~~~~python
import matplotlib as mpl, matplotlib.pyplot as plt
mpl.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7.5,
                     'ytick.labelsize':7.5,'legend.fontsize':7,'axes.spines.top':False,'axes.spines.right':False,
                     'figure.dpi':300,'savefig.bbox':'tight'})
ORANGE='#FFA800'; GREY='#9e9e9e'; BLUE='#2c7fb8'

fig,axes=plt.subplots(1,3,figsize=(10.2,3.3))
# Panel A: 4-way
ax=axes[0]
labs=["mouse DCR\nregion","human DCR\nregion","human DCR\nexon","human DCR exon\n(sub-interval)"]
vals=[res['1_mouse_DCR_region'].strict.mean()*100, res['2_human_DCR_region'].strict.mean()*100,
      res['3_human_DCR_exon'].strict.mean()*100, sub.dual_sub.mean()*100]
cols=[ORANGE,ORANGE,GREY,BLUE]
b=ax.bar(range(4),vals,color=cols,edgecolor='k',linewidth=0.5)
for i,v in enumerate(vals): ax.text(i,v+2,f"{v:.1f}%",ha='center',fontsize=7.5,weight='bold')
ax.set_xticks(range(4)); ax.set_xticklabels(labs,fontsize=6.8)
ax.set_ylabel("% scored dual-coding"); ax.set_ylim(0,112)
ax.set_title("A  Same rule, different objects",loc='left',weight='bold')
ax.text(0.5,-0.34,"full-length ORF rule            permissive",transform=ax.transAxes,ha='center',fontsize=6.5,style='italic')

# Panel B: geometry effect
ax=axes[1]
ax.plot(tab.cov*100,tab.pass_rate*100,'o-',color=GREY,lw=1.5,ms=5,label='strict full-length rule')
sub_tab=m2.assign(dec=pd.qcut(m2.ratio,10,labels=False,duplicates='drop')).groupby('dec').agg(
    cov=('ratio','median'),pr=('dual_sub','mean')).reset_index()
ax.plot(sub_tab.cov*100,sub_tab.pr*100,'s-',color=BLUE,lw=1.5,ms=4,label='sub-interval rule')
ax.set_xlabel("DCR coverage of the exon (%)"); ax.set_ylabel("% scored dual-coding")
ax.set_ylim(0,105); ax.legend(loc='center right',frameon=False)
ax.set_title("B  The 48% is exon geometry",loc='left',weight='bold')
ax.text(0.03,0.06,f"strict: r={r:.2f}, p={p:.0e}\nsub-interval: r={r2:.2f}",transform=ax.transAxes,fontsize=6.5)

# Panel C: discrimination sweep
ax=axes[2]
fr=[0.2,0.3,0.5,0.7,0.9,1.0]
a=[(fd>=t).mean()*100 for t in fr[:-1]]+[m.strict.mean()*100]
bb=[(fn>=t).mean()*100 for t in fr[:-1]]+[41.0]
ax.plot(fr,a,'o-',color=ORANGE,lw=1.6,ms=5,label='DCR exon')
ax.plot(fr,bb,'o-',color=GREY,lw=1.6,ms=5,label='non-DCR exon')
ax.set_xlabel("required dual-coding span\n(fraction of exon length)"); ax.set_ylabel("% scored dual-coding")
ax.set_ylim(0,105); ax.legend(frameon=False,loc='lower left')
ax.set_title("C  Potential barely discriminates",loc='left',weight='bold')
ax.annotate("gap only ~5-11 pts\nat any threshold",xy=(0.62,0.72),xycoords='axes fraction',fontsize=6.8,style='italic')
fig.suptitle("Reconciling the dual-coding measurements: 82% (region) vs 48% (containing exon)",
             x=0.005,ha='left',fontsize=10,weight='bold')
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(f"{WS}/tierA/definition_reconciliation.png",dpi=300); plt.close(fig)

recon=pd.DataFrame({
 'measurement':['mouse DCR region (strict)','human DCR region (strict)','human DCR exon (strict)',
                'human DCR exon (sub-interval >=21nt)','human DCR exon, DCR covers <50%','human DCR exon, DCR covers >=95%',
                'non-DCR exon (strict, Tier A)','median max dual-coding overlap, DCR exon','median max dual-coding overlap, non-DCR exon'],
 'n':[1599,1739,3324,3324,1016,1108,4000,3324,4000],
 'value':[f"{res['1_mouse_DCR_region'].strict.mean()*100:.1f}%",f"{res['2_human_DCR_region'].strict.mean()*100:.1f}%",
          f"{res['3_human_DCR_exon'].strict.mean()*100:.1f}%",f"{sub.dual_sub.mean()*100:.1f}%",
          f"{m[m.ratio<0.5].strict.mean()*100:.1f}%",f"{m[m.ratio>=0.95].strict.mean()*100:.1f}%","~41%",
          f"{np.median(ov_dcr):.0f} nt",f"{np.median(ov_nd):.0f} nt"],
 'interpretation':['conservation of dual-coding in mouse orthologs (PUBLISHED)',
                   'definitional sanity check - validates the rule',
                   'geometry-confounded; NOT "fraction of DCR exons that are dual-coding"',
                   'permissive rule; near-saturated, non-discriminative',
                   'strict rule fails when DCR is a small part of the exon',
                   'strict rule mostly passes when exon is nearly all DCR',
                   'Tier A matched comparator',
                   'dual-coding potential is near-universal in coding exons',
                   'essentially identical to DCR exons']})
recon.to_csv(f"{WS}/tierA/definition_reconciliation.csv",index=False)
print(recon.to_string(index=False))
~~~~

### [msg 385] python · env: dcr-tierA · Rendering reconciliation figure and table

~~~~python
tab=tab.rename(columns={'cov':'covm'})
sub_tab=m2.assign(dec=pd.qcut(m2.ratio,10,labels=False,duplicates='drop')).groupby('dec').agg(
    covm=('ratio','median'),pr=('dual_sub','mean')).reset_index()

fig,axes=plt.subplots(1,3,figsize=(10.2,3.3))
ax=axes[0]
labs=["mouse DCR\nregion","human DCR\nregion","human DCR\nexon","human DCR exon\n(sub-interval)"]
vals=[res['1_mouse_DCR_region'].strict.mean()*100, res['2_human_DCR_region'].strict.mean()*100,
      res['3_human_DCR_exon'].strict.mean()*100, sub.dual_sub.mean()*100]
ax.bar(range(4),vals,color=[ORANGE,ORANGE,GREY,BLUE],edgecolor='k',linewidth=0.5)
for i,v in enumerate(vals): ax.text(i,v+2,f"{v:.1f}%",ha='center',fontsize=7.5,weight='bold')
ax.set_xticks(range(4)); ax.set_xticklabels(labs,fontsize=6.8)
ax.set_ylabel("% scored dual-coding"); ax.set_ylim(0,112)
ax.set_title("A  Same rule, different objects",loc='left',weight='bold')

ax=axes[1]
ax.plot(tab.covm*100,tab.pass_rate*100,'o-',color=GREY,lw=1.5,ms=5,label='strict full-length rule')
ax.plot(sub_tab.covm*100,sub_tab.pr*100,'s-',color=BLUE,lw=1.5,ms=4,label='sub-interval rule')
ax.set_xlabel("DCR coverage of the exon (%)"); ax.set_ylabel("% scored dual-coding")
ax.set_ylim(0,105); ax.legend(loc='center right',frameon=False)
ax.set_title("B  The 48% is exon geometry",loc='left',weight='bold')
ax.text(0.03,0.06,f"strict: r={r:.2f}, p={p:.0e}\nsub-interval: r={r2:.2f}",transform=ax.transAxes,fontsize=6.5)

ax=axes[2]
fr=[0.2,0.3,0.5,0.7,0.9,1.0]
a=[(fd>=t).mean()*100 for t in fr[:-1]]+[m.strict.mean()*100]
bb=[(fn>=t).mean()*100 for t in fr[:-1]]+[41.0]
ax.plot(fr,a,'o-',color=ORANGE,lw=1.6,ms=5,label='DCR exon')
ax.plot(fr,bb,'o-',color=GREY,lw=1.6,ms=5,label='non-DCR exon')
ax.set_xlabel("required dual-coding span\n(fraction of exon length)"); ax.set_ylabel("% scored dual-coding")
ax.set_ylim(0,105); ax.legend(frameon=False,loc='lower left')
ax.set_title("C  Potential barely discriminates",loc='left',weight='bold')
ax.annotate("gap only ~5-11 pts\nat any threshold",xy=(0.40,0.70),xycoords='axes fraction',fontsize=6.8,style='italic')
fig.suptitle("Reconciling the dual-coding measurements: 82% (region) vs 48% (containing exon)",
             x=0.005,ha='left',fontsize=10,weight='bold')
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(f"{WS}/tierA/definition_reconciliation.png",dpi=300); plt.close(fig)
print("figure saved")

recon=pd.DataFrame({
 'measurement':['mouse DCR region (strict)','human DCR region (strict)','human DCR exon (strict)',
                'human DCR exon (sub-interval >=21nt)','human DCR exon, DCR covers <50% of exon','human DCR exon, DCR covers >=95% of exon',
                'non-DCR exon (strict, Tier A)','median max dual-coding overlap, DCR exon','median max dual-coding overlap, non-DCR exon'],
 'n':[1599,1739,3324,3324,1016,1108,4000,3324,4000],
 'value':[f"{res['1_mouse_DCR_region'].strict.mean()*100:.1f}%",f"{res['2_human_DCR_region'].strict.mean()*100:.1f}%",
          f"{res['3_human_DCR_exon'].strict.mean()*100:.1f}%",f"{sub.dual_sub.mean()*100:.1f}%",
          f"{m[m.ratio<0.5].strict.mean()*100:.1f}%",f"{m[m.ratio>=0.95].strict.mean()*100:.1f}%","~41%",
          f"{np.median(ov_dcr):.0f} nt",f"{np.median(ov_nd):.0f} nt"],
 'interpretation':['conservation of dual-coding in mouse orthologs (PUBLISHED, correct)',
                   'definitional sanity check - validates the rule',
                   'geometry-confounded; NOT "fraction of DCR exons that are dual-coding"',
                   'permissive rule; near-saturated, non-discriminative',
                   'strict rule fails when DCR is a small part of the exon',
                   'strict rule mostly passes when exon is nearly all DCR',
                   'Tier A matched comparator',
                   'dual-coding potential is near-universal in coding exons',
                   'essentially identical to DCR exons']})
recon.to_csv(f"{WS}/tierA/definition_reconciliation.csv",index=False)
print(recon.to_string(index=False))
~~~~

### [msg 387] python · env: dcr-tierA · Fixing panel A label spacing

~~~~python
fig,axes=plt.subplots(1,3,figsize=(10.6,3.5))
ax=axes[0]
labs=["mouse\nDCR region","human\nDCR region","human\nDCR exon","human DCR exon\nsub-interval"]
ax.bar(range(4),vals,color=[ORANGE,ORANGE,GREY,BLUE],edgecolor='k',linewidth=0.5,width=0.68)
for i,v in enumerate(vals): ax.text(i,v+2.5,f"{v:.1f}%",ha='center',fontsize=7.5,weight='bold')
ax.set_xticks(range(4)); ax.set_xticklabels(labs,fontsize=6.5)
ax.set_ylabel("% scored dual-coding"); ax.set_ylim(0,115); ax.set_xlim(-0.6,3.6)
ax.set_title("A  Same rule, different objects",loc='left',weight='bold')

ax=axes[1]
ax.plot(tab.covm*100,tab.pass_rate*100,'o-',color=GREY,lw=1.5,ms=5,label='strict full-length rule')
ax.plot(sub_tab.covm*100,sub_tab.pr*100,'s-',color=BLUE,lw=1.5,ms=4,label='sub-interval rule')
ax.set_xlabel("DCR coverage of the exon (%)"); ax.set_ylabel("% scored dual-coding")
ax.set_ylim(0,108); ax.set_xlim(0,120); ax.legend(loc='center right',frameon=False)
ax.set_title("B  The 48% is exon geometry",loc='left',weight='bold')
ax.text(0.03,0.05,f"strict: r={r:.2f}, p={p:.0e}\nsub-interval: r={r2:.2f}",transform=ax.transAxes,fontsize=6.5)

ax=axes[2]
ax.plot(fr,a,'o-',color=ORANGE,lw=1.6,ms=5,label='DCR exon')
ax.plot(fr,bb,'o-',color=GREY,lw=1.6,ms=5,label='non-DCR exon')
ax.set_xlabel("required dual-coding span\n(fraction of exon length)"); ax.set_ylabel("% scored dual-coding")
ax.set_ylim(0,108); ax.legend(frameon=False,loc='lower left')
ax.set_title("C  Potential barely discriminates",loc='left',weight='bold')
ax.annotate("gap only 5-11 pts\nat any threshold",xy=(0.34,0.62),xycoords='axes fraction',fontsize=6.8,style='italic')
fig.suptitle("Reconciling the dual-coding measurements: 82% (region) vs 48% (containing exon)",
             x=0.005,ha='left',fontsize=10,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/definition_reconciliation.png",dpi=300); plt.close(fig)
print("saved")
~~~~
