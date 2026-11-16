# 11_annotation_depth_check — as executed, not refactored

**Annotation-depth (transcript count) check of DCR discovery and the X' gap.**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[451..481]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 11 (msg 451 … msg 471) |
| Language / environment | bash (env: dcr-tierA); python (env: dcr-tierA) |
| Saved outputs | `annotation_bias.png`, `annotation_depth_check.csv`, `dcr_discovery_vs_annotation.csv`, `xprime_lifttool_comparison.csv` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 5 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

---

### [msg 451] bash · env: dcr-tierA · Inspecting transcript annotation depth

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
X="$WS/tierA/xprime"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
# Transcripts per gene: does DCR status track annotation depth?
head -2 "$REF/correspondance_ENSG_ENST_v114.hg38.14"
echo "---"
awk '{print $1}' "$REF/correspondance_ENSG_ENST_v114.hg38.14" | sort | uniq -c | awk '{print $1}' | sort -n | uniq -c | tail -5
echo "=== total genes / transcripts ==="
awk '{print $1}' "$REF/correspondance_ENSG_ENST_v114.hg38.14" | sort -u | wc -l
wc -l < "$REF/correspondance_ENSG_ENST_v114.hg38.14"
~~~~

### [msg 453] python · env: dcr-tierA · Loading transcript counts per gene

~~~~python
REFD="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
corr=pd.read_csv(f"{REFD}/correspondance_ENSG_ENST_v114.hg38.14",sep="\t")
corr.columns=['ENSG','ENST']
ntx=corr.groupby('ENSG').ENST.nunique()
# map exons -> ENST -> ENSG -> n transcripts
enst2ensg=dict(zip(corr.ENST,corr.ENSG))
def enst_of(name): return name.split('_cds_')[0].split('.')[0]
# DCR genes vs non-DCR genes annotation depth
dcrbed=pd.read_csv(f"{D}/output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed",
                   sep="\t",header=None,dtype=str)
print("DCR bed cols:",dcrbed.shape[1]); print(dcrbed.head(2).to_string())
~~~~

### [msg 455] python · env: dcr-tierA · Inspecting DCR bed gene column

~~~~python
DOO="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
dcrbed=pd.read_csv(f"{DOO}/output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed",
                   sep="\t",header=None,dtype=str)
print("cols:",dcrbed.shape[1])
print(dcrbed.head(3).to_string())
~~~~

### [msg 457] python · env: dcr-tierA · Inspecting GENCODE bed name field

~~~~python
dcr_genes=set(dcrbed[5].str.upper())
print("DCR genes:",len(dcr_genes))
# gene symbol -> need ENSG. Use the GENCODE bed to map symbol->ENST->ENSG
gc=pd.read_csv(f"{REFD}/hg38.GENECODE.V48.bed",sep="\t",header=None,dtype=str,usecols=[3],names=['n'])
print(gc.head(3).to_string())
~~~~

### [msg 459] python · env: dcr-tierA · Comparing annotation depth DCR vs non-DCR

~~~~python
Nsub2['enst']=Nsub2.name.str.split('_cds_').str[0].str.split('.').str[0]
Nsub2['ensg']=Nsub2.enst.map(enst2ensg)
Nsub2['ntx']=Nsub2.ensg.map(ntx)
print("N intervals with transcript-count annotation: %d/%d"%(Nsub2.ntx.notna().sum(),len(Nsub2)))
# DCR exon side
dexon=pd.read_csv(f"{DOO}/output/02_Conservation/DCR_exons_DCR_genes_UNIQ.bed",sep="\t",header=None,dtype=str)
dexon['enst']=dexon[3].str.split('_cds_').str[0].str.split('.').str[0]
dexon['ensg']=dexon.enst.map(enst2ensg); dexon['ntx']=dexon.ensg.map(ntx)
print("\n=== ANNOTATION DEPTH: transcripts per gene ===")
print(f"DCR-containing exons' genes  : median {dexon.ntx.median():.0f} transcripts (n={dexon.ntx.notna().sum()})")
print(f"non-DCR exons' genes         : median {Nsub2.ntx.median():.0f} transcripts (n={Nsub2.ntx.notna().sum()})")
from scipy.stats import mannwhitneyu
u,pv=mannwhitneyu(dexon.ntx.dropna(),Nsub2.ntx.dropna())
n1,n2=dexon.ntx.notna().sum(),Nsub2.ntx.notna().sum()
print(f"MWU p={pv:.3e}, rank-biserial r={1-2*u/(n1*n2):.3f}")
~~~~

### [msg 461] python · env: dcr-tierA · Testing DCR discovery vs annotation depth

~~~~python
# Direct annotation-bias test: does DCR discovery scale with transcript count?
allg=pd.DataFrame({'ntx':ntx})
dcr_ensg=set(dexon.ensg.dropna())
allg['isDCR']=allg.index.isin(dcr_ensg)
allg=allg[allg.ntx>=1]
q=pd.qcut(allg.ntx,8,duplicates='drop')
disc=allg.groupby(q,observed=True).agg(ngene=('isDCR','size'),pDCR=('isDCR','mean'),medtx=('ntx','median'))
disc['pDCR']*=100
print("=== DCR discovery rate vs annotation depth (all genes) ===")
print(disc.to_string(float_format=lambda v:f"{v:.2f}"))
r_ann,p_ann=stats.spearmanr(allg.ntx,allg.isDCR.astype(int))
print(f"\nSpearman(n_transcripts, is_DCR_gene) = {r_ann:.3f}, p={p_ann:.2e}")
~~~~

### [msg 463] python · env: dcr-tierA · Stratifying gap by annotation depth

~~~~python
# Joint length + annotation-depth matched comparison, UCSC lift
Nu['nid2']=Nu.nid
Nu2=Nu.merge(Nsub2[['ntx','ensg']].reset_index().rename(columns={'index':'nid2'}),on='nid2',how='left')
Du['enst']=[hD_map and h.split('::')[0] for h in Du.h]
# D side: attach ntx via the DCR exon table (gene-level), using dceID->gene from dcrbed
dcrbed.columns=['chr','start','end','dceID','strand','gene','nfr']
sym2ensg={}
for e,g in zip(dexon.ensg,dexon[3]): pass
# map DCR region -> gene symbol -> ensg via dexon gene set is unreliable; use ntx of DCR genes directly
gene_ntx=dexon.groupby('ensg').ntx.first()
# Instead: annotate D by the ntx of its gene through overlap with DCR exons
print("D intervals:",len(Du),"  N intervals:",len(Nu2))
# Stratify N by annotation depth and compute length-matched N' within each stratum
Nu2['txbin']=pd.cut(Nu2.ntx,[0,2,4,10,10**4],labels=['1-2','3-4','5-10','>10'])
out=[]
for tb,g in Nu2.groupby('txbin',observed=True):
    gg=g.copy(); gg['lb']=pd.cut(gg.hlen,bins)
    rn_=gg.groupby('lb',observed=True).strict.mean()*100
    com=[b for b in rd.index if b in rn_.index and wd[b]>=20]
    if not com: continue
    wv=wd[com]/wd[com].sum()
    dm_=float((rd.reindex(com)*wv).sum()); nm_=float((rn_.reindex(com)*wv).sum())
    out.append((tb,len(gg),dm_,nm_,dm_-nm_))
ann=pd.DataFrame(out,columns=['transcripts_per_gene','n','D_prime','N_prime','gap'])
print("\n=== D'-N' gap stratified by annotation depth of the N gene ===")
print(ann.to_string(index=False,float_format=lambda v:f"{v:.1f}"))
~~~~

### [msg 465] python · env: dcr-tierA · Quantifying ascertainment undercount

~~~~python
# If discovery rate in the best-annotated stratum applied genome-wide, how many DCR genes would there be?
top=disc.iloc[-1]; obs_dcr_genes=allg.isDCR.sum()
proj=allg.shape[0]*top.pDCR/100
print(f"observed DCR genes                     : {obs_dcr_genes:,}")
print(f"projected if all genes annotated as deeply as the top stratum: {proj:,.0f}")
print(f"implied undercount factor              : {proj/obs_dcr_genes:.1f}x")
print()
# restrict to genes that COULD yield a DCR (>=2 transcripts) - the fair denominator
elig=allg[allg.ntx>=2]
print(f"genes with >=2 transcripts (DCR-eligible): {len(elig):,} of {len(allg):,} ({len(elig)/len(allg)*100:.1f}%)")
print(f"DCR rate among eligible genes           : {elig.isDCR.mean()*100:.2f}%")
print(f"DCR rate among >10-transcript genes     : {top.pDCR:.2f}%")
# CMH on UCSC data with annotation strata
tabs2=[]
for tb,g in Nu2.groupby('txbin',observed=True):
    for b in cmp.index:
        d=Du[Du.lb==b]; n=g[g.hlen.between(b.left,b.right)]
        if len(d)<20 or len(n)<20: continue
        tabs2.append([[d.strict.sum(),len(d)-d.strict.sum()],[n.strict.sum(),len(n)-n.strict.sum()]])
t2=np.array(tabs2,dtype=float)
a2,b2,c2,d2=t2[:,0,0],t2[:,0,1],t2[:,1,0],t2[:,1,1]; nT2=t2.sum(axis=(1,2))
or2=(a2*d2/nT2).sum()/(b2*c2/nT2).sum()
E2=((a2+b2)*(a2+c2)/nT2).sum(); V2=((a2+b2)*(c2+d2)*(a2+c2)*(b2+d2)/(nT2**2*(nT2-1))).sum()
chi_2=(abs(a2.sum()-E2)-0.5)**2/V2
print(f"\nCMH stratified by BOTH length and annotation depth: OR={or2:.3f}, p={1-chi2.cdf(chi_2,1):.3e}, strata={len(t2)}")
ann.to_csv(f"{X}/annotation_depth_check.csv",index=False)
disc.to_csv(f"{X}/dcr_discovery_vs_annotation.csv")
~~~~

### [msg 467] python · env: dcr-tierA · Building annotation-bias figure

~~~~python
fig,axes=plt.subplots(1,3,figsize=(11,3.6))
ax=axes[0]
xx=np.arange(len(disc))
ax.bar(xx,disc.pDCR,color=ORANGE,edgecolor='k',linewidth=0.5,width=0.62)
for i,v in enumerate(disc.pDCR): ax.text(i,v+0.12,f"{v:.2f}%",ha='center',fontsize=7.5,weight='bold')
ax.set_xticks(xx); ax.set_xticklabels(['1-2','3-4','5-10','>10'],fontsize=7.5)
ax.set_xlabel("annotated transcripts per gene"); ax.set_ylabel("% of genes containing a DCR")
ax.set_ylim(0,7.2)
ax.text(0.03,0.88,f"Spearman r={r_ann:.2f}\np < 1e-300",transform=ax.transAxes,fontsize=7)
ax.set_title("A  DCR discovery is annotation-limited",loc='left',weight='bold')

ax=axes[1]
ax.bar([0,1],[obs_dcr_genes,proj],color=[GREY,BLUE],edgecolor='k',linewidth=0.5,width=0.5)
for i,v in enumerate([obs_dcr_genes,proj]): ax.text(i,v+90,f"{v:,.0f}",ha='center',fontsize=8.5,weight='bold')
ax.set_xticks([0,1]); ax.set_xticklabels(["observed","projected at\nuniform deep\nannotation"],fontsize=7.5)
ax.set_ylabel("genes containing a DCR"); ax.set_ylim(0,6100)
ax.text(0.5,0.55,f"{proj/obs_dcr_genes:.1f}× undercount",transform=ax.transAxes,ha='center',
        fontsize=10,weight='bold',style='italic',color=BLUE)
ax.set_title("B  Current DCR catalogue is a floor",loc='left',weight='bold')

ax=axes[2]
xx=np.arange(len(ann))
ax.bar(xx,ann.gap,color=BLUE,edgecolor='k',linewidth=0.5,width=0.6)
for i,v in enumerate(ann.gap): ax.text(i,v+0.16,f"{v:.1f}",ha='center',fontsize=7.5,weight='bold')
ax.set_xticks(xx); ax.set_xticklabels(ann.transcripts_per_gene,fontsize=7.5)
ax.set_xlabel("annotated transcripts per gene (N set)")
ax.set_ylabel("D' − N' gap (pts, length-matched)")
ax.set_ylim(0,8.2)
ax.axhline(4.5,ls='--',lw=0.9,color=GREY)
ax.text(3.3,4.75,"overall 4.5",fontsize=6.8,color=GREY,ha='right')
ax.set_title("C  ...but it doesn't rescue X'",loc='left',weight='bold')
fig.suptitle("Annotation depth drives DCR ascertainment (4× undercount) but does not explain the weak conservation signal",
             x=0.005,ha='left',fontsize=9.8,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/annotation_bias.png",dpi=300); plt.close(fig)
print("saved")
~~~~

### [msg 469] python · env: dcr-tierA · Fixing annotation figure label

~~~~python
ax=axes[1]
fig,axes=plt.subplots(1,3,figsize=(11,3.6))
ax=axes[0]
xx=np.arange(len(disc))
ax.bar(xx,disc.pDCR,color=ORANGE,edgecolor='k',linewidth=0.5,width=0.62)
for i,v in enumerate(disc.pDCR): ax.text(i,v+0.12,f"{v:.2f}%",ha='center',fontsize=7.5,weight='bold')
ax.set_xticks(xx); ax.set_xticklabels(['1-2','3-4','5-10','>10'],fontsize=7.5)
ax.set_xlabel("annotated transcripts per gene"); ax.set_ylabel("% of genes containing a DCR")
ax.set_ylim(0,7.2)
ax.text(0.03,0.88,f"Spearman r={r_ann:.2f}\np < 1e-300",transform=ax.transAxes,fontsize=7)
ax.set_title("A  DCR discovery is annotation-limited",loc='left',weight='bold')

ax=axes[1]
ax.bar([0,1],[obs_dcr_genes,proj],color=[GREY,BLUE],edgecolor='k',linewidth=0.5,width=0.5)
for i,v in enumerate([obs_dcr_genes,proj]): ax.text(i,v+110,f"{v:,.0f}",ha='center',fontsize=8.5,weight='bold')
ax.set_xticks([0,1]); ax.set_xticklabels(["observed","projected at\nuniform deep\nannotation"],fontsize=7.5)
ax.set_ylabel("genes containing a DCR"); ax.set_ylim(0,6600)
ax.annotate("",xy=(0.72,obs_dcr_genes),xytext=(0.72,proj),arrowprops=dict(arrowstyle='<->',lw=1.2,color=BLUE))
ax.text(0.60,(obs_dcr_genes+proj)/2,f"{proj/obs_dcr_genes:.1f}×\nundercount",fontsize=8.5,weight='bold',
        style='italic',color=BLUE,ha='right',va='center')
ax.set_title("B  Current DCR catalogue is a floor",loc='left',weight='bold')

ax=axes[2]
xx=np.arange(len(ann))
ax.bar(xx,ann.gap,color=BLUE,edgecolor='k',linewidth=0.5,width=0.6)
for i,v in enumerate(ann.gap): ax.text(i,v+0.16,f"{v:.1f}",ha='center',fontsize=7.5,weight='bold')
ax.set_xticks(xx); ax.set_xticklabels(ann.transcripts_per_gene,fontsize=7.5)
ax.set_xlabel("annotated transcripts per gene (N set)")
ax.set_ylabel("D' − N' gap (pts, length-matched)")
ax.set_ylim(0,8.2)
ax.axhline(4.5,ls='--',lw=0.9,color=GREY)
ax.text(3.42,4.75,"overall 4.5",fontsize=6.8,color=GREY,ha='right')
ax.set_title("C  ...but it doesn't rescue X'",loc='left',weight='bold')
fig.suptitle("Annotation depth drives DCR ascertainment (4× undercount) but does not explain the weak conservation signal",
             x=0.005,ha='left',fontsize=9.8,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/annotation_bias.png",dpi=300); plt.close(fig)
print("saved")
~~~~

### [msg 471] python · env: dcr-tierA · Saving lift tool comparison

~~~~python
ucsc=pd.DataFrame([{'tool':'UCSC liftOver','D_lifted':len(Du),'N_lifted':len(Nu),
  'D_prime_raw':Du.strict.mean()*100,'N_prime_raw':Nu.strict.mean()*100,
  'D_prime_matched':Du_m,'N_prime_matched':Nu_m,'gap':Du_m-Nu_m},
 {'tool':'pyliftover','D_lifted':1629,'N_lifted':32691,'D_prime_raw':81.5,'N_prime_raw':69.9,
  'D_prime_matched':81.8,'N_prime_matched':77.5,'gap':4.3}])
ucsc.to_csv(f"{WS}/tierA/xprime_lifttool_comparison.csv",index=False)
print(ucsc.to_string(index=False,float_format=lambda v:f"{v:.1f}"))
~~~~
