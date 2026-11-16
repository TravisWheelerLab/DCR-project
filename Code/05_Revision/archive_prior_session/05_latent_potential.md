# 05_latent_potential — as executed, not refactored

**Latent dual-coding potential in GC+length-matched non-DCR human exons.**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[267..278]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 4 (msg 267 … msg 273) |
| Language / environment | python (env: dcr-tierA) |
| Saved outputs | `latent_potential_nonDCR.png`, `latent_potential_nonDCR.csv` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 0 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

---

### [msg 267] python · env: dcr-tierA · Loading human exon sequences

~~~~python
# Human DCR exons (observed) and non-DCR human exons. Need sequences + the dual-coding rule applied
# to HUMAN exon sequences directly. Use hg38-extracted FASTAs already made.
def read_fasta_list(path):
    """Return list of (header, seq)."""
    out=[]; hdr=None; cur=[]
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if hdr is not None: out.append((hdr,''.join(cur).upper()))
                hdr=line[1:].strip(); cur=[]
            else: cur.append(line.strip())
        if hdr is not None: out.append((hdr,''.join(cur).upper()))
    return out

dcr_h  = read_fasta_list(f"{GC}/DCR_exons_DCR_genes_UNIQ.fa")
nDCR_dg= read_fasta_list(f"{GC}/nonDCR_exon_DCR_gene_UNIQ.fa")
nDCR_ng= read_fasta_list(f"{GC}/nonDCR_exon_nonDCR_gene_UNIQ.fa")
print("DCR human exons:",len(dcr_h),"| nonDCR in DCR genes:",len(nDCR_dg),"| nonDCR in nonDCR genes:",len(nDCR_ng))

def gcval(s):
    g=s.count('G')+s.count('C');n=g+s.count('A')+s.count('T');return 100.0*g/n if n else np.nan
def feats(lst):
    return np.array([[len(s),gcval(s)] for _,s in lst])
dcr_f=feats(dcr_h)
print("DCR human exon: median len",np.median(dcr_f[:,0]),"median GC",round(np.median(dcr_f[:,1]),1))
~~~~

### [msg 269] python · env: dcr-tierA · Scoring latent potential in matched non-DCR exons

~~~~python
# First: observed DCR human exons — what fraction have the dual-coding POTENTIAL by the esl rule?
dcr_dc=np.array([is_dualcoding_seq(s) for _,s in dcr_h])
print(f"Human DCR exons with dual-coding potential (esl rule): {100*dcr_dc.mean():.1f}%")

# NN GC+length matching of non-DCR pools to DCR exons, then score latent potential
def match_and_score(pool_lst, pool_feats, target_feats, n_draw=None, seed=1):
    mu=target_feats.mean(0); sd=target_feats.std(0)
    tree=cKDTree((pool_feats-mu)/sd)
    q=(target_feats-mu)/sd
    k=1
    dd,ii=tree.query(q,k=k)
    idx=ii if k==1 else ii[:,0]
    seqs=[pool_lst[j][1] for j in idx]
    flags=np.array([is_dualcoding_seq(s) for s in seqs])
    mf=np.array([[len(s),gcval(s)] for s in seqs])
    return flags, mf, idx

ndg_f=feats(nDCR_dg); nng_f=feats(nDCR_ng)
# match non-DCR in DCR genes
f_dg,mf_dg,_=match_and_score(nDCR_dg,ndg_f,dcr_f)
f_ng,mf_ng,_=match_and_score(nDCR_ng,nng_f,dcr_f)
print(f"\nnon-DCR exon in DCR genes (GC+len matched to DCR): latent potential {100*f_dg.mean():.1f}%")
print(f"  matched GC median {np.median(mf_dg[:,1]):.1f} (DCR {np.median(dcr_f[:,1]):.1f}); len {np.median(mf_dg[:,0]):.0f} (DCR {np.median(dcr_f[:,0]):.0f})")
print(f"non-DCR exon in non-DCR genes (GC+len matched to DCR): latent potential {100*f_ng.mean():.1f}%")
print(f"  matched GC median {np.median(mf_ng[:,1]):.1f}; len {np.median(mf_ng[:,0]):.0f}")
~~~~

### [msg 271] python · env: dcr-tierA · Bootstrapping latent-potential fractions

~~~~python
# Bootstrap CIs over DCR exons (resample DCR exons, re-match). Use k=10 NN pool, random pick per rep.
def matched_pool(pool_lst,pool_feats,target_feats,k=10):
    mu=target_feats.mean(0);sd=target_feats.std(0)
    tree=cKDTree((pool_feats-mu)/sd)
    dd,ii=tree.query((target_feats-mu)/sd,k=k)
    return ii  # (n_target, k) indices into pool
ii_dg=matched_pool(nDCR_dg,ndg_f,dcr_f,10)
ii_ng=matched_pool(nDCR_ng,nng_f,dcr_f,10)
# precompute dual-coding flags for the union of matched pool indices (cache)
def cache_flags(pool_lst, ii):
    uniq=np.unique(ii)
    return {j:is_dualcoding_seq(pool_lst[j][1]) for j in uniq}
fl_dg=cache_flags(nDCR_dg,ii_dg); fl_ng=cache_flags(nDCR_ng,ii_ng)

rng=np.random.default_rng(3)
B=500; nt=len(dcr_f)
boot_dcr=[];boot_dg=[];boot_ng=[]
for b in range(B):
    rows=rng.integers(0,nt,nt)
    boot_dcr.append(dcr_dc[rows].mean())
    kd=rng.integers(0,10,nt); kn=rng.integers(0,10,nt)
    boot_dg.append(np.mean([fl_dg[ii_dg[r,kd[i]]] for i,r in enumerate(rows)]))
    boot_ng.append(np.mean([fl_ng[ii_ng[r,kn[i]]] for i,r in enumerate(rows)]))
def ci(x):x=np.array(x)*100;return round(x.mean(),1),round(np.percentile(x,2.5),1),round(np.percentile(x,97.5),1)
print("DCR exons dual-coding potential:      %.1f%% (CI %.1f-%.1f)"%ci(boot_dcr))
print("non-DCR exon, DCR gene (matched):     %.1f%% (CI %.1f-%.1f)"%ci(boot_dg))
print("non-DCR exon, non-DCR gene (matched): %.1f%% (CI %.1f-%.1f)"%ci(boot_ng))
lat=pd.DataFrame({'class':['DCR_exon','nonDCR_exon_DCR_gene_matched','nonDCR_exon_nonDCR_gene_matched'],
 'latent_potential_pct':[ci(boot_dcr)[0],ci(boot_dg)[0],ci(boot_ng)[0]],
 'ci_low':[ci(boot_dcr)[1],ci(boot_dg)[1],ci(boot_ng)[1]],
 'ci_high':[ci(boot_dcr)[2],ci(boot_dg)[2],ci(boot_ng)[2]]})
lat.to_csv(f"{WS}/tierA/latent_potential_nonDCR.csv",index=False)
~~~~

### [msg 273] python · env: dcr-tierA · Rendering latent-potential figure

~~~~python
fig,ax=plt.subplots(figsize=(4.6,3.6))
labs=["DCR exon\n(observed)","non-DCR exon,\nDCR gene\n(GC+len matched)","non-DCR exon,\nnon-DCR gene\n(GC+len matched)"]
vals=[lat.latent_potential_pct[i] for i in range(3)]
los=[vals[i]-lat.ci_low[i] for i in range(3)];his=[lat.ci_high[i]-vals[i] for i in range(3)]
cols_=["#FFA800","#FFF0BD","#cccccc"]
ax.bar([0,1,2],vals,yerr=[los,his],color=cols_,width=0.62,capsize=3,
       error_kw=dict(lw=1),edgecolor='#888',linewidth=0.5)
for i,v in enumerate(vals):ax.text(i,v+2,f"{v:.1f}%",ha='center',fontsize=7.5)
ax.set_xticks([0,1,2]);ax.set_xticklabels(labs,fontsize=6.5)
ax.set_ylabel("% exons with latent dual-coding potential")
ax.set_ylim(0,60)
ax.axhline(41,color='#c0392b',ls='--',lw=0.8,alpha=0.6)
ax.set_title("R1-M2a: latent dual-coding potential is common\nin ordinary exons \u2014 necessary but not sufficient",loc="left",fontsize=8)
fig.tight_layout();fig.savefig(f"{WS}/tierA/latent_potential_nonDCR.png",dpi=300)
print("saved")
plt.close(fig)
print(lat.to_string(index=False))
~~~~
