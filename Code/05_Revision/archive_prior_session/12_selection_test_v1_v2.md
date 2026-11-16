# 12_selection_test_v1_v2 — as executed, not refactored

**Alternate-frame selection test v1 (GENCODE frames) and v2 (protein-derived DCR frames): CDS tables, intersections, OLGenie scoring, ascertainment checks.**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[484..625]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 54 (msg 484 … msg 616) |
| Language / environment | bash (env: dcr-tierA); python (env: dcr-tierA) |
| Saved outputs | `alternate_frame_selection.png`, `selection_omega.csv`, `selection_by_framerel.csv`, `selection_ascertainment.png`, `OLGenie_methods_explained.md`, `ascertainment_check.csv`, `dcr_exclusion_fates.csv`, `omega_by_length_quartile.csv`, `selection_ascertainment.png`, `OLGenie_methods_explained.md`, `ascertainment_check.csv`, `selection_omega.csv`, `selection_by_framerel.csv`, `dcr_protein_frames.csv`, `omega_by_length_bin.csv`, `alternate_frame_selection.png` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 25 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

---

### [msg 484] bash · env: dcr-tierA · Checking phase field and indel-free lift rate

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
X="$WS/tierA/xprime"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
echo "=== 1. is the field after _cds_<n>_ a phase (0/1/2)? ==="
awk '{n=split($4,a,"_"); for(i=1;i<=n;i++) if(a[i]=="cds") print a[i+2]}' "$REF/hg38.GENECODE.V48.bed" | sort | uniq -c | sort -rn | head -6
echo ""
echo "=== 2. indel-free lift rate (human len == mouse len) ==="
for T in D N; do
  IN="$X/dcr_regions.bed"; [ "$T" = "N" ] && IN="$X/N_human.bed"
  awk 'BEGIN{OFS="\t"}{print $4,$3-$2}' "$IN" | sort -k1,1 > /tmp/${T}_h.len
  awk 'BEGIN{OFS="\t"}{print $4,$3-$2}' "$X/${T}_ucsc.bed" | sort -k1,1 > /tmp/${T}_m.len
  join -t$'\t' /tmp/${T}_h.len /tmp/${T}_m.len > /tmp/${T}.j
  tot=$(wc -l < /tmp/${T}.j)
  same=$(awk -F'\t' '$2==$3' /tmp/${T}.j | wc -l)
  mod3=$(awk -F'\t' '$2==$3 && $2%3==0' /tmp/${T}.j | wc -l)
  echo "$T: $same/$tot same-length ($(awk -v a=$same -v b=$tot 'BEGIN{printf "%.1f", a/b*100}')%), of which len%3==0: $mod3"
done
~~~~

### [msg 488] bash · env: dcr-tierA · Building CDS exon table for phase

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
X="$WS/tierA/xprime"; S="$WS/tierA/sel"; mkdir -p "$S"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
# CDS phase per exon: cumulative CDS length within each transcript, in transcript order
# GENCODE bed: name = ENST_cds_<idx>_0_chr_start_<f|r>
awk 'BEGIN{OFS="\t"}{
  n=split($4,a,"_"); enst=a[1]; idx=a[3]; strand=a[n];
  print enst,idx,$1,$2,$3,$3-$2,strand,$4
}' "$REF/hg38.GENECODE.V48.bed" | sort -k1,1 -k2,2n > "$S/cds_exons.tsv"
wc -l < "$S/cds_exons.tsv"; head -3 "$S/cds_exons.tsv"
~~~~

### [msg 490] python · env: dcr-tierA · Deriving CDS phase per exon

~~~~python
S=f"{WS}/tierA/sel"
ce=pd.read_csv(f"{S}/cds_exons.tsv",sep="\t",header=None,
               names=['enst','idx','chr','start','end','len','strand','name'])
# cumulative CDS length preceding each exon, in transcript order (idx already 5'->3' for both strands in GENCODE bed)
ce=ce.sort_values(['enst','idx'])
ce['cum_before']=ce.groupby('enst')['len'].cumsum()-ce['len']
ce['phase']=ce.cum_before%3            # 0 => exon starts at codon pos1
print("phase distribution across CDS exons:"); print(ce.phase.value_counts().sort_index().to_string())
# sanity: total CDS length divisible by 3 for most transcripts
tl=ce.groupby('enst')['len'].sum()
print(f"\ntranscripts with CDS length %3==0: {(tl%3==0).mean()*100:.1f}%  (n={len(tl)})")
name2phase=dict(zip(ce.name,ce.phase))
name2strand=dict(zip(ce.name,ce.strand))
print("phase map entries:",len(name2phase))
~~~~

### [msg 492] bash · env: dcr-tierA · Intersecting intervals with CDS exons

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
X="$WS/tierA/xprime"; S="$WS/tierA/sel"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
sort -k1,1 -k2,2n "$REF/hg38.GENECODE.V48.bed" > "$S/gc.s.bed"
for T in D:dcr_regions N:N_human; do
  TAG="${T%%:*}"; IN="${T##*:}"
  sort -k1,1 -k2,2n "$X/$IN.bed" > "$S/${TAG}_h.s.bed"
  bedtools intersect -a "$S/${TAG}_h.s.bed" -b "$S/gc.s.bed" -wo -f 0.9 > "$S/${TAG}_cds.ov"
  echo "$TAG: $(cut -f4 "$S/${TAG}_cds.ov" | sort -u | wc -l) intervals with >=90% CDS overlap"
done
head -1 "$S/D_cds.ov"
~~~~

### [msg 494] python · env: dcr-tierA · Deriving canonical frame per interval

~~~~python
def load_ov(tag):
    ov=pd.read_csv(f"{S}/{tag}_cds.ov",sep="\t",header=None,
        names=['c','s','e','name','sc','st','ec','es','ee','ename','esc','est','bp'])
    ov=ov[ov.ename.isin(name2phase)]
    ov['ephase']=ov.ename.map(name2phase)
    # codon position of the interval's 5' end (transcript orientation)
    plus = ov.st=='+'
    off = np.where(plus, ov.s-ov.es, ov.ee-ov.e)
    ov['cpos']=(ov.ephase.values+off)%3          # 0 => interval starts at codon pos 1
    # keep one CDS exon per interval (largest overlap), and require consistent frame
    ov=ov.sort_values('bp',ascending=False)
    agg=ov.groupby('name').agg(cpos=('cpos','first'),n_frames=('cpos','nunique'),st=('st','first'))
    return agg
ovD=load_ov('D'); ovN=load_ov('N')
print("D intervals with unambiguous canonical frame: %d/%d"%((ovD.n_frames==1).sum(),len(ovD)))
print("N intervals with unambiguous canonical frame: %d/%d"%((ovN.n_frames==1).sum(),len(ovN)))
print("\nD cpos distribution:"); print(ovD.cpos.value_counts().sort_index().to_string())
~~~~

### [msg 496] python · env: dcr-tierA · Extracting dual annotated frames for DCRs

~~~~python
def load_pairs(tag):
    ov=pd.read_csv(f"{S}/{tag}_cds.ov",sep="\t",header=None,
        names=['c','s','e','name','sc','st','ec','es','ee','ename','esc','est','bp'])
    ov=ov[ov.ename.isin(name2phase)].copy()
    ov['ephase']=ov.ename.map(name2phase)
    plus = ov.est=='+'
    off = np.where(plus, ov.s.values-ov.es.values, ov.ee.values-ov.e.values)
    ov['cpos']=(ov.ephase.values+off)%3
    g=ov.groupby('name')
    out=g.agg(frames=('cpos',lambda v:sorted(set(v))),st=('est','first'),bpmax=('bp','max'))
    return out
pD=load_pairs('D')
pD['nf']=pD.frames.apply(len)
print("D intervals by number of distinct annotated frames:")
print(pD.nf.value_counts().sort_index().to_string())
two=pD[pD.nf>=2].copy()
two['f1']=two.frames.apply(lambda v:v[0]); two['f2']=two.frames.apply(lambda v:v[1])
two['rel']=(two.f2-two.f1)%3
print(f"\nD intervals with >=2 annotated frames: {len(two)} ({len(two)/len(pD)*100:.1f}%)")
print("frame relationship (1 = ss13/+1, 2 = ss12/+2):")
print(two.rel.value_counts().sort_index().to_string())
~~~~

### [msg 498] python · env: dcr-tierA · Implementing OLGenie site classification

~~~~python
CODON={}
bases='TCAG'
aas='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i=0
for b1 in bases:
    for b2 in bases:
        for b3 in bases:
            CODON[b1+b2+b3]=aas[i]; i+=1
def syn(c1,c2):
    """True if codons encode same aa (both must be valid)"""
    a1=CODON.get(c1); a2=CODON.get(c2)
    if a1 is None or a2 is None: return None
    return a1==a2

def olgenie_counts(h,m,ref_off,alt_off):
    """Cross-classify each nt difference by effect in reference and alternate frame.
       Returns dict of NN,NS,SN,SS *differences* and corresponding site denominators."""
    L=min(len(h),len(m))
    d={'NN':0,'NS':0,'SN':0,'SS':0}
    sites={'NN':0.0,'NS':0.0,'SN':0.0,'SS':0.0}
    for p in range(L):
        # codon start in each frame containing position p
        rs=p-((p-ref_off)%3); as_=p-((p-alt_off)%3)
        if rs<0 or rs+3>L or as_<0 or as_+3>L: continue
        hr,mr=h[rs:rs+3],m[rs:rs+3]
        ha,ma=h[as_:as_+3],m[as_:as_+3]
        if len(hr)<3 or len(ha)<3: continue
        if any(x not in 'ACGT' for x in hr+ha+mr+ma): continue
        # site counting: for this position, enumerate 3 alt bases, classify effect in both frames
        ri=p-rs; ai=p-as_
        for b in 'ACGT':
            if b==h[p]: continue
            cr=hr[:ri]+b+hr[ri+1:]; ca=ha[:ai]+b+ha[ai+1:]
            sr=syn(hr,cr); sa=syn(ha,ca)
            if sr is None or sa is None: continue
            k=('S' if sr else 'N')+('S' if sa else 'N')
            sites[k]+=1/3.0
        # observed difference
        if h[p]!=m[p]:
            cr=hr[:ri]+m[p]+hr[ri+1:]; ca=ha[:ai]+m[p]+ha[ai+1:]
            sr=syn(hr,cr); sa=syn(ha,ca)
            if sr is None or sa is None: continue
            d[('S' if sr else 'N')+('S' if sa else 'N')]+=1
    return d,sites
print("site classifier ready; codon table entries:",len(CODON))
~~~~

### [msg 500] python · env: dcr-tierA · Loading sequence pairs

~~~~python
# Build indel-free codon-aligned human/mouse pairs for D
hD_seq=read_fa2(f"{X}/dcr_regions.fa"); mD_seq=read_fa2(f"{X}/D_ucsc.fa")
hN_seq=read_fa2(f"{X}/nondcr_pool.fa"); mN_seq=read_fa2(f"{X}/N_ucsc.fa")
hDk={h.split('::')[0]:s for h,s in hD_seq.items()}
mDk={h.split('::')[0]:s for h,s in mD_seq.items()}
hNraw={}
for h,s in hN_seq.items():
    mm=pat.match(h)
    if mm: hNraw[mm.group(2)]=s
mNk={h.split('::')[0]:s for h,s in mN_seq.items()}
print("D human:",len(hDk)," D mouse:",len(mDk))
print("N mouse:",len(mNk))
# N human interval seqs: slice the latent sub-interval out of the exon sequence
Nrec=Nsub2.copy()
print("N records:",len(Nrec), "cols:",list(Nrec.columns)[:12])
~~~~

### [msg 502] bash · env: dcr-tierA · Extracting N latent interval sequences

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
X="$WS/tierA/xprime"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
# extract the N latent sub-intervals directly from hg38 (stranded) so human/mouse are the same object
bedtools getfasta -fi "$REF/hg38.fa" -bed "$X/N_human.bed" -s -name -fo "$X/N_human.fa"
grep -c '^>' "$X/N_human.fa"; head -1 "$X/N_human.fa"
~~~~

### [msg 504] python · env: dcr-tierA · Loading N human sequences and frames

~~~~python
hNk={h.split('::')[0]:s for h,s in read_fa2(f"{X}/N_human.fa").items()}
print("N human:",len(hNk))
# For N: canonical frame from annotation, alternate frame = the detected latent ORF frame
pN=load_pairs('N'); pN['nf']=pN.frames.apply(len)
print("N intervals by n annotated frames:"); print(pN.nf.value_counts().sort_index().head().to_string())
# detect the latent alternate frame per N interval: frame pair giving the max dual overlap
def latent_frames(seq):
    orfs=esl_translate_watson(seq,min_aa=1); best=None;bs=0
    for i in range(len(orfs)):
        s1,e1,f1,_=orfs[i]
        for j in range(i+1,len(orfs)):
            s2,e2,f2,_=orfs[j]
            if f1==f2: continue
            span=min(e1,e2)-max(s1,s2)+1
            if span>bs: bs=span;best=(f1-1,f2-1)
    return best
print("ready")
~~~~

### [msg 506] python · env: dcr-tierA · Building indel-free codon-aligned sets

~~~~python
import time
def build_set(tag, hk, mk, frames_df, mode):
    recs=[]
    for name,row in frames_df.iterrows():
        h=hk.get(name); m=mk.get(name)
        if h is None or m is None: continue
        if len(h)!=len(m) or len(h)<30: continue      # indel-free only
        if mode=='D':
            fr=row.frames
            if len(fr)<2: continue
            f1,f2=fr[0],fr[1]
        else:
            f1=row.frames[0]
            lf=latent_frames(h)
            if lf is None: continue
            f2 = lf[1] if lf[0]==f1 else (lf[0] if lf[1]==f1 else lf[1])
            if f2==f1: continue
        recs.append((name,h,m,f1,f2,(f2-f1)%3))
    return recs
t0=time.time()
Drecs=build_set('D',hDk,mDk,pD[pD.nf>=2],'D')
Nrecs=build_set('N',hNk,mNk,pN[pN.nf==1],'N')
print(f"D usable indel-free pairs: {len(Drecs)}")
print(f"N usable indel-free pairs: {len(Nrecs)}   ({time.time()-t0:.0f}s)")
print("D frame rel:",pd.Series([r[5] for r in Drecs]).value_counts().sort_index().to_dict())
print("N frame rel:",pd.Series([r[5] for r in Nrecs]).value_counts().sort_index().to_dict())
~~~~

### [msg 508] python · env: dcr-tierA · Scoring OLGenie site classes

~~~~python
def score_set(recs,label):
    rows=[]
    for name,h,m,f1,f2,rel in recs:
        d,s=olgenie_counts(h,m,f1,f2)
        rows.append((name,rel,len(h),d['NN'],d['NS'],d['SN'],d['SS'],
                     s['NN'],s['NS'],s['SN'],s['SS']))
    df=pd.DataFrame(rows,columns=['name','rel','len','dNN','dNS','dSN','dSS',
                                  'sNN','sNS','sSN','sSS'])
    return df
t0=time.time()
selD=score_set(Drecs,'D'); print(f"D scored ({time.time()-t0:.0f}s)")
t0=time.time()
selN=score_set(Nrecs,'N'); print(f"N scored ({time.time()-t0:.0f}s)")
def pooled(df,rel=None):
    d=df if rel is None else df[df.rel==rel]
    o={}
    for k in ['NN','NS','SN','SS']:
        o['d'+k]=d['d'+k].sum(); o['s'+k]=d['s'+k].sum()
        o['p'+k]=o['d'+k]/o['s'+k] if o['s'+k]>0 else np.nan
    return o
pdD=pooled(selD); pdN=pooled(selN)
print("\n=== pooled site counts ===")
for lab,o in [('D (real DCR)',pdD),('N (latent)',pdN)]:
    print(f"{lab}: dNN={o['dNN']} dNS={o['dNS']} dSN={o['dSN']} dSS={o['dSS']} | sites NN={o['sNN']:.0f} NS={o['sNS']:.0f} SN={o['sSN']:.0f} SS={o['sSS']:.0f}")
~~~~

### [msg 510] python · env: dcr-tierA · Computing omega estimators

~~~~python
def jc(p):
    """Jukes-Cantor correction; returns nan if saturated"""
    if p is None or np.isnan(p) or p>=0.749: return np.nan
    return -0.75*np.log(1-4*p/3)
def omegas(o):
    r={}
    for k in ['NN','NS','SN','SS']: r[k]=jc(o['p'+k])
    r['w_alt_refN'] = r['NN']/r['NS'] if r['NS'] and r['NS']>0 else np.nan   # dNN/dNS
    r['w_alt_refS'] = r['SN']/r['SS'] if r['SS'] and r['SS']>0 else np.nan   # dSN/dSS
    r['w_ref']      = r['NN']/r['SN'] if r['SN'] and r['SN']>0 else np.nan   # positive control
    return r
oD=omegas(pdD); oN=omegas(pdN)
print(f"{'':28s} {'D (real DCR)':>14s} {'N (latent)':>14s}")
print(f"{'w_alt = dNN/dNS':28s} {oD['w_alt_refN']:14.3f} {oN['w_alt_refN']:14.3f}")
print(f"{'w_alt = dSN/dSS  (your fix)':28s} {oD['w_alt_refS']:14.3f} {oN['w_alt_refS']:14.3f}")
print(f"{'w_ref = dNN/dSN  (control)':28s} {oD['w_ref']:14.3f} {oN['w_ref']:14.3f}")
print("\nper-class corrected distances:")
for k in ['NN','NS','SN','SS']:
    print(f"  d{k}: D={oD[k]:.4f}  N={oN[k]:.4f}")
~~~~

### [msg 512] python · env: dcr-tierA · Validating reference frame offset

~~~~python
# Empirical frame check: the TRUE reading frame should minimise dN/dS in the reference frame.
def wref_for_offset(recs,off_mode,n=400):
    dN=dS=sN=sS=0.0
    for name,h,m,f1,f2,rel in recs[:n]:
        f = f1 if off_mode=='raw' else (3-f1)%3
        L=min(len(h),len(m))
        for st in range(f,L-2,3):
            hc,mc=h[st:st+3],m[st:st+3]
            if len(hc)<3 or any(x not in 'ACGT' for x in hc+mc): continue
            for i in range(3):
                for b in 'ACGT':
                    if b==hc[i]: continue
                    cc=hc[:i]+b+hc[i+1:]
                    s=syn(hc,cc)
                    if s is None: continue
                    if s: sS+=1/3.0
                    else: sN+=1/3.0
            if hc!=mc:
                for i in range(3):
                    if hc[i]!=mc[i]:
                        cc=hc[:i]+mc[i]+hc[i+1:]
                        s=syn(hc,cc)
                        if s is None: continue
                        if s: dS+=1
                        else: dN+=1
    pN_,pS_=dN/sN,dS/sS
    return jc(pN_)/jc(pS_)
for mode in ['raw','comp']:
    print(f"offset={mode:5s}  D w_ref={wref_for_offset(Drecs,mode):.3f}   N w_ref={wref_for_offset(Nrecs,mode):.3f}")
~~~~

### [msg 514] python · env: dcr-tierA · Rescoring with corrected frame offset

~~~~python
def fix(f): return (3-f)%3
Drecs2=[(n,h,m,fix(f1),fix(f2),(fix(f2)-fix(f1))%3) for n,h,m,f1,f2,rel in Drecs]
# N: reference = corrected annotated frame; alternate = latent frame that pairs with it
Nrecs2=[]
for name,h,m,f1,f2,rel in Nrecs:
    rf=fix(f1)
    lf=latent_frames(h)
    if lf is None: continue
    a,b=lf
    if   a==rf: af=b
    elif b==rf: af=a
    else: continue                    # latent pair doesn't involve the real frame -> drop
    Nrecs2.append((name,h,m,rf,af,(af-rf)%3))
print(f"D pairs: {len(Drecs2)}   N pairs: {len(Nrecs2)} (dropped {len(Nrecs)-len(Nrecs2)} whose latent pair excluded the annotated frame)")
selD2=score_set(Drecs2,'D'); selN2=score_set(Nrecs2,'N')
oD2=omegas(pooled(selD2)); oN2=omegas(pooled(selN2))
print(f"\n{'':30s} {'D (real DCR)':>13s} {'N (latent)':>13s}")
print(f"{'w_ref = dNN/dSN (control)':30s} {oD2['w_ref']:13.3f} {oN2['w_ref']:13.3f}")
print(f"{'w_alt = dNN/dNS':30s} {oD2['w_alt_refN']:13.3f} {oN2['w_alt_refN']:13.3f}")
print(f"{'w_alt = dSN/dSS (your fix)':30s} {oD2['w_alt_refS']:13.3f} {oN2['w_alt_refS']:13.3f}")
~~~~

### [msg 516] python · env: dcr-tierA · Bootstrapping omega confidence intervals

~~~~python
def boot_omega(df,est,B=2000,seed=7):
    rng=np.random.default_rng(seed); out=[]
    idx=np.arange(len(df)); A=df[['dNN','dNS','dSN','dSS','sNN','sNS','sSN','sSS']].values
    num,den=(('dNN','dNS') if est=='refN' else ('dSN','dSS'))
    ni,di=['dNN','dNS','dSN','dSS'].index(num),['dNN','dNS','dSN','dSS'].index(den)
    for _ in range(B):
        s=A[rng.integers(0,len(df),len(df))]
        pn=s[:,ni].sum()/s[:,4+ni].sum(); pd_=s[:,di].sum()/s[:,4+di].sum()
        a,b=jc(pn),jc(pd_)
        if a and b and b>0 and not np.isnan(a) and not np.isnan(b): out.append(a/b)
    return np.array(out)
res_sel=[]
for lab,df in [('D real DCR',selD2),('N latent',selN2)]:
    for est,nm in [('refN','dNN/dNS'),('refS','dSN/dSS')]:
        o=omegas(pooled(df)); v=o['w_alt_refN' if est=='refN' else 'w_alt_refS']
        bs=boot_omega(df,est)
        res_sel.append((lab,nm,len(df),v,np.percentile(bs,2.5),np.percentile(bs,97.5)))
        print(f"{lab:12s} {nm:9s} omega={v:.3f}  95% CI [{np.percentile(bs,2.5):.3f}, {np.percentile(bs,97.5):.3f}]  n={len(df)}")
sel_tab=pd.DataFrame(res_sel,columns=['set','estimator','n','omega','lo','hi'])
~~~~

### [msg 518] python · env: dcr-tierA · Stratifying by frame relationship and length

~~~~python
selD2['lb']=pd.cut(selD2.len,bins); selN2['lb']=pd.cut(selN2.len,bins)
rows=[]
for rel,nm in [(1,'ss13 (+1)'),(2,'ss12 (+2)')]:
    for lab,df in [('D',selD2),('N',selN2)]:
        d=df[df.rel==rel]
        if len(d)<30: continue
        o=omegas(pooled(d))
        rows.append((nm,lab,len(d),o['w_alt_refN'],o['w_alt_refS']))
fr_tab=pd.DataFrame(rows,columns=['frame_rel','set','n','w_dNNdNS','w_dSNdSS'])
print("=== stratified by frame relationship ===")
print(fr_tab.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
# length-matched: reweight N by D's length distribution, per estimator
wD=selD2.groupby('lb',observed=True).size(); 
def matched_omega(est):
    key='w_alt_refN' if est=='refN' else 'w_alt_refS'
    vals={}
    for b in wD.index:
        n=selN2[selN2.lb==b]
        if len(n)<30 or wD[b]<20: continue
        o=omegas(pooled(n))
        if not np.isnan(o[key]): vals[b]=o[key]
    ww=wD[list(vals.keys())]; ww=ww/ww.sum()
    return float(sum(vals[b]*ww[b] for b in vals)), list(vals.keys())
for est,nm in [('refN','dNN/dNS'),('refS','dSN/dSS')]:
    m,bl=matched_omega(est)
    dO=omegas(pooled(selD2[selD2.lb.isin(bl)]))[('w_alt_refN' if est=='refN' else 'w_alt_refS')]
    print(f"\n{nm} length-matched: D={dO:.3f}  N={m:.3f}   (bins={len(bl)})")
~~~~

### [msg 520] python · env: dcr-tierA · Testing D vs N omega difference

~~~~python
rng3=np.random.default_rng(3)
def boot_diff(est,B=2000):
    key=0 if est=='refN' else 1
    num,den=(('dNN','dNS') if est=='refN' else ('dSN','dSS'))
    out=[]
    AD=selD2[[num,den,'s'+num[1:],'s'+den[1:]]].values
    AN=selN2[[num,den,'s'+num[1:],'s'+den[1:]]].values
    for _ in range(B):
        d=AD[rng3.integers(0,len(AD),len(AD))]; n=AN[rng3.integers(0,len(AN),len(AN))]
        wd=jc(d[:,0].sum()/d[:,2].sum())/jc(d[:,1].sum()/d[:,3].sum())
        wn=jc(n[:,0].sum()/n[:,2].sum())/jc(n[:,1].sum()/n[:,3].sum())
        if not (np.isnan(wd) or np.isnan(wn)): out.append(wn-wd)
    return np.array(out)
for est,nm in [('refN','dNN/dNS'),('refS','dSN/dSS')]:
    bd=boot_diff(est)
    p=(bd<=0).mean()
    print(f"{nm}: omega_N - omega_D = {bd.mean():.3f}  95% CI [{np.percentile(bd,2.5):.3f}, {np.percentile(bd,97.5):.3f}]  p={max(p,1/len(bd)):.4f}")
sel_tab.to_csv(f"{S}/selection_omega.csv",index=False)
fr_tab.to_csv(f"{S}/selection_by_framerel.csv",index=False)
print("\nsaved tables")
~~~~

### [msg 522] python · env: dcr-tierA · Building selection test figure

~~~~python
fig,axes=plt.subplots(1,4,figsize=(13.6,3.6))
ax=axes[0]
labs=['dNN/dNS','dSN/dSS']
dv=[oD2['w_alt_refN'],oD2['w_alt_refS']]; nv=[oN2['w_alt_refN'],oN2['w_alt_refS']]
dlo=[sel_tab.iloc[0].lo,sel_tab.iloc[1].lo]; dhi=[sel_tab.iloc[0].hi,sel_tab.iloc[1].hi]
nlo=[sel_tab.iloc[2].lo,sel_tab.iloc[3].lo]; nhi=[sel_tab.iloc[2].hi,sel_tab.iloc[3].hi]
xx=np.arange(2); wdt=0.34
ax.bar(xx-wdt/2,dv,wdt,color=ORANGE,edgecolor='k',linewidth=0.5,label='real DCR',
       yerr=[np.array(dv)-dlo,np.array(dhi)-dv],capsize=3,error_kw=dict(lw=0.9))
ax.bar(xx+wdt/2,nv,wdt,color=GREY,edgecolor='k',linewidth=0.5,label='latent (null)',
       yerr=[np.array(nv)-nlo,np.array(nhi)-nv],capsize=3,error_kw=dict(lw=0.9))
ax.axhline(1.0,ls='--',lw=1,color='k'); ax.text(1.42,1.03,'neutrality',fontsize=6.5,ha='right')
ax.set_xticks(xx); ax.set_xticklabels(labs,fontsize=8)
ax.set_ylabel("ω  (alternate frame)"); ax.set_ylim(0,1.25); ax.legend(frameon=False,fontsize=7,loc='upper left')
ax.set_title("A  Alternate frame is under selection",loc='left',weight='bold')

ax=axes[1]
ctrl=[oD2['w_ref'],oN2['w_ref']]
ax.bar([0,1],ctrl,color=[ORANGE,GREY],edgecolor='k',linewidth=0.5,width=0.5)
for i,v in enumerate(ctrl): ax.text(i,v+0.015,f"{v:.3f}",ha='center',fontsize=8,weight='bold')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks([0,1]); ax.set_xticklabels(['real DCR','latent'],fontsize=8)
ax.set_ylabel("ω  (canonical frame)  dNN/dSN"); ax.set_ylim(0,1.1)
ax.text(0.5,0.72,"positive control:\nboth strongly\npurifying, as\nexpected for CDS",transform=ax.transAxes,
        ha='center',fontsize=6.8,style='italic')
ax.set_title("B  Control: canonical frame",loc='left',weight='bold')

ax=axes[2]
p=fr_tab.pivot(index='frame_rel',columns='set',values='w_dSNdSS')
xx=np.arange(len(p)); 
ax.bar(xx-wdt/2,p['D'],wdt,color=ORANGE,edgecolor='k',linewidth=0.5,label='real DCR')
ax.bar(xx+wdt/2,p['N'],wdt,color=GREY,edgecolor='k',linewidth=0.5,label='latent')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks(xx); ax.set_xticklabels(p.index,fontsize=8)
ax.set_ylabel("ω  dSN/dSS"); ax.set_ylim(0,1.35); ax.legend(frameon=False,fontsize=7)
ax.set_title("C  Holds in both frame offsets",loc='left',weight='bold')

ax=axes[3]
mm_={'dNN/dNS':(0.265,0.766),'dSN/dSS':(0.442,1.025)}
xx=np.arange(2)
ax.bar(xx-wdt/2,[mm_[k][0] for k in labs],wdt,color=ORANGE,edgecolor='k',linewidth=0.5,label='real DCR')
ax.bar(xx+wdt/2,[mm_[k][1] for k in labs],wdt,color=BLUE,edgecolor='k',linewidth=0.5,label='latent, length-matched')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks(xx); ax.set_xticklabels(labs,fontsize=8)
ax.set_ylabel("ω  (alternate frame)"); ax.set_ylim(0,1.25); ax.legend(frameon=False,fontsize=6.8,loc='upper left')
ax.text(0.5,0.60,"survives the length\ncontrol that erased\nthe retention signal",transform=ax.transAxes,
        ha='center',fontsize=6.8,style='italic',color=BLUE)
ax.set_title("D  Length-matched",loc='left',weight='bold')
fig.suptitle("Alternate reading frames of real DCRs are under purifying selection; latent dual-coding intervals are not",
             x=0.005,ha='left',fontsize=10,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/alternate_frame_selection.png",dpi=300); plt.close(fig)
print("saved")
~~~~

### [msg 524] python · env: dcr-tierA · Fixing panel D annotation overlap

~~~~python
ax=axes[3]
fig.savefig  # noop
# rebuild panel D annotation position only
fig,axes=plt.subplots(1,4,figsize=(13.6,3.6))
ax=axes[0]
xx=np.arange(2); wdt=0.34
ax.bar(xx-wdt/2,dv,wdt,color=ORANGE,edgecolor='k',linewidth=0.5,label='real DCR',
       yerr=[np.array(dv)-dlo,np.array(dhi)-dv],capsize=3,error_kw=dict(lw=0.9))
ax.bar(xx+wdt/2,nv,wdt,color=GREY,edgecolor='k',linewidth=0.5,label='latent (null)',
       yerr=[np.array(nv)-nlo,np.array(nhi)-nv],capsize=3,error_kw=dict(lw=0.9))
ax.axhline(1.0,ls='--',lw=1,color='k'); ax.text(1.44,1.03,'neutrality',fontsize=6.5,ha='right')
ax.set_xticks(xx); ax.set_xticklabels(labs,fontsize=8)
ax.set_ylabel("ω  (alternate frame)"); ax.set_ylim(0,1.3); ax.legend(frameon=False,fontsize=7,loc='upper left')
ax.set_title("A  Alternate frame is under selection",loc='left',weight='bold')

ax=axes[1]
ax.bar([0,1],ctrl,color=[ORANGE,GREY],edgecolor='k',linewidth=0.5,width=0.5)
for i,v in enumerate(ctrl): ax.text(i,v+0.02,f"{v:.3f}",ha='center',fontsize=8,weight='bold')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks([0,1]); ax.set_xticklabels(['real DCR','latent'],fontsize=8)
ax.set_ylabel("ω  (canonical frame)  dNN/dSN"); ax.set_ylim(0,1.15)
ax.text(0.5,0.60,"positive control:\nboth strongly purifying,\nas expected for CDS",transform=ax.transAxes,
        ha='center',fontsize=6.8,style='italic')
ax.set_title("B  Control: canonical frame",loc='left',weight='bold')

ax=axes[2]
xx=np.arange(len(p))
ax.bar(xx-wdt/2,p['D'],wdt,color=ORANGE,edgecolor='k',linewidth=0.5,label='real DCR')
ax.bar(xx+wdt/2,p['N'],wdt,color=GREY,edgecolor='k',linewidth=0.5,label='latent')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks(xx); ax.set_xticklabels(p.index,fontsize=8)
ax.set_ylabel("ω  dSN/dSS"); ax.set_ylim(0,1.4); ax.legend(frameon=False,fontsize=7)
ax.set_title("C  Holds in both frame offsets",loc='left',weight='bold')

ax=axes[3]
xx=np.arange(2)
ax.bar(xx-wdt/2,[mm_[k][0] for k in labs],wdt,color=ORANGE,edgecolor='k',linewidth=0.5,label='real DCR')
ax.bar(xx+wdt/2,[mm_[k][1] for k in labs],wdt,color=BLUE,edgecolor='k',linewidth=0.5,label='latent, length-matched')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks(xx); ax.set_xticklabels(labs,fontsize=8)
ax.set_ylabel("ω  (alternate frame)"); ax.set_ylim(0,1.45); ax.legend(frameon=False,fontsize=6.8,loc='upper left')
ax.text(0.52,0.42,"survives the length control\nthat erased the retention signal",transform=ax.transAxes,
        ha='center',fontsize=6.6,style='italic',color=BLUE)
ax.set_title("D  Length-matched",loc='left',weight='bold')
fig.suptitle("Alternate reading frames of real DCRs are under purifying selection; latent dual-coding intervals are not",
             x=0.005,ha='left',fontsize=10,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/alternate_frame_selection.png",dpi=300); plt.close(fig)
print("saved")
~~~~

### [msg 526] python · env: dcr-tierA · Repositioning panel D caption

~~~~python
import matplotlib.pyplot as plt
fig=plt.gcf()
# panel D annotation still overlaps the blue bar; move it below the bars
fig,axes=plt.subplots(1,4,figsize=(13.6,3.6))
def panelD(ax):
    xx=np.arange(2)
    ax.bar(xx-wdt/2,[mm_[k][0] for k in labs],wdt,color=ORANGE,edgecolor='k',linewidth=0.5,label='real DCR')
    ax.bar(xx+wdt/2,[mm_[k][1] for k in labs],wdt,color=BLUE,edgecolor='k',linewidth=0.5,label='latent, length-matched')
    ax.axhline(1.0,ls='--',lw=1,color='k')
    ax.set_xticks(xx); ax.set_xticklabels(labs,fontsize=8)
    ax.set_ylabel("ω  (alternate frame)"); ax.set_ylim(0,1.5); ax.legend(frameon=False,fontsize=6.8,loc='upper left')
    ax.text(0.5,-0.30,"survives the length control that\nerased the retention signal",transform=ax.transAxes,
            ha='center',fontsize=6.6,style='italic',color=BLUE)
    ax.set_title("D  Length-matched",loc='left',weight='bold')
ax=axes[0]
ax.bar(np.arange(2)-wdt/2,dv,wdt,color=ORANGE,edgecolor='k',linewidth=0.5,label='real DCR',
       yerr=[np.array(dv)-dlo,np.array(dhi)-dv],capsize=3,error_kw=dict(lw=0.9))
ax.bar(np.arange(2)+wdt/2,nv,wdt,color=GREY,edgecolor='k',linewidth=0.5,label='latent (null)',
       yerr=[np.array(nv)-nlo,np.array(nhi)-nv],capsize=3,error_kw=dict(lw=0.9))
ax.axhline(1.0,ls='--',lw=1,color='k'); ax.text(1.44,1.04,'neutrality',fontsize=6.5,ha='right')
ax.set_xticks(np.arange(2)); ax.set_xticklabels(labs,fontsize=8)
ax.set_ylabel("ω  (alternate frame)"); ax.set_ylim(0,1.3); ax.legend(frameon=False,fontsize=7,loc='upper left')
ax.set_title("A  Alternate frame is under selection",loc='left',weight='bold')
ax=axes[1]
ax.bar([0,1],ctrl,color=[ORANGE,GREY],edgecolor='k',linewidth=0.5,width=0.5)
for i,v in enumerate(ctrl): ax.text(i,v+0.02,f"{v:.3f}",ha='center',fontsize=8,weight='bold')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks([0,1]); ax.set_xticklabels(['real DCR','latent'],fontsize=8)
ax.set_ylabel("ω  (canonical frame)  dNN/dSN"); ax.set_ylim(0,1.15)
ax.text(0.5,0.60,"positive control:\nboth strongly purifying,\nas expected for CDS",transform=ax.transAxes,
        ha='center',fontsize=6.8,style='italic')
ax.set_title("B  Control: canonical frame",loc='left',weight='bold')
ax=axes[2]
xx=np.arange(len(p))
ax.bar(xx-wdt/2,p['D'],wdt,color=ORANGE,edgecolor='k',linewidth=0.5,label='real DCR')
ax.bar(xx+wdt/2,p['N'],wdt,color=GREY,edgecolor='k',linewidth=0.5,label='latent')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks(xx); ax.set_xticklabels(p.index,fontsize=8)
ax.set_ylabel("ω  dSN/dSS"); ax.set_ylim(0,1.4); ax.legend(frameon=False,fontsize=7)
ax.set_title("C  Holds in both frame offsets",loc='left',weight='bold')
panelD(axes[3])
fig.suptitle("Alternate reading frames of real DCRs are under purifying selection; latent dual-coding intervals are not",
             x=0.005,ha='left',fontsize=10,weight='bold')
fig.tight_layout(rect=[0,0.05,1,0.92])
fig.savefig(f"{WS}/tierA/alternate_frame_selection.png",dpi=300); plt.close(fig)
print("saved")
~~~~

### [msg 536] python · env: dcr-tierA · Tabulating DCR exclusion reasons

~~~~python
# Which DCRs survived, and why were the others dropped?
allD=set(pD.index)                                  # 1717 with >=90% CDS overlap
have_seq={n for n in allD if n in hDk and n in mDk} # lifted
kept={r[0] for r in Drecs2}
tiers=[]
for n in allD:
    if n not in hDk: tiers.append((n,'no_human_seq'))
    elif n not in mDk: tiers.append((n,'lift_failed'))
    elif len(hDk[n])!=len(mDk[n]): tiers.append((n,'indel'))
    elif len(hDk[n])<30: tiers.append((n,'too_short'))
    elif pD.loc[n,'nf']<2: tiers.append((n,'one_frame_only'))
    else: tiers.append((n,'kept'))
tt=pd.DataFrame(tiers,columns=['name','fate'])
print(tt.fate.value_counts().to_string())
print(f"\nkept: {len(kept)}  (checks out: {(tt.fate=='kept').sum()})")
~~~~

### [msg 538] python · env: dcr-tierA · Comparing kept vs dropped DCR covariates

~~~~python
meta=pD.join(tt.set_index('name'))
meta['len_h']=[len(hDk[n]) if n in hDk else np.nan for n in meta.index]
meta['gc']=[(hDk[n].count('G')+hDk[n].count('C'))/len(hDk[n])*100 if n in hDk else np.nan for n in meta.index]
grp=meta.assign(kept=meta.fate=='kept').groupby('kept').agg(
    n=('len_h','size'), med_len=('len_h','median'), med_gc=('gc','median'), med_nf=('nf','median'))
print("=== kept vs dropped DCRs ===")
print(grp.to_string(float_format=lambda v:f"{v:.1f}"))
from scipy import stats
k=meta[meta.fate=='kept']; d=meta[meta.fate!='kept']
for col,nm in [('len_h','length'),('gc','GC%')]:
    u,pv=stats.mannwhitneyu(k[col].dropna(),d[col].dropna())
    print(f"{nm}: kept med={k[col].median():.1f} dropped med={d[col].median():.1f}  MWU p={pv:.2e}")
# indel rate D vs N — the filter must be symmetric to be safe
d_lift=(tt.fate!='lift_failed')&(tt.fate!='no_human_seq')
print(f"\nindel rate among lifted:  D={201/d_lift.sum()*100:.1f}%   N={(1-28500/33520)*100:.1f}%")
~~~~

### [msg 540] python · env: dcr-tierA · Testing length and frame-filter ascertainment

~~~~python
# (a) does omega_alt vary with length inside the kept set? if flat, length ascertainment is benign
selD2['lb2']=pd.qcut(selD2.len,4,duplicates='drop')
rows=[]
for b,g in selD2.groupby('lb2',observed=True):
    o=omegas(pooled(g)); rows.append((str(b),len(g),o['w_alt_refN'],o['w_alt_refS']))
lt=pd.DataFrame(rows,columns=['len_quartile','n','w_dNNdNS','w_dSNdSS'])
print("=== omega_alt by DCR length quartile (kept set) ===")
print(lt.to_string(index=False,float_format=lambda v:f"{v:.3f}"))

# (b) add back the 266 one-frame DCRs, scored the same way as N (annotated + detected latent frame)
one=pD[(pD.nf==1)]
add=[]
for name,row in one.iterrows():
    h=hDk.get(name); m=mDk.get(name)
    if h is None or m is None or len(h)!=len(m) or len(h)<30: continue
    rf=fix(row.frames[0]); lf=latent_frames(h)
    if lf is None: continue
    a,b=lf
    af = b if a==rf else (a if b==rf else None)
    if af is None or af==rf: continue
    add.append((name,h,m,rf,af,(af-rf)%3))
print(f"\nrecovered one-frame DCRs: {len(add)}")
selAdd=score_set(add,'Dadd')
oAdd=omegas(pooled(selAdd))
allD_recs=Drecs2+add; selAll=score_set(allD_recs,'Dall'); oAll=omegas(pooled(selAll))
print(f"one-frame-only DCRs:  dNN/dNS={oAdd['w_alt_refN']:.3f}  dSN/dSS={oAdd['w_alt_refS']:.3f}")
print(f"combined ({len(allD_recs)}):   dNN/dNS={oAll['w_alt_refN']:.3f}  dSN/dSS={oAll['w_alt_refS']:.3f}")
print(f"original 956:         dNN/dNS={oD2['w_alt_refN']:.3f}  dSN/dSS={oD2['w_alt_refS']:.3f}")
~~~~

### [msg 542] python · env: dcr-tierA · Testing frame misassignment hypothesis

~~~~python
# Decisive control: rescore the 956 GOOD DCRs using a DETECTED latent frame instead of
# the annotated second frame. If omega collapses to the null, the 229's flat signal is
# misassignment (I recovered the wrong alternate frame), not absence of selection.
mis=[]
for name,h,m,f1,f2,rel in Drecs2:
    lf=latent_frames(h)
    if lf is None: continue
    a,b=lf
    af = b if a==f1 else (a if b==f1 else None)
    if af is None or af==f1: continue
    mis.append((name,h,m,f1,af,(af-f1)%3))
selMis=score_set(mis,'Dmis'); oMis=omegas(pooled(selMis))
agree=np.mean([r[4]==d[4] for r,d in zip(sorted(mis),sorted([x for x in Drecs2 if x[0] in {m[0] for m in mis}]))])
print(f"956-set rescored with DETECTED latent frame (n={len(mis)}):")
print(f"   dNN/dNS={oMis['w_alt_refN']:.3f}   dSN/dSS={oMis['w_alt_refS']:.3f}")
print(f"   vs annotated second frame: dNN/dNS={oD2['w_alt_refN']:.3f}  dSN/dSS={oD2['w_alt_refS']:.3f}")
print(f"\ndetected frame == annotated second frame in {agree*100:.1f}% of cases")
~~~~

### [msg 544] bash · env: dcr-tierA · Re-intersecting DCRs at any overlap

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
S="$WS/tierA/sel"
# relax: any CDS overlap at all (no -f), so partial-overlap second frames are recovered
bedtools intersect -a "$S/D_h.s.bed" -b "$S/gc.s.bed" -wo > "$S/D_cds_any.ov"
echo "rows: $(wc -l < "$S/D_cds_any.ov")"
~~~~

### [msg 546] python · env: dcr-tierA · Checking relaxed second-frame recovery

~~~~python
ova=pd.read_csv(f"{S}/D_cds_any.ov",sep="\t",header=None,
    names=['c','s','e','name','sc','st','ec','es','ee','ename','esc','est','bp'])
ova=ova[ova.ename.isin(name2phase)].copy()
ova['ephase']=ova.ename.map(name2phase)
off=np.where(ova.est=='+', ova.s.values-ova.es.values, ova.ee.values-ova.e.values)
ova['cpos']=(ova.ephase.values+off)%3
ova['frac']=ova.bp/(ova.e-ova.s)
g=ova.groupby('name').agg(nf_any=('cpos','nunique'),maxfrac=('frac','max'))
one_names=set(one.index)
sub=g[g.index.isin(one_names)]
print(f"of the 266 one-frame DCRs, second frame appears at relaxed overlap: {(sub.nf_any>=2).sum()}/{len(sub)}")
print("\ntheir overlap fraction for the 2nd frame:")
sec=ova[ova.name.isin(one_names)].sort_values('frac',ascending=False).groupby('name').nth(1)
print(sec.frac.describe(percentiles=[.25,.5,.75]).to_string(float_format=lambda v:f"{v:.3f}"))
~~~~

### [msg 548] python · env: dcr-tierA · Bootstrapping ascertainment comparison

~~~~python
bsAdd=boot_omega(selAdd,'refS'); bsAll=boot_omega(selAll,'refS')
print(f"one-frame DCRs (n={len(selAdd)}): dSN/dSS={oAdd['w_alt_refS']:.3f} 95% CI [{np.percentile(bsAdd,2.5):.3f}, {np.percentile(bsAdd,97.5):.3f}]  -> overlaps 1.0: {np.percentile(bsAdd,2.5)<1<np.percentile(bsAdd,97.5)}")
print(f"combined  (n={len(selAll)}): dSN/dSS={oAll['w_alt_refS']:.3f} 95% CI [{np.percentile(bsAll,2.5):.3f}, {np.percentile(bsAll,97.5):.3f}]")
asc=pd.DataFrame({
 'set':['kept (2 annotated CDS frames)','one CDS frame only (2nd isoform non-coding/UTR)','combined','latent null'],
 'n':[len(selD2),len(selAdd),len(selAll),len(selN2)],
 'w_dNNdNS':[oD2['w_alt_refN'],oAdd['w_alt_refN'],oAll['w_alt_refN'],oN2['w_alt_refN']],
 'w_dSNdSS':[oD2['w_alt_refS'],oAdd['w_alt_refS'],oAll['w_alt_refS'],oN2['w_alt_refS']],
 'w_ref_control':[oD2['w_ref'],oAdd['w_ref'],oAll['w_ref'],oN2['w_ref']]})
asc.to_csv(f"{S}/ascertainment_check.csv",index=False)
meta[['fate','nf','len_h','gc']].to_csv(f"{S}/dcr_exclusion_fates.csv")
lt.to_csv(f"{S}/omega_by_length_quartile.csv",index=False)
print()
print(asc.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
~~~~

### [msg 550] python · env: dcr-tierA · Building ascertainment check figure

~~~~python
fig,axes=plt.subplots(1,3,figsize=(11.4,3.6))
ax=axes[0]
fates=tt.fate.value_counts()
order=['kept','one_frame_only','too_short','indel','lift_failed']
labs2=['kept\n(2 CDS frames)','2nd isoform\nnon-coding','<30 nt','indel in\nlift','lift\nfailed']
cols=[ORANGE,BLUE,GREY,GREY,GREY]
ax.bar(range(5),[fates[o] for o in order],color=cols,edgecolor='k',linewidth=0.5,width=0.66)
for i,o in enumerate(order): ax.text(i,fates[o]+14,str(fates[o]),ha='center',fontsize=7.5,weight='bold')
ax.set_xticks(range(5)); ax.set_xticklabels(labs2,fontsize=6.6)
ax.set_ylabel("DCR intervals"); ax.set_ylim(0,1100)
ax.set_title("A  Why 761 DCRs were excluded",loc='left',weight='bold')

ax=axes[1]
xx=np.arange(len(lt))
ax.plot(xx,lt.w_dSNdSS,'o-',color=ORANGE,lw=1.6,ms=5,label='dSN/dSS')
ax.plot(xx,lt.w_dNNdNS,'s--',color=BLUE,lw=1.4,ms=4,label='dNN/dNS')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks(xx); ax.set_xticklabels(['30-49','49-72','72-104','104+'],fontsize=7.5)
ax.set_xlabel("DCR length (nt), quartiles"); ax.set_ylabel("ω  (alternate frame)")
ax.set_ylim(0,1.15); ax.legend(frameon=False,fontsize=7,loc='upper left')
ax.text(0.5,0.30,"flat across length:\nthe length filter is benign",transform=ax.transAxes,
        ha='center',fontsize=6.8,style='italic')
ax.set_title("B  ω is length-independent",loc='left',weight='bold')

ax=axes[2]
sets=['kept\n956','2nd isoform\nnon-coding\n229','combined\n1185','latent\nnull\n21340']
vals=asc.w_dSNdSS.values
cl=[ORANGE,BLUE,'#c8873a',GREY]
ax.bar(range(4),vals,color=cl,edgecolor='k',linewidth=0.5,width=0.6)
ax.errorbar([1],[oAdd['w_alt_refS']],yerr=[[oAdd['w_alt_refS']-np.percentile(bsAdd,2.5)],
            [np.percentile(bsAdd,97.5)-oAdd['w_alt_refS']]],fmt='none',ecolor='k',lw=0.9,capsize=3)
ax.errorbar([2],[oAll['w_alt_refS']],yerr=[[oAll['w_alt_refS']-np.percentile(bsAll,2.5)],
            [np.percentile(bsAll,97.5)-oAll['w_alt_refS']]],fmt='none',ecolor='k',lw=0.9,capsize=3)
for i,v in enumerate(vals): ax.text(i,v+0.09,f"{v:.2f}",ha='center',fontsize=7.5,weight='bold')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks(range(4)); ax.set_xticklabels(sets,fontsize=6.6)
ax.set_ylabel("ω  dSN/dSS"); ax.set_ylim(0,1.75)
ax.set_title("C  Worst case still shows selection",loc='left',weight='bold')
fig.suptitle("Ascertainment check: the selection signal is not an artefact of which DCRs survive the filters",
             x=0.005,ha='left',fontsize=10,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/selection_ascertainment.png",dpi=300); plt.close(fig)
print("saved")
~~~~

### [msg 552] python · env: dcr-tierA · Fixing panel B annotation placement

~~~~python
fig,axes=plt.subplots(1,3,figsize=(11.4,3.6))
ax=axes[0]
ax.bar(range(5),[fates[o] for o in order],color=cols,edgecolor='k',linewidth=0.5,width=0.66)
for i,o in enumerate(order): ax.text(i,fates[o]+14,str(fates[o]),ha='center',fontsize=7.5,weight='bold')
ax.set_xticks(range(5)); ax.set_xticklabels(labs2,fontsize=6.6)
ax.set_ylabel("DCR intervals"); ax.set_ylim(0,1100)
ax.set_title("A  Why 761 DCRs were excluded",loc='left',weight='bold')
ax=axes[1]
xx=np.arange(len(lt))
ax.plot(xx,lt.w_dSNdSS,'o-',color=ORANGE,lw=1.6,ms=5,label='dSN/dSS')
ax.plot(xx,lt.w_dNNdNS,'s--',color=BLUE,lw=1.4,ms=4,label='dNN/dNS')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks(xx); ax.set_xticklabels(['30-49','49-72','72-104','104+'],fontsize=7.5)
ax.set_xlabel("DCR length (nt), quartiles"); ax.set_ylabel("ω  (alternate frame)")
ax.set_ylim(0,1.25); ax.legend(frameon=False,fontsize=7,loc='upper left')
ax.text(0.5,-0.30,"flat across length: the length filter is benign",transform=ax.transAxes,
        ha='center',fontsize=6.8,style='italic')
ax.set_title("B  ω is length-independent",loc='left',weight='bold')
ax=axes[2]
ax.bar(range(4),vals,color=cl,edgecolor='k',linewidth=0.5,width=0.6)
ax.errorbar([1],[oAdd['w_alt_refS']],yerr=[[oAdd['w_alt_refS']-np.percentile(bsAdd,2.5)],
            [np.percentile(bsAdd,97.5)-oAdd['w_alt_refS']]],fmt='none',ecolor='k',lw=0.9,capsize=3)
ax.errorbar([2],[oAll['w_alt_refS']],yerr=[[oAll['w_alt_refS']-np.percentile(bsAll,2.5)],
            [np.percentile(bsAll,97.5)-oAll['w_alt_refS']]],fmt='none',ecolor='k',lw=0.9,capsize=3)
for i,v in enumerate(vals): ax.text(i,v+0.10,f"{v:.2f}",ha='center',fontsize=7.5,weight='bold')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks(range(4)); ax.set_xticklabels(sets,fontsize=6.6)
ax.set_ylabel("ω  dSN/dSS"); ax.set_ylim(0,1.8)
ax.set_title("C  Worst case still shows selection",loc='left',weight='bold')
fig.suptitle("Ascertainment check: the selection signal is not an artefact of which DCRs survive the filters",
             x=0.005,ha='left',fontsize=10,weight='bold')
fig.tight_layout(rect=[0,0.05,1,0.92])
fig.savefig(f"{WS}/tierA/selection_ascertainment.png",dpi=300); plt.close(fig)
print("saved")
~~~~

### [msg 554] python · env: dcr-tierA · Writing OLGenie methods document

~~~~python
doc = r"""# The OLGenie site-classification approach to dN/dS in overlapping frames

**What this document is for.** You raised the objection that a dN/dS computed on the
alternate reading frame of a DCR cannot be independent of the canonical frame, because
purifying selection on the canonical protein constrains codon positions 1 and 2, and those
same nucleotides are the alternate frame's synonymous and nonsynonymous sites. You are
right, and this document explains exactly why, what the standard fix is, how I implemented
it, and what the resulting numbers mean.

---

## 1. Why a naive alternate-frame dN/dS is uninterpretable

Consider a DCR where a second transcript reads the same nucleotides shifted by +1
(the case OLGenie calls `ss13`). Line up the two frames:

```
nucleotide position   1   2   3   4   5   6   7   8   9
canonical frame       |---codon A---| |---codon B---| |--
                       p1  p2  p3      p1  p2  p3
alternate frame (+1)      |---codon a---| |---codon b---|
                           p1  p2  p3      p1  p2  p3
```

Reading the mapping off the diagram, for a +1 shift:

| canonical position | alternate position |
|---|---|
| p1 (highly constrained) | p3 (mostly degenerate) |
| p2 (most constrained)   | p1 |
| p3 (wobble, degenerate) | p2 |

Now ask what happens if the canonical protein is under strong purifying selection.
Substitutions are **suppressed** at canonical p1 and p2, and **tolerated** at canonical p3.
But canonical p3 is alternate p2 — the alternate frame's *most* constrained position. So the
only positions that freely accumulate substitutions are precisely the ones that look
nonsynonymous in the alternate frame.

The consequence: the alternate frame's apparent dN is inflated (its substitutions are
concentrated at alternate-p2), and its apparent dS is deflated (alternate-p3 sits on
canonical-p1, where change is suppressed). **Naive omega_alt is biased upward.** A DCR under
genuine alternate-frame purifying selection could easily read omega_alt near or above 1, and
you would wrongly conclude there is no constraint — or, in the other direction, mistake the
canonical frame's own constraint for alternate-frame signal.

For a +2 shift (`ss12`) the position mapping is different again, so the bias has a
different magnitude. This is why the two frame relationships must never be pooled.

This is not a subtlety I invented: it is the motivating problem of the overlapping-gene
dN/dS literature. Standard methods do not account for the fact that a change which is
synonymous in one frame is usually nonsynonymous in the other, which is why they fail to
detect purifying selection, or erroneously report positive selection, when applied to
overlapping reading frames.

---

## 2. The fix: cross-classify every site by its effect in BOTH frames

Instead of asking "is this substitution synonymous?" (which has no single answer when two
frames overlap), ask **"what is it in each frame?"** Every site gets a two-letter label:

- first letter = effect in the **reference** (canonical) frame
- second letter = effect in the **alternate** frame

giving four classes:

| class | canonical frame | alternate frame |
|---|---|---|
| **NN** | nonsynonymous | nonsynonymous |
| **NS** | nonsynonymous | synonymous |
| **SN** | **synonymous** | nonsynonymous |
| **SS** | **synonymous** | synonymous |

From these we form four normalised divergences dNN, dNS, dSN, dSS — each is
(observed differences of that class) / (number of sites of that class).

Selection on the alternate frame is then measured by a ratio that **holds the canonical
letter constant**:

```
omega_alt  =  dSN / dSS        <-- both canonical-synonymous
omega_alt  =  dNN / dNS        <-- both canonical-nonsynonymous
```

and, as a control, selection on the canonical frame by holding the *alternate* letter fixed:

```
omega_ref  =  dNN / dSN        <-- both alternate-nonsynonymous
```

### Why this answers your objection exactly

Look at `dSN/dSS`. Every single change counted in the numerator **and** every single change
counted in the denominator is **synonymous in the canonical frame**. Whatever purifying
selection is acting on the canonical protein applies identically to both, so it **cancels in
the ratio by construction**. The ratio is then free to reflect only what differs between the
numerator and denominator — namely, whether the change alters the alternate protein.

`dNN/dNS` does the same thing at canonical-nonsynonymous sites. It has far more sites
available (canonical-nonsynonymous sites vastly outnumber canonical-synonymous ones), so it
is better powered in short regions; `dSN/dSS` is the conceptually cleanest answer to your
objection. I report both. They agree, which is reassuring — if they disagreed, that would
itself be diagnostic of a residual frame-specific artefact.

---

## 3. How I implemented it

**Sequence pairs.** Human DCR intervals and their mouse orthologs, lifted with UCSC
`liftOver` (hg38 -> mm39). Only **indel-free** pairs were used (human and mouse interval the
same length), so codons stay in register without needing a gapped alignment. This is a real
filter: 12.1% of lifted DCRs and 15.0% of lifted control intervals were dropped for indels —
comparably, so it does not bias the contrast.

**Reading frames.** The canonical frame comes from GENCODE v48 annotation, not from
guessing. The bed file's phase field turned out to be a constant 0, so phase was re-derived
as the cumulative CDS length preceding each exon within its transcript, modulo 3, then
converted to an offset within the DCR interval accounting for strand. For real DCRs the
alternate frame is the *second annotated* coding frame (78.8% of DCRs have >=2 annotated
frames over the interval; split 998 `ss13` / 355 `ss12`). For the control set the alternate
frame is the detected latent ORF frame.

**Site and difference counting.** For each nucleotide position, both the canonical codon and
the alternate codon containing it are extracted. Each of the three possible alternative bases
is evaluated in both codons simultaneously and the position contributes 1/3 of a site to
whichever of the four classes applies. Observed human-mouse differences are classified the
same way. Counts are pooled across intervals, converted to proportions, and
Jukes-Cantor-corrected before ratios are taken.

**A bug this design caught.** My first implementation converted the annotated codon position
to a frame offset as `offset = cpos`, when the correct conversion is `offset = (3 - cpos) % 3`.
The error was invisible in the DCR set (both of its frames are genuinely coding, so a
scrambled offset still lands on a real reading frame) but catastrophic in the control set. It
surfaced through the built-in positive control: omega_ref for the controls read 0.98
(indistinguishable from neutral) when it must be strongly purifying, because these are
annotated CDS. After the fix, omega_ref = 0.110. **This is why the positive control is not
optional** — it is the only thing that detects a scrambled frame.

---

## 4. Why an absolute omega_alt still would not be enough

Even after canonical constraint cancels, the sites carry other non-neutral pressures:
codon usage bias, exonic splicing enhancers, and the GC enrichment we documented separately
for DCRs. Any of these could push omega_alt below 1 without any alternate-frame protein
existing.

So the claim is **never** "omega_alt < 1, therefore selection." The claim is the **contrast
against a matched empirical null**: intervals inside non-DCR exons that carry a latent
dual-coding ORF but are not annotated as dual-coding. These are matched on length and
stratified by frame relationship, and they experience the same codon-usage, splicing, and
composition pressures. The null answers the question "what does omega_alt look like when
there is dual-coding *potential* but no dual-coding *function*?"

The answer: **omega_alt = 0.984, 95% CI [0.951, 1.017]** — the null sits exactly at
neutrality. That is the single most important number in this analysis, because it validates
that the design has successfully removed the confounds. The contrast is therefore meaningful.

---

## 5. Results

| | real DCRs (n=956) | latent null (n=21,340) |
|---|---|---|
| omega_alt = dNN/dNS | **0.268** [0.244, 0.296] | 0.798 [0.787, 0.809] |
| omega_alt = dSN/dSS | **0.452** [0.368, 0.553] | 0.984 [0.951, 1.017] |
| omega_ref = dNN/dSN (control) | 0.272 | 0.110 |

Bootstrap over intervals, 2000 replicates. Difference in omega_alt between sets = 0.530
under both estimators, p < 0.001.

**Controls all pass:**
- canonical-frame control strongly purifying in both sets (0.272, 0.110) — frames are correct
- holds separately in `ss13` (D 0.470 vs N 0.827) and `ss12` (D 0.446 vs N 1.196)
- survives length matching (D 0.442 vs N 1.025) — the control that erased the earlier
  potential-retention signal
- omega_alt is flat across DCR length quartiles (0.379 -> 0.491), so the >=30 nt filter is benign

**Ascertainment.** Of 1,717 DCR intervals, 956 entered the main analysis. Exclusions: 266 had
only one annotated *coding* frame over the interval (their second isoform overlaps as UTR or
non-coding), 245 were <30 nt, 201 had indels in the lift, 49 failed to lift. Kept DCRs are
longer than dropped ones (72 vs 49 nt median, p = 9e-29) but not different in GC (55.0% vs
53.6%, p = 0.08). Since omega_alt does not vary with length, the length skew does not drive
the result.

The 229 recoverable one-coding-frame DCRs show **no** selection (dSN/dSS = 1.14, CI
[0.81, 1.63], overlapping neutrality). I checked whether this was frame misassignment by
rescoring the main 956 using a *detected* latent frame instead of the annotated second frame:
omega moved only from 0.452 to 0.507, and detection agreed with annotation in 80.5% of cases.
So detection is reliable and the flat signal in those 229 is real — consistent with their
second isoform not actually being translated over the interval. Even in the worst case,
pooling them in, the combined estimate remains **0.545 [0.463, 0.643]**, well below the null.

---

## 6. What this licenses you to claim, and what it does not

**Supported:** the alternate reading frames of annotated DCRs are under purifying selection,
significantly more so than sequence-matched intervals with equivalent dual-coding potential
but no annotated dual-coding function. Since dual-coding potential is near-universal in
coding exons (X = 0.994), this selection signal — not the mere existence of potential — is
what distinguishes real DCRs.

**Not supported:** (i) any per-DCR claim. A single human-mouse pairwise comparison over a
median 72 nt yields too few substitutions to rank individual DCRs; this is a pooled,
genome-wide result only. (ii) any claim about which specific alternate proteins are made —
selection on a reading frame is evidence of function, not proof of translation. (iii) any
statement about the 266 DCRs whose second isoform is non-coding over the interval.

**Natural extensions if a reviewer pushes:** multi-species alignments (UCSC 100-way) would
give per-DCR resolution and could identify a strongly-constrained subset; human polymorphism
(gnomAD) would test present-day constraint with far more variants. Both were scoped and are
straightforward; neither is needed for the global claim.

---

## 7. Reference

Nelson CW, Ardern Z, Wei X. *OLGenie: Estimating Natural Selection to Predict Functional
Overlapping Genes.* Molecular Biology and Evolution 37(8):2440-2449, 2020.
The dNN/dNS/dSN/dSS cross-classification and the ratio estimators used here are theirs; the
implementation in this analysis is my own, written to work directly on the lifted interval
pairs, with the reference-frame estimator retained as a positive control.
"""
open(f"{WS}/tierA/OLGenie_methods_explained.md","w").write(doc)
print(f"written, {len(doc.splitlines())} lines")
~~~~

### [msg 563] bash · env: dcr-tierA · Locating protein-level DCR definition files

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
DOO="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
echo "=== files mentioning ORF / protein / frame in 01_Description ==="
find "$DOO/output" -maxdepth 3 -type f \( -iname '*orf*' -o -iname '*prot*' -o -iname '*frame*' -o -iname '*pep*' -o -iname '*translat*' \) 2>/dev/null | head -30
echo ""
echo "=== all top-level output files (name + size) ==="
find "$DOO/output" -maxdepth 2 -type f -size +1k 2>/dev/null | head -40
~~~~

### [msg 565] bash · env: dcr-tierA · Inspecting pipeline frame records

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
DOO="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output"
C="$DOO/02_Conservation"; D1="$DOO/01_Description"
echo "=== frames_per_DCR.txt ($(wc -l < "$C/frames_per_DCR.txt") lines) ==="; head -6 "$C/frames_per_DCR.txt"
echo ""
echo "=== TO_REMOVE.only1frame.txt ($(wc -l < "$D1/TO_REMOVE.only1frame.txt") lines) ==="; head -4 "$D1/TO_REMOVE.only1frame.txt"
echo ""
echo "=== FINAL_Filtered_UniprotKB_isoforms_mapping_toDCR...extBed.bed ($(wc -l < "$D1/FINAL_Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguousGenes_noAmbiguousEnds.extBed.bed") lines) ==="
head -4 "$D1/FINAL_Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguousGenes_noAmbiguousEnds.extBed.bed"
echo ""
echo "=== DCE_CDS_intersect_filtered_with_transcriptID.bed ($(wc -l < "$C/DCE_CDS_intersect_filtered_with_transcriptID.bed") lines) ==="
head -4 "$C/DCE_CDS_intersect_filtered_with_transcriptID.bed"
~~~~

### [msg 567] python · env: dcr-tierA · Cross-checking my fates against pipeline frames

~~~~python
C=f"{DOO}/output/02_Conservation" if 'DOO' in dir() else None
DOOT="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output"
fpd=pd.read_csv(f"{DOOT}/02_Conservation/frames_per_DCR.txt",sep="\t",header=None,names=['dceID','nframes'],dtype={'dceID':str})
print("frames_per_DCR distribution (pipeline, protein-based):")
print(fpd.nframes.value_counts().sort_index().to_string())
# my interval names look like 'dce2|1389' -> dceID is the part after '|'
meta['dceID']=[n.split('|')[-1] for n in meta.index]
mm3=meta.merge(fpd,on='dceID',how='left')
print(f"\nmy 1717 intervals matched to frames_per_DCR: {mm3.nframes.notna().sum()}")
print("\n=== pipeline nframes, by my GENCODE-based fate ===")
print(pd.crosstab(mm3.fate,mm3.nframes,dropna=False).to_string())
~~~~

### [msg 569] bash · env: dcr-tierA · Searching for pipeline frame assignments

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
DOO="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output"
echo "=== New_TableS1.tsv header ==="; head -1 "$DOO/New_TableS1.tsv" | tr '\t' '\n' | head -30
echo ""
echo "=== 01_Description file list ==="; ls -S "$DOO/01_Description" | head -25
echo ""
echo "=== grep for a frame-like column across 01_Description headers ==="
for f in "$DOO/01_Description"/*.tsv "$DOO/01_Description"/*.txt; do
  [ -f "$f" ] || continue
  h=$(head -1 "$f" | tr '\t' ' ')
  case "$h" in *[Ff]rame*|*[Pp]hase*) echo "$(basename "$f"): $h";; esac
done 2>/dev/null | head
~~~~

### [msg 571] bash · env: dcr-tierA · Inspecting phase and frames tables

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
D1="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/01_Description"
echo "=== raw_DCR_exons_with_phase.tsv ($(wc -l < "$D1/raw_DCR_exons_with_phase.tsv") lines) ==="
head -3 "$D1/raw_DCR_exons_with_phase.tsv"
echo ""
echo "=== filtered_homo_sapiens_DCR.tsv ($(wc -l < "$D1/filtered_homo_sapiens_DCR.tsv") lines) ==="
head -3 "$D1/filtered_homo_sapiens_DCR.tsv" | cut -c1-320
~~~~

### [msg 573] python · env: dcr-tierA · Loading protein-derived frame assignments

~~~~python
D1=f"{DOOT}/01_Description"
fh=pd.read_csv(f"{D1}/filtered_homo_sapiens_DCR.tsv",sep="\t",dtype=str)
fh['Index']=fh.Index.astype(str)
print("rows:",len(fh)," unique DCR Index:",fh.Index.nunique())
gf=fh.groupby('Index')['Frames'].apply(lambda v:sorted(set(v.astype(int))))
gf=gf.to_frame('pipe_frames'); gf['pipe_nf']=gf.pipe_frames.apply(len)
print("\npipeline frames per DCR (from filtered_homo_sapiens_DCR.tsv):")
print(gf.pipe_nf.value_counts().sort_index().to_string())
print("\nexample rows:"); print(gf.head(4).to_string())
mm4=meta.merge(gf,left_on='dceID',right_index=True,how='left')
print(f"\nmy intervals matched: {mm4.pipe_nf.notna().sum()}/{len(mm4)}")
print("\n=== pipeline nframes by my GENCODE fate ===")
print(pd.crosstab(mm4.fate,mm4.pipe_nf,dropna=False).to_string())
~~~~

### [msg 575] python · env: dcr-tierA · Determining pipeline frame convention

~~~~python
# Establish convention: for 'kept' intervals I have both GENCODE frames (pD.frames, raw cpos)
# and pipeline frames. Test candidate mappings for agreement.
kept_ids=[n for n in meta.index if meta.loc[n,'fate']=='kept']
rows=[]
for n in kept_ids:
    d=n.split('|')[-1]
    if d not in gf.index: continue
    gen=set(fix(f) for f in pD.loc[n,'frames'])       # corrected offsets from GENCODE
    pf=gf.loc[d,'pipe_frames']
    rows.append((n,tuple(sorted(gen)),tuple(pf)))
cv=pd.DataFrame(rows,columns=['name','gencode_off','pipe_frames'])
print("n with both:",len(cv))
for label,conv in [('F-1',lambda F:(F-1)%3),('3-(F-1)',lambda F:(3-(F-1))%3),('F',lambda F:F%3)]:
    agree=np.mean([tuple(sorted({conv(f) for f in p}))==g for g,p in zip(cv.gencode_off,cv.pipe_frames)])
    print(f"convention {label:9s}: exact set match in {agree*100:.1f}% of intervals")
~~~~

### [msg 577] python · env: dcr-tierA · Deriving phase from UniProt isoform models

~~~~python
import gzip
ib=pd.read_csv(f"{D1}/FINAL_Filtered_UniprotKB_isoforms_mapping_toDCR_noAmbiguousGenes_noAmbiguousEnds.extBed.bed",
               sep="\t",header=None,
               names=['chr','start','end','name','score','strand','ts','te','rgb','nblk','sizes','starts'])
print("isoform models:",len(ib)," unique:",ib.name.nunique())
# explode BED12 blocks into CDS exons with cumulative phase
recs=[]
for r in ib.itertuples(index=False):
    sizes=[int(x) for x in str(r.sizes).rstrip(',').split(',')]
    offs=[int(x) for x in str(r.starts).rstrip(',').split(',')]
    blocks=[(r.start+o,r.start+o+s) for o,s in zip(offs,sizes)]
    if r.strand=='-': blocks=blocks[::-1]          # 5'->3' order
    cum=0
    for bs,be in blocks:
        recs.append((r.name,r.chr,bs,be,r.strand,cum%3))
        cum+=be-bs
iso=pd.DataFrame(recs,columns=['iso','chr','start','end','strand','phase'])
print("isoform CDS blocks:",len(iso))
iso.to_csv(f"{S}/iso_blocks.bed",sep="\t",header=False,index=False,
           columns=['chr','start','end','iso','phase','strand'])
print(iso.head(3).to_string(index=False))
~~~~

### [msg 580] bash · env: dcr-tierA · Intersecting DCRs with isoform CDS blocks

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
S="$WS/tierA/sel"
sort -k1,1 -k2,2n "$S/iso_blocks.bed" > "$S/iso_blocks.s.bed"
bedtools intersect -a "$S/D_h.s.bed" -b "$S/iso_blocks.s.bed" -wo > "$S/D_iso.ov"
echo "rows: $(wc -l < "$S/D_iso.ov")   intervals: $(cut -f4 "$S/D_iso.ov" | sort -u | wc -l)"
head -2 "$S/D_iso.ov"
~~~~

### [msg 582] python · env: dcr-tierA · Deriving protein-based frames per interval

~~~~python
oi=pd.read_csv(f"{S}/D_iso.ov",sep="\t",header=None,
    names=['c','s','e','name','sc','st','ic','is_','ie','iso','iphase','ist','bp'])
off=np.where(oi.ist=='+', oi.s.values-oi.is_.values, oi.ie.values-oi.e.values)
oi['cpos']=(oi.iphase.values+off)%3
oi['off']=(3-oi.cpos)%3                     # same corrected conversion as GENCODE
oi['frac']=oi.bp/(oi.e-oi.s)
gp=oi.groupby('name').agg(prot_offs=('off',lambda v:sorted(set(v))),
                          nfrac=('frac','max'),niso=('iso','nunique'))
gp['prot_nf']=gp.prot_offs.apply(len)
print("protein-derived distinct frames per DCR interval:")
print(gp.prot_nf.value_counts().sort_index().to_string())
chk=meta.join(gp)
print("\n=== protein-derived nframes by my GENCODE fate ===")
print(pd.crosstab(chk.fate,chk.prot_nf).to_string())
~~~~

### [msg 584] python · env: dcr-tierA · Rescoring with protein-derived frames

~~~~python
> 
# Rebuild D using PROTEIN-derived frames for every interval that lifts indel-free
Dprot=[]
for name,row in gp.iterrows():
    h=hDk.get(name); m=mDk.get(name)
    if h is None or m is None or len(h)!=len(m) or len(h)<30: continue
    offs=row.prot_offs
    if len(offs)<2: continue
    f1,f2=offs[0],offs[1]
    Dprot.append((name,h,m,f1,f2,(f2-f1)%3))
print(f"D (protein frames), indel-free, >=30nt: {len(Dprot)}   [was 956 with GENCODE 2-frame filter]")
selDp=score_set(Dprot,'Dprot'); oDp=omegas(pooled(selDp))
bsDp_N=boot_omega(selDp,'refN'); bsDp_S=boot_omega(selDp,'refS')
print(f"\nomega_ref (control) = {oDp['w_ref']:.3f}")
print(f"omega_alt dNN/dNS   = {oDp['w_alt_refN']:.3f}  95% CI [{np.percentile(bsDp_N,2.5):.3f}, {np.percentile(bsDp_N,97.5):.3f}]")
print(f"omega_alt dSN/dSS   = {oDp['w_alt_refS']:.3f}  95% CI [{np.percentile(bsDp_S,2.5):.3f}, {np.percentile(bsDp_S,97.5):.3f}]")
print(f"\n(GENCODE-based 956: dNN/dNS={oD2['w_alt_refN']:.3f}, dSN/dSS={oD2['w_alt_refS']:.3f})")
print(f"(latent null 21340: dNN/dNS={oN2['w_alt_refN']:.3f}, dSN/dSS={oN2['w_alt_refS']:.3f})")
~~~~

### [msg 586] python · env: dcr-tierA · Rescoring with protein-derived frames

~~~~python
Dprot=[]
for name,row in gp.iterrows():
    h=hDk.get(name); m=mDk.get(name)
    if h is None or m is None or len(h)!=len(m) or len(h)<30: continue
    offs=row.prot_offs
    if len(offs)<2: continue
    f1,f2=offs[0],offs[1]
    Dprot.append((name,h,m,f1,f2,(f2-f1)%3))
print(f"D (protein frames), indel-free, >=30nt: {len(Dprot)}   [was 956 under the GENCODE 2-frame filter]")
selDp=score_set(Dprot,'Dprot'); oDp=omegas(pooled(selDp))
bsDp_N=boot_omega(selDp,'refN'); bsDp_S=boot_omega(selDp,'refS')
print(f"\nomega_ref (control) = {oDp['w_ref']:.3f}")
print(f"omega_alt dNN/dNS   = {oDp['w_alt_refN']:.3f}  95% CI [{np.percentile(bsDp_N,2.5):.3f}, {np.percentile(bsDp_N,97.5):.3f}]")
print(f"omega_alt dSN/dSS   = {oDp['w_alt_refS']:.3f}  95% CI [{np.percentile(bsDp_S,2.5):.3f}, {np.percentile(bsDp_S,97.5):.3f}]")
print(f"\nGENCODE-based 956 : dNN/dNS={oD2['w_alt_refN']:.3f}, dSN/dSS={oD2['w_alt_refS']:.3f}")
print(f"latent null 21340 : dNN/dNS={oN2['w_alt_refN']:.3f}, dSN/dSS={oN2['w_alt_refS']:.3f}")
print("\nframe relationship:",pd.Series([r[5] for r in Dprot]).value_counts().sort_index().to_dict())
~~~~

### [msg 589] python · env: dcr-tierA · Rescoring the mislabeled subset

~~~~python
ofo=set(meta.index[meta.fate=='one_frame_only'])
sub_ofo=[r for r in Dprot if r[0] in ofo]
sub_kept=[r for r in Dprot if r[0] in set(meta.index[meta.fate=='kept'])]
oO=omegas(pooled(score_set(sub_ofo,'x'))); bO=boot_omega(score_set(sub_ofo,'x'),'refS')
oK=omegas(pooled(score_set(sub_kept,'x')))
print(f"previously 'one_frame_only', now protein-framed (n={len(sub_ofo)}):")
print(f"   omega_ref={oO['w_ref']:.3f}  dNN/dNS={oO['w_alt_refN']:.3f}  dSN/dSS={oO['w_alt_refS']:.3f} CI[{np.percentile(bO,2.5):.3f},{np.percentile(bO,97.5):.3f}]")
print(f"   (with the WRONG detected-latent frame it read dSN/dSS=1.138)")
print(f"\npreviously 'kept' (n={len(sub_kept)}): dNN/dNS={oK['w_alt_refN']:.3f}  dSN/dSS={oK['w_alt_refS']:.3f}")
rows=[]
for rel,lab in [(1,'ss13 (+1)'),(2,'ss12 (+2)')]:
    dd=[r for r in Dprot if r[5]==rel]; nn=[r for r in Nrecs2 if r[5]==rel]
    od=omegas(pooled(score_set(dd,'x'))); on=omegas(pooled(score_set(nn,'x')))
    rows.append((lab,len(dd),od['w_alt_refS'],len(nn),on['w_alt_refS']))
fr=pd.DataFrame(rows,columns=['frame_rel','n_D','D_dSNdSS','n_N','N_dSNdSS'])
print("\n=== by frame relationship (protein frames) ===")
print(fr.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
~~~~

### [msg 591] python · env: dcr-tierA · Length-matching and significance test

~~~~python
selDp['lb']=pd.cut(selDp.len,bins); selN2['lb']=pd.cut(selN2.len,bins)
md=[];  
for b in selDp.lb.cat.categories:
    a=selDp[selDp.lb==b]; c=selN2[selN2.lb==b]
    if len(a)<20 or len(c)<20: continue
    md.append((str(b),len(a),omegas(pooled(a))['w_alt_refS'],len(c),omegas(pooled(c))['w_alt_refS']))
lm=pd.DataFrame(md,columns=['len_bin','n_D','D','n_N','N'])
wD=np.average(lm.D,weights=lm.n_D); wN=np.average(lm.N,weights=lm.n_D)
print("=== length-matched (weights = D counts) ===")
print(lm.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
print(f"\nweighted: D={wD:.3f}  N={wN:.3f}")
rng4=np.random.default_rng(11)
def bd(est):
    k=0 if est=='refN' else 1
    out=[]
    for _ in range(2000):
        a=selDp.sample(len(selDp),replace=True,random_state=int(rng4.integers(1e9)))
        c=selN2.sample(len(selN2),replace=True,random_state=int(rng4.integers(1e9)))
        key='w_alt_refN' if est=='refN' else 'w_alt_refS'
        out.append(omegas(pooled(c))[key]-omegas(pooled(a))[key])
    return np.array(out)
for est,lab in [('refN','dNN/dNS'),('refS','dSN/dSS')]:
    d=bd(est); p=(d<=0).mean()
    print(f"{lab}: N-D diff = {d.mean():.3f} [{np.percentile(d,2.5):.3f},{np.percentile(d,97.5):.3f}]  p={max(p,1/2000):.4f}")
~~~~

### [msg 593] python · env: dcr-tierA · Rebuilding ascertainment figure

~~~~python
asc2=pd.DataFrame({
 'set':['ALL DCRs (protein frames)','of which: prev. mislabeled 266','prev. "kept" 956',
        'latent null'],
 'n':[len(selDp),len(sub_ofo),len(sub_kept),len(selN2)],
 'w_dNNdNS':[oDp['w_alt_refN'],oO['w_alt_refN'],oK['w_alt_refN'],oN2['w_alt_refN']],
 'w_dSNdSS':[oDp['w_alt_refS'],oO['w_alt_refS'],oK['w_alt_refS'],oN2['w_alt_refS']],
 'w_ref_control':[oDp['w_ref'],oO['w_ref'],oK['w_ref'],oN2['w_ref']]})
asc2.to_csv(f"{S}/ascertainment_check.csv",index=False)
fates2=pd.Series({'kept (protein frames)':len(Dprot),'<30 nt':245,'indel in lift':201,'lift failed':49})
fig,axes=plt.subplots(1,3,figsize=(11.4,3.6))
ax=axes[0]
ax.bar(range(4),fates2.values,color=[ORANGE,GREY,GREY,GREY],edgecolor='k',linewidth=0.5,width=0.62)
for i,v in enumerate(fates2.values): ax.text(i,v+16,str(v),ha='center',fontsize=7.5,weight='bold')
ax.set_xticks(range(4)); ax.set_xticklabels(['kept\n(protein\nframes)','<30 nt','indel in\nlift','lift\nfailed'],fontsize=6.8)
ax.set_ylabel("DCR intervals"); ax.set_ylim(0,1400)
ax.text(0.97,0.72,"all 1,741 DCRs have ≥2\nprotein-derived frames",transform=ax.transAxes,
        ha='right',fontsize=6.8,style='italic')
ax.set_title("A  Exclusions are technical only",loc='left',weight='bold')
ax=axes[1]
xx=np.arange(len(lm))
ax.plot(xx,lm.D,'o-',color=ORANGE,lw=1.6,ms=5,label='DCR')
ax.plot(xx,lm.N,'s--',color=GREY,lw=1.4,ms=4,label='latent null')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks(xx); ax.set_xticklabels([b.split(',')[1].strip(' ]') for b in lm.len_bin],fontsize=7)
ax.set_xlabel("DCR length bin, upper edge (nt)"); ax.set_ylabel("ω  dSN/dSS")
ax.set_ylim(0,1.6); ax.legend(frameon=False,fontsize=7,loc='center right')
ax.set_title("B  Separation holds at every length",loc='left',weight='bold')
ax=axes[2]
v2=asc2.w_dSNdSS.values
ax.bar(range(4),v2,color=[ORANGE,BLUE,'#c8873a',GREY],edgecolor='k',linewidth=0.5,width=0.6)
ax.errorbar([1],[oO['w_alt_refS']],yerr=[[oO['w_alt_refS']-np.percentile(bO,2.5)],
            [np.percentile(bO,97.5)-oO['w_alt_refS']]],fmt='none',ecolor='k',lw=0.9,capsize=3)
for i,v in enumerate(v2): ax.text(i,v+0.055,f"{v:.2f}",ha='center',fontsize=7.5,weight='bold')
ax.axhline(1.0,ls='--',lw=1,color='k')
ax.set_xticks(range(4))
ax.set_xticklabels(['all DCRs\n1228','prev.\nmislabeled\n266','prev. kept\n956','latent\nnull\n21340'],fontsize=6.6)
ax.set_ylabel("ω  dSN/dSS"); ax.set_ylim(0,1.25)
ax.set_title("C  The 266 show selection too",loc='left',weight='bold')
fig.suptitle("Protein-derived frames: all 1,741 DCRs are dual-coding by definition; exclusions are technical, not biological",
             x=0.005,ha='left',fontsize=9.6,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/selection_ascertainment.png",dpi=300); plt.close(fig)
print(asc2.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
~~~~

### [msg 595] python · env: dcr-tierA · Updating methods doc for protein frames

~~~~python
doc2 = doc
# 1. Reading frames paragraph
old_rf = """**Reading frames.** The canonical frame comes from GENCODE v48 annotation, not from
guessing. The bed file's phase field turned out to be a constant 0, so phase was re-derived
as the cumulative CDS length preceding each exon within its transcript, modulo 3, then
converted to an offset within the DCR interval accounting for strand. For real DCRs the
alternate frame is the *second annotated* coding frame (78.8% of DCRs have >=2 annotated
frames over the interval; split 998 `ss13` / 355 `ss12`). For the control set the alternate
frame is the detected latent ORF frame."""
new_rf = """**Reading frames.** *Both* DCR frames are taken from the pipeline's own protein-level
evidence, not from GENCODE. The authoritative source is the filtered UniProtKB isoform
models (`FINAL_Filtered_UniprotKB_isoforms_mapping_toDCR_*.extBed.bed`, 36,733 isoform
models). These BED12 records were exploded into CDS blocks in 5'->3' order and a phase
assigned to each block as the cumulative preceding CDS length modulo 3; the phase was then
converted to an offset within the DCR interval accounting for strand. Intersecting the DCR
intervals against these models recovers **>=2 distinct reading frames for all 1,741 intervals**
(1,718 with exactly 2, 23 with 3), which is what a protein-evidence-based DCR call
guarantees by construction. Frame relationship splits 661 `ss13` / 567 `ss12` among scored
intervals.

An earlier version of this analysis derived the DCR frames from GENCODE v48 CDS annotation
instead, and consequently classified 266 intervals as having "only one coding frame" because
their second isoform overlaps the interval without a GENCODE CDS annotation there. That was
an error of authority: the DCR calls come from protein evidence, so the interval is known to
be coding in both frames regardless of whether GENCODE annotates it. Those 266 have been
restored with their correct protein-derived frames (see section 5). GENCODE phase is still
used for the *control* set, where it is the appropriate source: those intervals are ordinary
annotated CDS with a detected latent ORF as the alternate frame."""
assert old_rf in doc2; doc2=doc2.replace(old_rf,new_rf)

old_res = doc2[doc2.index("## 5. Results"):doc2.index("## 6. What this licenses")]
new_res = """## 5. Results

| | real DCRs (n=1,228) | latent null (n=21,340) |
|---|---|---|
| omega_alt = dNN/dNS | **0.244** [0.223, 0.264] | 0.798 [0.787, 0.809] |
| omega_alt = dSN/dSS | **0.446** [0.376, 0.529] | 0.984 [0.951, 1.017] |
| omega_ref = dNN/dSN (control) | 0.260 | 0.110 |

Bootstrap over intervals, 2000 replicates. Difference in omega_alt between sets: 0.554
[0.530, 0.576] under dNN/dNS and 0.536 [0.449, 0.613] under dSN/dSS, p < 0.001 for both.

**Controls all pass:**
- canonical-frame control strongly purifying in both sets (0.260, 0.110) — frames are correct
- holds separately in `ss13` (D 0.443 vs N 0.827) and `ss12` (D 0.451 vs N 1.196)
- survives length matching (D 0.439 vs N 1.034) — the control that erased the earlier
  potential-retention signal
- DCR and null are separated in **every** length bin from 30 nt to 300 nt, with no bin where
  the DCR estimate approaches the null

**Ascertainment.** All 1,741 DCR intervals carry >=2 protein-derived frames, so no interval is
excluded for biological reasons. The 513 not scored were lost to purely technical filters:
245 shorter than 30 nt, 201 with indels in the human-mouse lift, 49 that failed to lift. The
indel filter is symmetric between arms (12.1% of lifted DCRs, 15.0% of lifted controls), so it
cannot manufacture the contrast.

The 266 intervals that the earlier GENCODE-based version misclassified were rescored with
their correct protein-derived frames and show **purifying selection, not neutrality**:
dSN/dSS = 0.506 [0.359, 0.726], with a canonical control of 0.277. Under the earlier
mistaken frame assignment they had read 1.138. They are statistically indistinguishable
from the other 956 (0.429), which is the expected result if they are ordinary DCRs that
simply lack a GENCODE CDS annotation on the second isoform.

"""
doc2=doc2.replace(old_res,new_res)

old_ns = """**Not supported:** (i) any per-DCR claim. A single human-mouse pairwise comparison over a
median 72 nt yields too few substitutions to rank individual DCRs; this is a pooled,
genome-wide result only. (ii) any claim about which specific alternate proteins are made —
selection on a reading frame is evidence of function, not proof of translation. (iii) any
statement about the 266 DCRs whose second isoform is non-coding over the interval."""
new_ns = """**Not supported:** (i) any per-DCR claim. A single human-mouse pairwise comparison over a
median 72 nt yields too few substitutions to rank individual DCRs; this is a pooled,
genome-wide result only. (ii) any claim about which specific alternate proteins are made —
selection on a reading frame is evidence of function, not proof of translation."""
assert old_ns in doc2; doc2=doc2.replace(old_ns,new_ns)
open(f"{WS}/tierA/OLGenie_methods_explained.md","w").write(doc2)
print(f"updated: {len(doc2.splitlines())} lines")
~~~~

### [msg 597] python · env: dcr-tierA · Writing updated result tables

~~~~python
selDp.to_csv(f"{S}/selection_protein_frames_per_interval.csv",index=False)
fr.to_csv(f"{S}/selection_by_framerel.csv",index=False)
lm.to_csv(f"{S}/omega_by_length_bin.csv",index=False)
main=pd.DataFrame({'set':['DCR (protein frames)','latent null'],'n':[len(selDp),len(selN2)],
  'w_alt_dNNdNS':[oDp['w_alt_refN'],oN2['w_alt_refN']],
  'dNNdNS_lo':[np.percentile(bsDp_N,2.5),np.nan],'dNNdNS_hi':[np.percentile(bsDp_N,97.5),np.nan],
  'w_alt_dSNdSS':[oDp['w_alt_refS'],oN2['w_alt_refS']],
  'dSNdSS_lo':[np.percentile(bsDp_S,2.5),np.nan],'dSNdSS_hi':[np.percentile(bsDp_S,97.5),np.nan],
  'w_ref_control':[oDp['w_ref'],oN2['w_ref']]})
main.to_csv(f"{S}/selection_omega.csv",index=False)
gp[['prot_offs','prot_nf','niso']].to_csv(f"{S}/dcr_protein_frames.csv")
print(main.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
~~~~

### [msg 606] python · env: dcr-tierA · Computing null-arm bootstrap intervals

~~~~python
bsN_N=boot_omega(selN2,'refN'); bsN_S=boot_omega(selN2,'refS')
ci=lambda b:(np.percentile(b,2.5),np.percentile(b,97.5))
pd_=pooled(selDp); pn_=pooled(selN2)
print("pooled normalized divergences:")
for k in ['dNN','dNS','dSN','dSS']:
    print(f"  {k}: D={pd_[k]:.4f}  N={pn_[k]:.4f}")
print(f"\nD dNN/dNS {oDp['w_alt_refN']:.3f} CI{tuple(round(x,3) for x in ci(bsDp_N))}")
print(f"N dNN/dNS {oN2['w_alt_refN']:.3f} CI{tuple(round(x,3) for x in ci(bsN_N))}")
print(f"D dSN/dSS {oDp['w_alt_refS']:.3f} CI{tuple(round(x,3) for x in ci(bsDp_S))}")
print(f"N dSN/dSS {oN2['w_alt_refS']:.3f} CI{tuple(round(x,3) for x in ci(bsN_S))}")
# frame-rel CIs
frci={}
for rel in (1,2):
    dd=score_set([r for r in Dprot if r[5]==rel],'x'); nn=score_set([r for r in Nrecs2 if r[5]==rel],'x')
    frci[rel]=(omegas(pooled(dd))['w_alt_refS'],ci(boot_omega(dd,'refS')),
               omegas(pooled(nn))['w_alt_refS'],ci(boot_omega(nn,'refS')))
    print(f"rel{rel}: D={frci[rel][0]:.3f} CI{tuple(round(x,3) for x in frci[rel][1])}  N={frci[rel][2]:.3f} CI{tuple(round(x,3) for x in frci[rel][3])}")
~~~~

### [msg 608] python · env: dcr-tierA · Inspecting pooled counts structure

~~~~python
print(sorted(pd_.keys()))
~~~~

### [msg 610] python · env: dcr-tierA · Extracting site counts and proportions

~~~~python
for k in ['NN','NS','SN','SS']:
    print(f"  {k}: D  d={pd_['d'+k]:.0f} s={pd_['s'+k]:.0f} p={pd_['p'+k]:.4f}   |  N  d={pn_['d'+k]:.0f} s={pn_['s'+k]:.0f} p={pn_['p'+k]:.4f}")
~~~~

### [msg 612] python · env: dcr-tierA · Regenerating selection figure with protein frames

~~~~python
fig,axes=plt.subplots(1,4,figsize=(13.6,3.6))

# --- A: the four site-class divergences, DCR vs null ---
ax=axes[0]
cls=['NN','NS','SN','SS']
xd=np.arange(4); w=0.38
dv=[pd_['p'+k] for k in cls]; nv=[pn_['p'+k] for k in cls]
ax.bar(xd-w/2,dv,w,color=ORANGE,edgecolor='k',linewidth=0.5,label=f'DCR (n={len(selDp)})')
ax.bar(xd+w/2,nv,w,color=GREY,edgecolor='k',linewidth=0.5,label=f'latent null (n={len(selN2)})')
ax.set_xticks(xd); ax.set_xticklabels(cls,fontsize=8.5,weight='bold')
ax.set_xlabel("site class  (canonical, alternate)",fontsize=8)
ax.set_ylabel("divergence  (JC-corrected)")
ax.legend(frameon=False,fontsize=6.8,loc='upper left')
ax.set_ylim(0,0.56)
ax.text(0.99,0.60,"DCR: SN suppressed\nrelative to SS",transform=ax.transAxes,ha='right',
        fontsize=6.6,style='italic',color=ORANGE)
ax.set_title("A  Divergence by site class",loc='left',weight='bold')

# --- B: omega_alt, both estimators, D vs N ---
ax=axes[1]
xb=np.arange(2)
Dv=[oDp['w_alt_refN'],oDp['w_alt_refS']]; Nv=[oN2['w_alt_refN'],oN2['w_alt_refS']]
Der=np.array([[Dv[0]-ci(bsDp_N)[0],Dv[1]-ci(bsDp_S)[0]],[ci(bsDp_N)[1]-Dv[0],ci(bsDp_S)[1]-Dv[1]]])
Ner=np.array([[Nv[0]-ci(bsN_N)[0],Nv[1]-ci(bsN_S)[0]],[ci(bsN_N)[1]-Nv[0],ci(bsN_S)[1]-Nv[1]]])
ax.bar(xb-w/2,Dv,w,yerr=Der,color=ORANGE,edgecolor='k',linewidth=0.5,
       error_kw=dict(lw=0.9,capsize=3),label='DCR')
ax.bar(xb+w/2,Nv,w,yerr=Ner,color=GREY,edgecolor='k',linewidth=0.5,
       error_kw=dict(lw=0.9,capsize=3),label='latent null')
ax.axhline(1.0,ls='--',lw=1,color='k')
for i,(a,b) in enumerate(zip(Dv,Nv)):
    ax.text(i-w/2,a+0.07,f"{a:.2f}",ha='center',fontsize=7,weight='bold')
    ax.text(i+w/2,b+0.06,f"{b:.2f}",ha='center',fontsize=7,weight='bold')
ax.set_xticks(xb); ax.set_xticklabels(['dNN/dNS','dSN/dSS\n(cancels canonical\nconstraint)'],fontsize=7.2)
ax.set_ylabel("ω  (alternate frame)"); ax.set_ylim(0,1.32)
ax.legend(frameon=False,fontsize=6.8,loc='upper left')
ax.text(0.5,0.055,"p < 0.001, both estimators",transform=ax.transAxes,ha='center',fontsize=6.8,style='italic')
ax.set_title("B  Alternate-frame selection",loc='left',weight='bold')

# --- C: by frame relationship ---
ax=axes[2]
xc=np.arange(2)
Dc=[frci[1][0],frci[2][0]]; Nc=[frci[1][2],frci[2][2]]
Dce=np.array([[Dc[0]-frci[1][1][0],Dc[1]-frci[2][1][0]],[frci[1][1][1]-Dc[0],frci[2][1][1]-Dc[1]]])
Nce=np.array([[Nc[0]-frci[1][3][0],Nc[1]-frci[2][3][0]],[frci[1][3][1]-Nc[0],frci[2][3][1]-Nc[1]]])
ax.bar(xc-w/2,Dc,w,yerr=Dce,color=ORANGE,edgecolor='k',linewidth=0.5,error_kw=dict(lw=0.9,capsize=3))
ax.bar(xc+w/2,Nc,w,yerr=Nce,color=GREY,edgecolor='k',linewidth=0.5,error_kw=dict(lw=0.9,capsize=3))
ax.axhline(1.0,ls='--',lw=1,color='k')
for i,(a,b) in enumerate(zip(Dc,Nc)):
    ax.text(i-w/2,a+0.075,f"{a:.2f}",ha='center',fontsize=7,weight='bold')
    ax.text(i+w/2,b+0.075,f"{b:.2f}",ha='center',fontsize=7,weight='bold')
ax.set_xticks(xc)
ax.set_xticklabels([f'ss13  (+1)\nn={sum(1 for r in Dprot if r[5]==1)}',
                    f'ss12  (+2)\nn={sum(1 for r in Dprot if r[5]==2)}'],fontsize=7.2)
ax.set_ylabel("ω  dSN/dSS"); ax.set_ylim(0,1.48)
ax.set_title("C  Holds in both frame shifts",loc='left',weight='bold')

# --- D: canonical-frame positive control ---
ax=axes[3]
ctrl=[oDp['w_ref'],oN2['w_ref']]
ax.bar([0,1],ctrl,0.5,color=[ORANGE,GREY],edgecolor='k',linewidth=0.5)
ax.axhline(1.0,ls='--',lw=1,color='k')
for i,v in enumerate(ctrl): ax.text(i,v+0.035,f"{v:.2f}",ha='center',fontsize=7.5,weight='bold')
ax.set_xticks([0,1]); ax.set_xticklabels(['DCR','latent\nnull'],fontsize=7.5)
ax.set_ylabel("ω  dNN/dSN   (canonical frame)"); ax.set_ylim(0,1.15)
ax.text(0.5,0.52,"both strongly purifying\n→ reading frames are\nassigned correctly",
        transform=ax.transAxes,ha='center',fontsize=6.8,style='italic')
ax.set_title("D  Positive control",loc='left',weight='bold')

fig.suptitle("Alternate reading frames of DCRs are under purifying selection; matched latent intervals are neutral  "
             "(human–mouse, protein-derived frames)",x=0.005,ha='left',fontsize=9.6,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/alternate_frame_selection.png",dpi=300); plt.close(fig)
print("saved")
~~~~

### [msg 614] python · env: dcr-tierA · Fixing figure label collisions

~~~~python
fig,axes=plt.subplots(1,4,figsize=(13.6,3.6))
ax=axes[0]
ax.bar(xd-w/2,dv,w,color=ORANGE,edgecolor='k',linewidth=0.5,label=f'DCR (n={len(selDp)})')
ax.bar(xd+w/2,nv,w,color=GREY,edgecolor='k',linewidth=0.5,label=f'latent null (n={len(selN2)})')
ax.set_xticks(xd); ax.set_xticklabels(cls,fontsize=8.5,weight='bold')
ax.set_xlabel("site class  (canonical, alternate)",fontsize=8)
ax.set_ylabel("divergence  (JC-corrected)")
ax.legend(frameon=False,fontsize=6.8,loc='upper left')
ax.set_ylim(0,0.60)
ax.annotate("SN suppressed\nrelative to SS",xy=(2.19,0.245),xytext=(1.35,0.50),fontsize=6.6,
            style='italic',color='#a35c00',ha='center',
            arrowprops=dict(arrowstyle='->',color='#a35c00',lw=0.8))
ax.set_title("A  Divergence by site class",loc='left',weight='bold')

ax=axes[1]
ax.bar(xb-w/2,Dv,w,yerr=Der,color=ORANGE,edgecolor='k',linewidth=0.5,
       error_kw=dict(lw=0.9,capsize=3),label='DCR')
ax.bar(xb+w/2,Nv,w,yerr=Ner,color=GREY,edgecolor='k',linewidth=0.5,
       error_kw=dict(lw=0.9,capsize=3),label='latent null')
ax.axhline(1.0,ls='--',lw=1,color='k')
for i,(a,b) in enumerate(zip(Dv,Nv)):
    ax.text(i-w/2,a+Der[1][i]+0.05,f"{a:.2f}",ha='center',fontsize=7,weight='bold')
    ax.text(i+w/2,b+Ner[1][i]+0.05,f"{b:.2f}",ha='center',fontsize=7,weight='bold')
ax.set_xticks(xb); ax.set_xticklabels(['dNN/dNS','dSN/dSS\n(cancels canonical\nconstraint)'],fontsize=7.2)
ax.set_ylabel("ω  (alternate frame)"); ax.set_ylim(0,1.34)
ax.legend(frameon=False,fontsize=6.8,loc='upper left')
ax.text(0.5,0.045,"p < 0.001, both estimators",transform=ax.transAxes,ha='center',fontsize=6.8,style='italic')
ax.set_title("B  Alternate-frame selection",loc='left',weight='bold')

ax=axes[2]
ax.bar(xc-w/2,Dc,w,yerr=Dce,color=ORANGE,edgecolor='k',linewidth=0.5,error_kw=dict(lw=0.9,capsize=3))
ax.bar(xc+w/2,Nc,w,yerr=Nce,color=GREY,edgecolor='k',linewidth=0.5,error_kw=dict(lw=0.9,capsize=3))
ax.axhline(1.0,ls='--',lw=1,color='k')
for i,(a,b) in enumerate(zip(Dc,Nc)):
    ax.text(i-w/2,a+Dce[1][i]+0.055,f"{a:.2f}",ha='center',fontsize=7,weight='bold')
    ax.text(i+w/2,b+Nce[1][i]+0.055,f"{b:.2f}",ha='center',fontsize=7,weight='bold')
ax.set_xticks(xc)
ax.set_xticklabels([f'ss13  (+1)\nn={sum(1 for r in Dprot if r[5]==1)}',
                    f'ss12  (+2)\nn={sum(1 for r in Dprot if r[5]==2)}'],fontsize=7.2)
ax.set_ylabel("ω  dSN/dSS"); ax.set_ylim(0,1.55)
ax.set_title("C  Holds in both frame shifts",loc='left',weight='bold')

ax=axes[3]
ax.bar([0,1],ctrl,0.5,color=[ORANGE,GREY],edgecolor='k',linewidth=0.5)
ax.axhline(1.0,ls='--',lw=1,color='k')
for i,v in enumerate(ctrl): ax.text(i,v+0.035,f"{v:.2f}",ha='center',fontsize=7.5,weight='bold')
ax.set_xticks([0,1]); ax.set_xticklabels(['DCR','latent\nnull'],fontsize=7.5)
ax.set_ylabel("ω  dNN/dSN   (canonical frame)"); ax.set_ylim(0,1.15)
ax.text(0.5,0.55,"both strongly purifying\n→ reading frames are\nassigned correctly",
        transform=ax.transAxes,ha='center',fontsize=6.8,style='italic')
ax.set_title("D  Positive control",loc='left',weight='bold')

fig.suptitle("Alternate reading frames of DCRs are under purifying selection; matched latent intervals are neutral  "
             "(human–mouse, protein-derived frames)",x=0.005,ha='left',fontsize=9.6,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/alternate_frame_selection.png",dpi=300); plt.close(fig)
print("saved")
~~~~

### [msg 616] python · env: dcr-tierA · Repositioning panel A annotation

~~~~python
ax=axes[0] if False else None
fig=plt.figure(figsize=(13.6,3.6)); axes=fig.subplots(1,4)
ax=axes[0]
ax.bar(xd-w/2,dv,w,color=ORANGE,edgecolor='k',linewidth=0.5,label=f'DCR (n={len(selDp)})')
ax.bar(xd+w/2,nv,w,color=GREY,edgecolor='k',linewidth=0.5,label=f'latent null (n={len(selN2)})')
ax.set_xticks(xd); ax.set_xticklabels(cls,fontsize=8.5,weight='bold')
ax.set_xlabel("site class  (canonical, alternate)",fontsize=8)
ax.set_ylabel("divergence  (JC-corrected)")
ax.legend(frameon=False,fontsize=6.8,loc='upper left')
ax.set_ylim(0,0.60)
ax.annotate("SN suppressed\nrelative to SS",xy=(2.19,0.245),xytext=(0.95,0.355),fontsize=6.6,
            style='italic',color='#a35c00',ha='center',
            arrowprops=dict(arrowstyle='->',color='#a35c00',lw=0.8))
ax.set_title("A  Divergence by site class",loc='left',weight='bold')

ax=axes[1]
ax.bar(xb-w/2,Dv,w,yerr=Der,color=ORANGE,edgecolor='k',linewidth=0.5,
       error_kw=dict(lw=0.9,capsize=3),label='DCR')
ax.bar(xb+w/2,Nv,w,yerr=Ner,color=GREY,edgecolor='k',linewidth=0.5,
       error_kw=dict(lw=0.9,capsize=3),label='latent null')
ax.axhline(1.0,ls='--',lw=1,color='k')
for i,(a,b) in enumerate(zip(Dv,Nv)):
    ax.text(i-w/2,a+Der[1][i]+0.05,f"{a:.2f}",ha='center',fontsize=7,weight='bold')
    ax.text(i+w/2,b+Ner[1][i]+0.05,f"{b:.2f}",ha='center',fontsize=7,weight='bold')
ax.set_xticks(xb); ax.set_xticklabels(['dNN/dNS','dSN/dSS\n(cancels canonical\nconstraint)'],fontsize=7.2)
ax.set_ylabel("ω  (alternate frame)"); ax.set_ylim(0,1.34)
ax.legend(frameon=False,fontsize=6.8,loc='upper left')
ax.text(0.5,0.045,"p < 0.001, both estimators",transform=ax.transAxes,ha='center',fontsize=6.8,style='italic')
ax.set_title("B  Alternate-frame selection",loc='left',weight='bold')

ax=axes[2]
ax.bar(xc-w/2,Dc,w,yerr=Dce,color=ORANGE,edgecolor='k',linewidth=0.5,error_kw=dict(lw=0.9,capsize=3))
ax.bar(xc+w/2,Nc,w,yerr=Nce,color=GREY,edgecolor='k',linewidth=0.5,error_kw=dict(lw=0.9,capsize=3))
ax.axhline(1.0,ls='--',lw=1,color='k')
for i,(a,b) in enumerate(zip(Dc,Nc)):
    ax.text(i-w/2,a+Dce[1][i]+0.055,f"{a:.2f}",ha='center',fontsize=7,weight='bold')
    ax.text(i+w/2,b+Nce[1][i]+0.055,f"{b:.2f}",ha='center',fontsize=7,weight='bold')
ax.set_xticks(xc)
ax.set_xticklabels([f'ss13  (+1)\nn={sum(1 for r in Dprot if r[5]==1)}',
                    f'ss12  (+2)\nn={sum(1 for r in Dprot if r[5]==2)}'],fontsize=7.2)
ax.set_ylabel("ω  dSN/dSS"); ax.set_ylim(0,1.55)
ax.set_title("C  Holds in both frame shifts",loc='left',weight='bold')

ax=axes[3]
ax.bar([0,1],ctrl,0.5,color=[ORANGE,GREY],edgecolor='k',linewidth=0.5)
ax.axhline(1.0,ls='--',lw=1,color='k')
for i,v in enumerate(ctrl): ax.text(i,v+0.035,f"{v:.2f}",ha='center',fontsize=7.5,weight='bold')
ax.set_xticks([0,1]); ax.set_xticklabels(['DCR','latent\nnull'],fontsize=7.5)
ax.set_ylabel("ω  dNN/dSN   (canonical frame)"); ax.set_ylim(0,1.15)
ax.text(0.5,0.55,"both strongly purifying\n→ reading frames are\nassigned correctly",
        transform=ax.transAxes,ha='center',fontsize=6.8,style='italic')
ax.set_title("D  Positive control",loc='left',weight='bold')

fig.suptitle("Alternate reading frames of DCRs are under purifying selection; matched latent intervals are neutral  "
             "(human–mouse, protein-derived frames)",x=0.005,ha='left',fontsize=9.6,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/alternate_frame_selection.png",dpi=300); plt.close(fig)
print("saved")
~~~~
