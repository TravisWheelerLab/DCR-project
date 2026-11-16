# 10_xprime_conservation_filter — as executed, not refactored

**Conservation as a specificity filter, X vs X' (pyliftover, then UCSC liftOver re-lift); builds dcr_regions.fa, N/D human and mouse FASTAs.**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[393..449]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 22 (msg 393 … msg 449) |
| Language / environment | bash (env: dcr-tierA); python (env: dcr-tierA) |
| Saved outputs | `xprime_conservation.png`, `xprime_length_strata.csv`, `xprime_summary.csv`, `xprime_robustness.csv` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 22 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

Environment operations in this span (not code cells):

* `[msg 399] manage_packages: {"environment": "dcr-tierA", "mode": "install", "packages": ["pyliftover"], "use_pip": true}`

---

### [msg 393] bash · env: dcr-tierA · Checking mm39 and liftOver availability

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
echo "=== mm39.fa status in both granted dirs ==="
for P in "$REF/mm39.fa" "$D/data/mm39.fa"; do
  printf "%s -> " "$P"
  if [ -e "$P" ]; then echo "EXISTS ($(du -h "$P" 2>/dev/null | cut -f1))"; 
  elif [ -L "$P" ]; then echo "DEAD SYMLINK -> $(readlink "$P")"; else echo "MISSING"; fi
done
echo ""
echo "=== any mm39 fasta anywhere visible ==="
ls -la "$REF" | grep -i mm39
echo "--- DOOVER/data ---"
ls -la "$D/data" 2>/dev/null | grep -i -E "mm39|chain" | head
echo ""
echo "=== liftOver / chain ==="
which liftOver CrossMap.py 2>/dev/null || echo "liftOver: MISSING"
ls -la "$REF/hg38ToMm39.over.chain" 2>/dev/null | awk '{print "chain:",$5,"bytes"}'
~~~~

### [msg 397] bash · env: dcr-tierA · Verifying mm39 genome readability

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
M="<WHEELER_LAB_DRIVE>/2024/07-02-DCE_Conservation/mm39.fa"
ls -la "$M" 2>&1 | head -2
if [ -r "$M" ]; then
  echo "READABLE: $(du -h "$M" | cut -f1)"
  head -1 "$M"
  ls -la "${M}.fai" 2>/dev/null | awk '{print "fai:",$5,"bytes"}'
else echo "STILL UNREADABLE"; fi
~~~~

### [msg 403] bash · env: dcr-tierA · Inspecting exon bed strand encoding

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
echo "=== original exon beds: format ==="
for F in DCR_exons_DCR_genes_UNIQ nonDCR_exon_DCR_gene_UNIQ nonDCR_exon_nonDCR_gene_UNIQ; do
  P=$(find "$D/output" -name "${F}.bed" | head -1)
  echo "--- $F : $P"
  head -2 "$P"
  awk '{print NF}' "$P" | sort -u | head -3 | tr '\n' ' '; echo "cols"
done
echo ""
echo "=== strand encoding in name field ==="
P=$(find "$D/output" -name "DCR_exons_DCR_genes_UNIQ.bed" | head -1)
awk '{n=split($4,a,"_"); print a[n]}' "$P" | sort | uniq -c | head
~~~~

### [msg 405] bash · env: dcr-tierA · Extracting stranded non-DCR exon pool

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
X="$WS/tierA/xprime"; mkdir -p "$X"
C="$D/output/02_Conservation"
# stranded beds from the _f/_r suffix; sample non-DCR pool for tractability
for F in nonDCR_exon_DCR_gene_UNIQ nonDCR_exon_nonDCR_gene_UNIQ; do
  awk 'BEGIN{OFS="\t"}{n=split($4,a,"_"); s=(a[n]=="r")?"-":"+"; print $1,$2,$3,$4,0,s}' "$C/$F.bed"
done | shuf -n 40000 --random-source=<(yes) > "$X/nondcr_pool.bed"
wc -l < "$X/nondcr_pool.bed"; head -2 "$X/nondcr_pool.bed"
bedtools getfasta -fi "$REF/hg38.fa" -bed "$X/nondcr_pool.bed" -s -name -fo "$X/nondcr_pool.fa"
grep -c '^>' "$X/nondcr_pool.fa"; head -1 "$X/nondcr_pool.fa"
~~~~

### [msg 407] bash · env: dcr-tierA · Extracting stranded exon and DCR sequences

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
X="$WS/tierA/xprime"; mkdir -p "$X"
C="$D/output/02_Conservation"
for F in nonDCR_exon_DCR_gene_UNIQ nonDCR_exon_nonDCR_gene_UNIQ; do
  awk -v F="$F" 'BEGIN{OFS="\t"}{n=split($4,a,"_"); s=(a[n]=="r")?"-":"+"; print $1,$2,$3,F"|"$4,0,s}' "$C/$F.bed"
done > "$X/nondcr_all.bed"
wc -l < "$X/nondcr_all.bed"
# deterministic subsample: every Nth line
awk 'NR%8==1' "$X/nondcr_all.bed" > "$X/nondcr_pool.bed"
wc -l < "$X/nondcr_pool.bed"; head -2 "$X/nondcr_pool.bed"
bedtools getfasta -fi "$REF/hg38.fa" -bed "$X/nondcr_pool.bed" -s -name -fo "$X/nondcr_pool.fa"
grep -c '^>' "$X/nondcr_pool.fa"; head -1 "$X/nondcr_pool.fa"
# stranded human DCR regions (D)
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"dce"$7"|"NR,0,$5}' "$D/output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed" > "$X/dcr_regions.bed"
bedtools getfasta -fi "$REF/hg38.fa" -bed "$X/dcr_regions.bed" -s -name -fo "$X/dcr_regions.fa"
grep -c '^>' "$X/dcr_regions.fa"
~~~~

### [msg 409] python · env: dcr-tierA · Identifying latent dual-coding intervals

~~~~python
import re
X=f"{WS}/tierA/xprime"
def read_fa2(path):
    out={};hdr=None;cur=[]
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if hdr is not None: out[hdr]=''.join(cur).upper()
                hdr=line[1:].strip(); cur=[]
            else: cur.append(line.strip())
        if hdr is not None: out[hdr]=''.join(cur).upper()
    return out

def best_dual_interval(seq, min_len=21):
    """Largest sub-interval covered by >=2 ORFs in different frames.
       Returns (start0, end_excl, span) in sequence coords, or None."""
    orfs=esl_translate_watson(seq,min_aa=1)
    best=None; bs=0
    for i in range(len(orfs)):
        s1,e1,f1,_=orfs[i]
        for j in range(i+1,len(orfs)):
            s2,e2,f2,_=orfs[j]
            if f1==f2: continue
            a=max(s1,s2); b=min(e1,e2)
            span=b-a+1
            if span>=min_len and span>bs:
                bs=span; best=(a-1,b)
    return (best[0],best[1],bs) if best else None

nd=read_fa2(f"{X}/nondcr_pool.fa")
print("non-DCR pool sequences:",len(nd))
# parse header -> genomic coords
pat=re.compile(r'^(.+?)\|(.+?)::(chr[^:]+):(\d+)-(\d+)\(([+-])\)$')
rows=[]
for h,s in nd.items():
    mm=pat.match(h)
    if not mm: continue
    cls,name,c,st,en,strand=mm.groups(); st=int(st); en=int(en)
    bd=best_dual_interval(s,21)
    if bd is None: continue
    a,b,span=bd
    # map sequence coords back to genome (seq is already strand-oriented)
    if strand=='+': gs,ge=st+a, st+b
    else:           gs,ge=en-b, en-a
    rows.append((c,gs,ge,cls,name,strand,span,len(s)))
N=pd.DataFrame(rows,columns=['chr','start','end','cls','name','strand','span','exlen'])
print("non-DCR exons with a latent dual-coding interval >=21nt:",len(N),f"({len(N)/len(nd)*100:.1f}%)")
print("median latent interval span:",N.span.median(),"nt   | DCR region median: 64 nt")
N.to_csv(f"{X}/N_latent_intervals.tsv",sep="\t",index=False)
~~~~

### [msg 411] python · env: dcr-tierA · Lifting N and D intervals to mouse

~~~~python
from pyliftover import LiftOver
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
lo=LiftOver(f"{REF}/hg38ToMm39.over.chain")
print("chain loaded")
def lift_iv(c,s,e,strand):
    a=lo.convert_coordinate(c,s,strand); b=lo.convert_coordinate(c,e-1,strand)
    if not a or not b: return None
    c1,p1,st1,_=a[0]; c2,p2,st2,_=b[0]
    if c1!=c2 or st1!=st2: return None
    lo_,hi=min(p1,p2),max(p1,p2)+1
    if hi-lo_< 15 or hi-lo_ > 5*(e-s)+50: return None
    return (c1,lo_,hi,st1)
# lift D (real DCR regions)
dcr_bed=pd.read_csv(f"{X}/dcr_regions.bed",sep="\t",header=None,names=['chr','start','end','name','score','strand'])
dl=[]
for t in dcr_bed.itertuples():
    r=lift_iv(t.chr,t.start,t.end,t.strand)
    dl.append((t.name,)+r if r else (t.name,None,None,None,None))
D_lift=pd.DataFrame(dl,columns=['name','mchr','mstart','mend','mstrand'])
print("D: %d/%d lifted (%.1f%%)"%(D_lift.mchr.notna().sum(),len(D_lift),D_lift.mchr.notna().mean()*100))
# lift N (latent intervals) 
nl=[]
for t in N.itertuples():
    r=lift_iv(t.chr,t.start,t.end,t.strand)
    nl.append((t.Index,)+r if r else (t.Index,None,None,None,None))
N_lift=pd.DataFrame(nl,columns=['idx','mchr','mstart','mend','mstrand'])
print("N: %d/%d lifted (%.1f%%)"%(N_lift.mchr.notna().sum(),len(N_lift),N_lift.mchr.notna().mean()*100))
~~~~

### [msg 413] python · env: dcr-tierA · Extracting mouse ortholog sequences

~~~~python
import subprocess, os
MM="<WHEELER_LAB_DRIVE>/2024/07-02-DCE_Conservation/mm39.fa"
def write_bed_and_extract(df, cols, out_prefix):
    ok=df[df.mchr.notna()].copy()
    bed=f"{X}/{out_prefix}.bed"; fa=f"{X}/{out_prefix}.fa"
    with open(bed,"w") as fh:
        for t in ok.itertuples():
            fh.write(f"{t.mchr}\t{int(t.mstart)}\t{int(t.mend)}\t{getattr(t,cols)}\t0\t{t.mstrand}\n")
    env=dict(os.environ); env["PATH"]="/opt/homebrew/bin:"+env.get("PATH","")
    r=subprocess.run(["bedtools","getfasta","-fi",MM,"-bed",bed,"-s","-name","-fo",fa],
                     capture_output=True,text=True,env=env)
    if r.returncode!=0: print("ERR",r.stderr[:300])
    return fa, len(ok)
faD,nD=write_bed_and_extract(D_lift,'name','D_mouse')
faN,nN=write_bed_and_extract(N_lift,'idx','N_mouse')
print("extracted D:",nD," N:",nN)
seqD=read_fa2(faD); seqN=read_fa2(faN)
print("seqs D:",len(seqD)," N:",len(seqN))
~~~~

### [msg 415] bash · env: dcr-tierA · Extracting mouse sequences in fresh process

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
X="$WS/tierA/xprime"
MM="<WHEELER_LAB_DRIVE>/2024/07-02-DCE_Conservation/mm39.fa"
ls -la "$MM" | awk '{print "mm39:",$5}'
for P in D_mouse N_mouse; do
  bedtools getfasta -fi "$MM" -bed "$X/$P.bed" -s -name -fo "$X/$P.fa" 2>&1 | head -2
  echo "$P: $(grep -c '^>' "$X/$P.fa" 2>/dev/null) seqs"
done
~~~~

### [msg 417] python · env: dcr-tierA · Scoring mouse dual-coding for N and D

~~~~python
seqD=read_fa2(f"{X}/D_mouse.fa"); seqN=read_fa2(f"{X}/N_mouse.fa")
print("loaded D:",len(seqD)," N:",len(seqN))
def is_dual_strict(s):
    L=len(s); orfs=esl_translate_watson(s,min_aa=1)
    return sum(1 for o in orfs if (o[0]-1-3)<=0 and (L-o[1])<=3)>1
def is_dual_sub(s,minlen=21):
    return best_dual_interval(s,minlen) is not None
Dm=pd.DataFrame([(h,len(s),is_dual_strict(s),is_dual_sub(s)) for h,s in seqD.items()],
                columns=['h','len','strict','sub'])
Nm=pd.DataFrame([(h,len(s),is_dual_strict(s),is_dual_sub(s)) for h,s in seqN.items()],
                columns=['h','len','strict','sub'])
print(f"\nD' (real DCRs, mouse):     strict {Dm.strict.mean()*100:.1f}%   sub {Dm['sub'].mean()*100:.1f}%   n={len(Dm)}")
print(f"N' (latent, mouse):        strict {Nm.strict.mean()*100:.1f}%   sub {Nm['sub'].mean()*100:.1f}%   n={len(Nm)}")
print(f"\n(published D mouse-ortholog strict = 82.4%; this D' = {Dm.strict.mean()*100:.1f}% - consistency check)")
~~~~

### [msg 419] python · env: dcr-tierA · Length-stratifying mouse retention

~~~~python
Dm['name']=[h.split('::')[0] for h in Dm.h]; Nm['idx']=[int(h.split('::')[0]) for h in Nm.h]
# attach human-side length & GC for matching
def gcf(s): return (s.count('G')+s.count('C'))/max(1,len(s))
hD=read_fa2(f"{X}/dcr_regions.fa"); hN=read_fa2(f"{X}/nondcr_pool.fa")
hD_map={h.split('::')[0]:v for h,v in hD.items()}
Dm['hlen']=[len(hD_map.get(n,'')) for n in Dm.name]; Dm['hgc']=[gcf(hD_map.get(n,'N')) for n in Dm.name]
Nsub=N.reset_index().rename(columns={'index':'idx'})
Nm=Nm.merge(Nsub[['idx','span','chr','start','end','strand']],on='idx',how='left')
Nm['hlen']=Nm.span
# GC of the latent interval itself
seqNh={}
for h,s in hN.items():
    mm=pat.match(h)
    if mm: seqNh[mm.group(2)]=s
Nm['hgc']=np.nan
tmp=Nsub.set_index('idx')
Nm['hgc']=[gcf(seqNh[tmp.loc[i,'name']][int(tmp.loc[i,'start'])-0:0] ) if False else np.nan for i in Nm.idx[:0]] if False else Nm['hgc']
# simpler: recompute GC from mouse-lifted human interval via genome-free approach - use span+exon GC proxy
print("D': strict %.1f%% (n=%d), median human len %d"%(Dm.strict.mean()*100,len(Dm),Dm.hlen.median()))
print("N': strict %.1f%% (n=%d), median human len %d"%(Nm.strict.mean()*100,len(Nm),Nm.hlen.median()))
# length-matched comparison
bins=[0,30,45,60,90,120,180,300,10**6]
Dm['lb']=pd.cut(Dm.hlen,bins); Nm['lb']=pd.cut(Nm.hlen,bins)
cmp=pd.DataFrame({'D_prime':Dm.groupby('lb',observed=True).strict.mean()*100,
                  'nD':Dm.groupby('lb',observed=True).strict.size(),
                  'N_prime':Nm.groupby('lb',observed=True).strict.mean()*100,
                  'nN':Nm.groupby('lb',observed=True).strict.size()})
cmp['gap']=cmp.D_prime-cmp.N_prime
print("\nLength-stratified mouse dual-coding retention:")
print(cmp.to_string(float_format=lambda v:f"{v:.1f}"))
~~~~

### [msg 421] python · env: dcr-tierA · Computing length-matched X and X-prime

~~~~python
# Length-matched N' rate: reweight N by D's length distribution
w=cmp.nD/cmp.nD.sum()
N_matched=float((cmp.N_prime*w).sum()); D_raw=float(Dm.strict.mean()*100)
print(f"D' retention                     : {D_raw:.1f}%")
print(f"N' retention (raw)               : {Nm.strict.mean()*100:.1f}%")
print(f"N' retention (length-matched to D): {N_matched:.1f}%")
print(f"--> length-matched gap = {D_raw-N_matched:.1f} pts; retention ratio = {D_raw/N_matched:.3f}x")

# Now X and X' at genome scale. N pool was a 1/8 subsample of 307,523 non-DCR exons.
SCALE=8
N_total = len(N)*SCALE                     # human exons with latent dual-coding potential
D_total = len(dcr_bed)                     # real DCR intervals
Xv = N_total/(N_total+D_total)
Np = N_lift.mchr.notna().mean()*(N_matched/100)*N_total   # lifted & retained (length-matched rate)
Dp = D_lift.mchr.notna().mean()*(D_raw/100)*D_total
Xp = Np/(Np+Dp)
print(f"\nN  = {N_total:,}   D  = {D_total:,}   ->  X  = {Xv:.5f}")
print(f"N' = {Np:,.0f}   D' = {Dp:,.0f}   ->  X' = {Xp:.5f}")
print(f"\nspecificity gain (1-X)/(1-X') = {(1-Xp)/(1-Xv):.3f}x   [1.0 = conservation adds nothing]")
# bootstrap the retention ratio
rng2=np.random.default_rng(11); ratios=[]
dv=Dm.strict.values.astype(float); nv=Nm.strict.values.astype(float); nlb=Nm.lb.values
for _ in range(1000):
    db=rng2.choice(dv,len(dv),replace=True)
    idx=rng2.integers(0,len(nv),len(nv)); nb=pd.DataFrame({'s':nv[idx],'lb':nlb[idx]})
    nm_=float((nb.groupby('lb',observed=True).s.mean()*100*w).reindex(w.index).fillna(0).sum())
    ratios.append(db.mean()*100/nm_ if nm_>0 else np.nan)
ratios=np.array(ratios); ratios=ratios[~np.isnan(ratios)]
print(f"retention ratio bootstrap 95% CI: {np.percentile(ratios,2.5):.3f} - {np.percentile(ratios,97.5):.3f}")
from scipy import stats as st
ct=np.array([[Dm.strict.sum(),len(Dm)-Dm.strict.sum()],[Nm.strict.sum(),len(Nm)-Nm.strict.sum()]])
print("Fisher OR=%.3f p=%.3e (unmatched)"%st.fisher_exact(ct))
~~~~

### [msg 424] python · env: dcr-tierA · Testing robustness of the X-prime gap

~~~~python
# Robustness: does the D'/N' gap grow under stricter definitions of "latent potential"?
checks={}
for minspan in [21,45,60,90]:
    sel=Nm[Nm.hlen>=minspan]
    if len(sel)<200: continue
    s2=sel.copy(); s2['lb']=pd.cut(s2.hlen,bins)
    r=s2.groupby('lb',observed=True).strict.mean()*100
    m=float((r*w).reindex(w.index).fillna(0).sum())
    checks[f"latent span>={minspan}"]=(len(sel),m,D_raw-m)
# also: restrict N to exons in DCR-containing genes only (same gene environment)
Nm['cls']=Nm.idx.map(Nsub.set_index('idx').cls)
for c in Nm.cls.dropna().unique():
    sel=Nm[Nm.cls==c].copy(); sel['lb']=pd.cut(sel.hlen,bins)
    r=sel.groupby('lb',observed=True).strict.mean()*100
    m=float((r*w).reindex(w.index).fillna(0).sum())
    checks[c.replace('_UNIQ','')]=(len(sel),m,D_raw-m)
print(f"D' = {D_raw:.1f}%   (reference)\n")
print(f"{'N-set definition':38s} {'n':>7s} {'N-prime':>9s} {'gap':>7s}")
for k,(n,m,g) in checks.items(): print(f"{k:38s} {n:7d} {m:8.1f}% {g:6.1f}")
~~~~

### [msg 426] python · env: dcr-tierA · Correcting length-matching renormalization

~~~~python
def matched_rate(sub_df):
    """length-matched rate with weights RENORMALIZED over bins present in both sets"""
    s=sub_df.copy(); s['lb']=pd.cut(s.hlen,bins)
    r=s.groupby('lb',observed=True).strict.mean()*100
    ww=w.reindex(r.index).dropna()
    r=r.reindex(ww.index)
    ww=ww/ww.sum()
    return float((r*ww).sum()), int(len(s)), list(ww.index.astype(str))
print(f"D' = {D_raw:.1f}%  (reference)\n")
print(f"{'N-set definition':34s} {'n':>7s} {'N-prime':>9s} {'gap':>6s}  bins used")
for minspan in [21,45,60,90]:
    sel=Nm[Nm.hlen>=minspan]
    m,n,bl=matched_rate(sel)
    print(f"{'latent span>='+str(minspan):34s} {n:7d} {m:8.1f}% {D_raw-m:5.1f}  {len(bl)}")
for c in ['nonDCR_exon_DCR_gene_UNIQ','nonDCR_exon_nonDCR_gene_UNIQ']:
    sel=Nm[Nm.cls==c]; m,n,bl=matched_rate(sel)
    print(f"{c.replace('_UNIQ',''):34s} {n:7d} {m:8.1f}% {D_raw-m:5.1f}  {len(bl)}")
~~~~

### [msg 428] python · env: dcr-tierA · Bin-matching both arms fairly

~~~~python
def paired_rates(sub_df):
    s=sub_df.copy(); s['lb']=pd.cut(s.hlen,bins)
    rn=s.groupby('lb',observed=True).strict.mean()*100
    rd=Dm.groupby('lb',observed=True).strict.mean()*100
    common=[b for b in rn.index if b in rd.index and cmp.loc[b,'nD']>=20]
    ww=cmp.loc[common,'nD']; ww=ww/ww.sum()
    return float((rd.reindex(common)*ww).sum()), float((rn.reindex(common)*ww).sum()), len(s)
print(f"{'N-set definition':30s} {'n':>7s} {'D-prime':>8s} {'N-prime':>8s} {'gap':>6s}")
rob=[]
for minspan in [21,45,60,90]:
    d_,n_,cnt=paired_rates(Nm[Nm.hlen>=minspan])
    print(f"{'latent span>='+str(minspan):30s} {cnt:7d} {d_:7.1f}% {n_:7.1f}% {d_-n_:5.1f}")
    rob.append(('span>=%d'%minspan,cnt,d_,n_,d_-n_))
for c,lab in [('nonDCR_exon_DCR_gene_UNIQ','in DCR genes'),('nonDCR_exon_nonDCR_gene_UNIQ','in non-DCR genes')]:
    d_,n_,cnt=paired_rates(Nm[Nm.cls==c])
    print(f"{lab:30s} {cnt:7d} {d_:7.1f}% {n_:7.1f}% {d_-n_:5.1f}")
    rob.append((lab,cnt,d_,n_,d_-n_))
pd.DataFrame(rob,columns=['N_definition','n','D_prime','N_prime','gap']).to_csv(f"{X}/xprime_robustness.csv",index=False)
~~~~

### [msg 431] python · env: dcr-tierA · Running stratified CMH test

~~~~python
from scipy.stats import chi2
# Cochran-Mantel-Haenszel across length strata
tabs=[]
for b in cmp.index:
    d=Dm[Dm.lb==b]; n=Nm[Nm.lb==b]
    if len(d)<20 or len(n)<20: continue
    tabs.append([[d.strict.sum(),len(d)-d.strict.sum()],[n.strict.sum(),len(n)-n.strict.sum()]])
tabs=np.array(tabs,dtype=float)
a=tabs[:,0,0]; b_=tabs[:,0,1]; c_=tabs[:,1,0]; d_=tabs[:,1,1]; nT=tabs.sum(axis=(1,2))
num=(a*d_/nT).sum(); den=(b_*c_/nT).sum(); or_mh=num/den
E=( (a+b_)*(a+c_)/nT ).sum()
V=( (a+b_)*(c_+d_)*(a+c_)*(b_+d_)/(nT**2*(nT-1)) ).sum()
cmh=(abs(a.sum()-E)-0.5)**2/V; p_cmh=1-chi2.cdf(cmh,1)
print(f"CMH (length-stratified): OR = {or_mh:.3f},  chi2 = {cmh:.1f},  p = {p_cmh:.3e},  strata = {len(tabs)}")
print(f"unstratified Fisher OR = 1.903  ->  {(1.903-or_mh)/(1.903-1)*100:.0f}% of the raw association was length confounding")
res_x={'X':Xv,'X_prime':Xp,'D_prime_pct':D_raw,'N_prime_raw_pct':float(Nm.strict.mean()*100),
       'N_prime_matched_pct':N_matched,'ratio':D_raw/N_matched,
       'ratio_lo':float(np.percentile(ratios,2.5)),'ratio_hi':float(np.percentile(ratios,97.5)),
       'OR_MH':or_mh,'p_CMH':p_cmh,'N_total':N_total,'D_total':D_total,
       'lift_rate_D':float(D_lift.mchr.notna().mean()),'lift_rate_N':float(N_lift.mchr.notna().mean())}
pd.DataFrame([res_x]).to_csv(f"{X}/xprime_summary.csv",index=False)
print("\nX  = %.4f    X' = %.4f    specificity gain = %.2fx"%(Xv,Xp,(1-Xp)/(1-Xv)))
~~~~

### [msg 433] python · env: dcr-tierA · Building X-prime figure

~~~~python
fig,axes=plt.subplots(1,3,figsize=(11,3.6))
ax=axes[0]
bl=[str(i) for i in cmp.index]
xp_=np.arange(len(cmp))
ax.plot(xp_,cmp.D_prime,'o-',color=ORANGE,lw=1.8,ms=6,label="D' real DCRs")
ax.plot(xp_,cmp.N_prime,'s-',color=GREY,lw=1.8,ms=5,label="N' latent potential")
ax.set_xticks(xp_); ax.set_xticklabels(['≤30','30-45','45-60','60-90','90-120','120-180','180-300','>300'],
                                        rotation=45,ha='right',fontsize=6.5)
ax.set_xlabel("human interval length (nt)"); ax.set_ylabel("% dual-coding retained in mouse")
ax.legend(frameon=False,fontsize=7); ax.set_ylim(0,108)
ax.set_title("A  Retention is length-driven",loc='left',weight='bold')

ax=axes[1]
bars=[D_raw,float(Nm.strict.mean()*100),N_matched]
ax.bar(range(3),bars,color=[ORANGE,GREY,BLUE],edgecolor='k',linewidth=0.5,width=0.62)
for i,v in enumerate(bars): ax.text(i,v+1.8,f"{v:.1f}%",ha='center',fontsize=8,weight='bold')
ax.set_xticks(range(3)); ax.set_xticklabels(["D'\nreal DCR","N'\nlatent (raw)","N'\nlength-matched"],fontsize=7)
ax.set_ylabel("% dual-coding retained in mouse"); ax.set_ylim(0,100)
ax.annotate("",xy=(2,N_matched+0.5),xytext=(2,D_raw-0.5),arrowprops=dict(arrowstyle='<->',lw=1.1,color='k'))
ax.text(2.06,(D_raw+N_matched)/2,f"{D_raw-N_matched:.1f} pts",fontsize=7,va='center')
ax.text(0.02,0.03,f"CMH OR={or_mh:.2f}, p={p_cmh:.0e}\n61% of raw signal was length",
        transform=ax.transAxes,fontsize=6.5)
ax.set_title("B  Length-matched gap is small",loc='left',weight='bold')

ax=axes[2]
ax.bar([0,1],[Xv*100,Xp*100],color=[GREY,BLUE],edgecolor='k',linewidth=0.5,width=0.5)
for i,v in enumerate([Xv*100,Xp*100]): ax.text(i,v+0.06,f"{v:.2f}%",ha='center',fontsize=8,weight='bold')
ax.set_xticks([0,1]); ax.set_xticklabels(["X\npotential only","X'\n+ mouse conservation"],fontsize=7)
ax.set_ylabel("% of dual-coding-capable regions\nthat are NOT real DCRs")
ax.set_ylim(98.5,100.2)
ax.text(0.5,0.28,f"specificity gain\n{(1-Xp)/(1-Xv):.2f}×",transform=ax.transAxes,ha='center',
        fontsize=9,weight='bold',style='italic')
ax.set_title("C  Conservation barely narrows it",loc='left',weight='bold')
fig.suptitle("Dual-coding potential is necessary but far from sufficient — and conservation adds only modest specificity",
             x=0.005,ha='left',fontsize=9.6,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/xprime_conservation.png",dpi=300); plt.close(fig)
cmp.to_csv(f"{WS}/tierA/xprime_length_strata.csv")
print("saved")
~~~~

### [msg 435] python · env: dcr-tierA · Rebuilding X-prime figure with robustness panel

~~~~python
fig,axes=plt.subplots(1,4,figsize=(13.6,3.6))
ax=axes[0]
xp_=np.arange(len(cmp))
ax.plot(xp_,cmp.D_prime,'o-',color=ORANGE,lw=1.8,ms=6,label="D' real DCRs")
ax.plot(xp_,cmp.N_prime,'s-',color=GREY,lw=1.8,ms=5,label="N' latent potential")
ax.set_xticks(xp_); ax.set_xticklabels(['≤30','30-45','45-60','60-90','90-120','120-180','180-300','>300'],
                                        rotation=45,ha='right',fontsize=6.5)
ax.set_xlabel("human interval length (nt)"); ax.set_ylabel("% dual-coding retained in mouse")
ax.legend(frameon=False,fontsize=7); ax.set_ylim(0,108)
ax.set_title("A  Retention is length-driven",loc='left',weight='bold')

ax=axes[1]
bars=[D_raw,float(Nm.strict.mean()*100),N_matched]
ax.bar(range(3),bars,color=[ORANGE,GREY,BLUE],edgecolor='k',linewidth=0.5,width=0.62)
for i,v in enumerate(bars): ax.text(i,v+2.2,f"{v:.1f}%",ha='center',fontsize=8,weight='bold')
ax.set_xticks(range(3)); ax.set_xticklabels(["D'\nreal DCR","N'\nlatent (raw)","N'\nlength-\nmatched"],fontsize=7)
ax.set_ylabel("% dual-coding retained in mouse"); ax.set_ylim(0,112)
ax.plot([0,0,2,2],[D_raw+6,D_raw+9,D_raw+9,N_matched+6],lw=0.9,color='k')
ax.text(1.0,D_raw+10.5,f"{D_raw-N_matched:.1f} pts",fontsize=7.5,ha='center',weight='bold')
ax.text(0.02,0.03,f"CMH OR={or_mh:.2f}, p={p_cmh:.0e}\n61% of raw signal was length",
        transform=ax.transAxes,fontsize=6.5)
ax.set_title("B  Length-matched gap is small",loc='left',weight='bold')

ax=axes[2]
rb=pd.read_csv(f"{X}/xprime_robustness.csv")
yy=np.arange(len(rb))
ax.barh(yy,rb.gap,color=BLUE,edgecolor='k',linewidth=0.5,height=0.6)
for i,v in enumerate(rb.gap): ax.text(v+0.25,i,f"{v:.1f}",va='center',fontsize=7)
ax.set_yticks(yy); ax.set_yticklabels(rb.N_definition,fontsize=6.8); ax.invert_yaxis()
ax.set_xlabel("D' − N' gap (pts, length-matched)"); ax.set_xlim(0,11)
ax.set_title("C  Stable across definitions",loc='left',weight='bold')

ax=axes[3]
ax.bar([0,1],[Xv*100,Xp*100],color=[GREY,BLUE],edgecolor='k',linewidth=0.5,width=0.5)
for i,v in enumerate([Xv*100,Xp*100]): ax.text(i,v+0.06,f"{v:.2f}%",ha='center',fontsize=8,weight='bold')
ax.set_xticks([0,1]); ax.set_xticklabels(["X\npotential only","X'\n+ mouse\nconservation"],fontsize=7)
ax.set_ylabel("% of dual-coding-capable regions\nthat are NOT real DCRs")
ax.set_ylim(98.5,100.3)
ax.text(0.5,0.30,f"specificity gain\n{(1-Xp)/(1-Xv):.2f}×",transform=ax.transAxes,ha='center',
        fontsize=9,weight='bold',style='italic')
ax.set_title("D  Conservation barely narrows it",loc='left',weight='bold')
fig.suptitle("Dual-coding potential is necessary but far from sufficient — and conservation of potential adds only modest specificity",
             x=0.005,ha='left',fontsize=10,weight='bold')
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(f"{WS}/tierA/xprime_conservation.png",dpi=300); plt.close(fig)
print("saved")
~~~~

### [msg 443] bash · env: dcr-tierA · Checking non-DCR exon contamination

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
C="$D/output/02_Conservation"
echo "=== chain files in DOOVER/data ==="
ls -la "$D/data" | grep -i chain
echo ""
echo "=== CONTAMINATION CHECK: do 'non-DCR' exons overlap real DCR intervals? ==="
DCRBED="$D/output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed"
sort -k1,1 -k2,2n "$DCRBED" > /tmp/dcr.s.bed
for F in nonDCR_exon_DCR_gene_UNIQ nonDCR_exon_nonDCR_gene_UNIQ; do
  sort -k1,1 -k2,2n "$C/$F.bed" > /tmp/q.bed
  tot=$(wc -l < /tmp/q.bed)
  ov=$(bedtools intersect -a /tmp/q.bed -b /tmp/dcr.s.bed -u | wc -l)
  echo "$F: $ov / $tot overlap a real DCR interval"
done
echo ""
echo "=== annotation depth: transcripts per gene available? ==="
ls "$REF" | grep -i -E "gencode|ENSG|ENST|gtf|bed" | head -20
~~~~

### [msg 445] bash · env: dcr-tierA · Verifying UCSC liftOver binary

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
which liftOver && liftOver 2>&1 | head -3
echo "---"
ls -la <HOME>/Documents/liftOver 2>/dev/null
~~~~

### [msg 447] bash · env: dcr-tierA · Re-lifting with UCSC liftOver

~~~~bash
export PATH=<HOME>/Documents:/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
X="$WS/tierA/xprime"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
MM="<WHEELER_LAB_DRIVE>/2024/07-02-DCE_Conservation/mm39.fa"
CH="$REF/hg38ToMm39.over.chain"
# N intervals as a bed (from the python-written tsv)
awk 'BEGIN{OFS="\t"} NR>1 {print $1,$2,$3,"N"NR-2,0,$6}' "$X/N_latent_intervals.tsv" > "$X/N_human.bed"
wc -l < "$X/N_human.bed"
for S in dcr_regions:D N_human:N; do
  IN="${S%%:*}"; TAG="${S##*:}"
  liftOver -minMatch=0.10 "$X/$IN.bed" "$CH" "$X/${TAG}_ucsc.bed" "$X/${TAG}_ucsc.unmapped" 2>&1 | tail -2
  echo "$TAG lifted: $(wc -l < "$X/${TAG}_ucsc.bed") / $(wc -l < "$X/$IN.bed")"
  bedtools getfasta -fi "$MM" -bed "$X/${TAG}_ucsc.bed" -s -name -fo "$X/${TAG}_ucsc.fa"
  echo "$TAG seqs: $(grep -c '^>' "$X/${TAG}_ucsc.fa")"
done
~~~~

### [msg 449] python · env: dcr-tierA · Rescoring with UCSC lift

~~~~python
sD=read_fa2(f"{X}/D_ucsc.fa"); sN=read_fa2(f"{X}/N_ucsc.fa")
Du=pd.DataFrame([(h,len(s),is_dual_strict(s)) for h,s in sD.items()],columns=['h','len','strict'])
Nu=pd.DataFrame([(h,len(s),is_dual_strict(s)) for h,s in sN.items()],columns=['h','len','strict'])
# human lengths
Du['hlen']=[len(hD_map.get(h.split('::')[0],'')) for h in Du.h]
Nsub2=Nsub.reset_index(drop=True)
Nu['nid']=[int(h.split('::')[0][1:]) for h in Nu.h]
Nu['hlen']=Nu.nid.map(Nsub2.span)
Du['lb']=pd.cut(Du.hlen,bins); Nu['lb']=pd.cut(Nu.hlen,bins)
rd=Du.groupby('lb',observed=True).strict.mean()*100
rn=Nu.groupby('lb',observed=True).strict.mean()*100
wd=Du.groupby('lb',observed=True).strict.size()
common=[b for b in rd.index if b in rn.index and wd[b]>=20]
ww2=wd[common]/wd[common].sum()
Du_m=float((rd.reindex(common)*ww2).sum()); Nu_m=float((rn.reindex(common)*ww2).sum())
print("=== UCSC liftOver (matches your pipeline) ===")
print(f"D' raw           : {Du.strict.mean()*100:.1f}%   (n={len(Du)})   [pyliftover: 81.5%]")
print(f"N' raw           : {Nu.strict.mean()*100:.1f}%   (n={len(Nu)})   [pyliftover: 69.9%]")
print(f"D' length-matched: {Du_m:.1f}%")
print(f"N' length-matched: {Nu_m:.1f}%")
print(f"length-matched gap = {Du_m-Nu_m:.1f} pts   [pyliftover: 4.3 pts]")
~~~~
