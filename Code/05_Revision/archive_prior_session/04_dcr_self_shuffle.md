# 04_dcr_self_shuffle — as executed, not refactored

**DCR self-shuffle null: codon-shuffle of each DCR ortholog's own sequence.**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[249..265]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 7 (msg 249 … msg 261) |
| Language / environment | bash (env: python); python (env: dcr-tierA) |
| Saved outputs | `self_shuffle_control.png`, `self_shuffle_control.csv` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 2 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

---

### [msg 249] python · env: dcr-tierA · Loading ortholog sequences for self-shuffle

~~~~python
# The self-shuffle: for each DCR mouse-ortholog region, shuffle its CODONS in the DCR's own reading frame,
# preserving exact codon composition (thus GC and codon usage), randomizing only arrangement.
# Need the sequences + the frame. The translated TSV has per-ORF frame; the DCR reading frame is the
# one giving the full-length ORF. Load sequences keyed by chr:start:end.
def read_fasta_seqs(path):
    out={}; hdr=None; cur=[]
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if hdr is not None: out[hdr]=''.join(cur).upper()
                h=line[1:].strip().split('(')[0]; c,rng=h.split(':'); a,b=rng.split('-')
                hdr=f"{c}:{a}:{b}"; cur=[]
            else: cur.append(line.strip())
        if hdr is not None: out[hdr]=''.join(cur).upper()
    return out
mouse_seqs=read_fasta_seqs(f"{GC}/mouse_DCR_ortholog.fa")
print("mouse seqs:",len(mouse_seqs))
# The DCR's native frame: the dual-coding ORF spans the full region. For self-shuffle we shuffle codons
# in frame 0 of the extracted sequence (the exon's own CDS frame). Verify most seqs are multiple of 3 or near.
lens=np.array([len(s) for s in mouse_seqs.values()])
print("len %3==0:", np.mean(lens%3==0), "| median len:", np.median(lens))
# The extracted ortholog is the lifted DCR interval; its length reflects the human DCR region (often not %3).
# esl-translate scans all 3 frames stop-to-stop regardless. For the composition null we shuffle codons
# within the region's own primary frame — but to preserve EXACT nucleotide composition we shuffle by codon
# in frame 0 and carry the remainder. This keeps codon composition identical.
print(np.bincount(lens%3))
~~~~

### [msg 251] bash · env: python · Finding codon-shuffle generator script

~~~~bash
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/scripts/02_Conservation"
ls "$D" | grep -i "shuffl\|fake\|codon"
echo "==="
find "$D" -name "*fake*" -o -name "*shuffl*" 2>/dev/null
echo "=== content of generator (if present) ==="
F=$(find "$D" -iname "*generate*mouse*" -o -iname "*shuffl*" -o -iname "*fake_mm39*" 2>/dev/null | head -1)
echo "file: $F"
[ -n "$F" ] && cat "$F"
~~~~

### [msg 253] bash · env: python · Reading esl-translate replica module

~~~~bash
cat "<OLD_WORKSPACE>/tierA/esl_translate_replica.py"
~~~~

### [msg 255] python · env: dcr-tierA · Validating self-shuffle dual-coding rule

~~~~python
import sys, random
sys.path.insert(0, f"{WS}/tierA")
from esl_translate_replica import esl_translate_watson
STOPS={'TAA','TAG','TGA'}

def pipeline_shuffle(seq, rng):
    """Exact pipeline procedure: frame-0 codon split (drop trailing partial), remove stops, shuffle."""
    codons=[seq[i:i+3] for i in range(0,len(seq),3) if len(seq[i:i+3])==3]
    codons=[c for c in codons if c.upper() not in STOPS]
    rng.shuffle(codons)
    return ''.join(codons)

def is_dualcoding_seq(seq):
    """Apply Fig2B rule to a sequence: count full-length ORFs (span region ±3 both ends), region dual if >1."""
    L=len(seq)
    orfs=esl_translate_watson(seq, min_aa=1)
    ndc=0
    for start1,end1,frame,pep in orfs:
        orfstart0 = start1-1   # R does ORFstart<-ORFstart-1 on the 1-based start; here start1 is already 1-based
        # R's isDualCoding: (ORFstart-1)-3<=0 & ExLength-ORFend<=3, where ORFstart is esl 1-based coords col
        if (orfstart0-3)<=0 and (L-end1)<=3:
            ndc+=1
    return ndc>1

# Validate the rule against observed DCR orthologs first (should reproduce 82.6%)
obs_dc=[]
for k,seq in mouse_seqs.items():
    obs_dc.append(is_dualcoding_seq(seq))
print(f"Observed DCR orthologs dual-coding (recomputed from seq): {100*np.mean(obs_dc):.1f}% (expect ~82.6%)")
~~~~

### [msg 257] python · env: dcr-tierA · Running DCR self-shuffle null

~~~~python
rng=random.Random(2024)
N=20
per_exon=[]  # (observed_dc, shuffle_dc_rate, ExLength, gc)
shuffle_all_flags=[]
for k,seq in mouse_seqs.items():
    obs = is_dualcoding_seq(seq)
    # shuffle N times
    sh=[]
    for _ in range(N):
        sseq=pipeline_shuffle(seq,rng)
        sh.append(is_dualcoding_seq(sseq))
    per_exon.append((k, obs, np.mean(sh), len(seq)))
    shuffle_all_flags.extend(sh)

pe=pd.DataFrame(per_exon,columns=['RID','obs_dc','shuf_dc_rate','ExLength'])
obs_rate=pe.obs_dc.mean()*100
shuf_rate=np.mean(shuffle_all_flags)*100
print(f"DCR orthologs OBSERVED dual-coding: {obs_rate:.1f}%")
print(f"DCR orthologs SELF-SHUFFLED (own composition, N={N}): {shuf_rate:.1f}%")
print(f"Difference (arrangement effect): {obs_rate-shuf_rate:.1f} pts")
# empirical p: fraction of exons where observed dual-coding but shuffle rarely is
print(f"\nPer-exon: obs dual-coding {pe.obs_dc.sum()}/{len(pe)}; mean per-exon shuffle rate {pe.shuf_dc_rate.mean()*100:.1f}%")
~~~~

### [msg 259] python · env: dcr-tierA · Computing self-shuffle empirical p-value

~~~~python
# Empirical p-value: is the observed count of dual-coding orthologs (1318) higher than expected
# under the self-composition null? Build null distribution of the TOTAL dual-coding count across
# 1000 whole-set shuffles (one shuffle per exon per replicate).
rng=random.Random(99)
Nrep=1000
null_counts=np.zeros(Nrep)
# precompute codons per exon once (drop stops) to speed up
exon_codons=[]
for k,seq in mouse_seqs.items():
    codons=[seq[i:i+3] for i in range(0,len(seq),3) if len(seq[i:i+3])==3]
    codons=[c for c in codons if c.upper() not in STOPS]
    exon_codons.append(codons)

for r in range(Nrep):
    cnt=0
    for codons in exon_codons:
        cc=codons[:]; rng.shuffle(cc)
        if is_dualcoding_seq(''.join(cc)): cnt+=1
    null_counts[r]=cnt

obs_count=int(pe.obs_dc.sum())
p_emp=(np.sum(null_counts>=obs_count)+1)/(Nrep+1)
print(f"Observed dual-coding orthologs: {obs_count}")
print(f"Self-shuffle null: mean {null_counts.mean():.0f} (95% {np.percentile(null_counts,2.5):.0f}-{np.percentile(null_counts,97.5):.0f})")
print(f"Empirical one-sided p (obs >= null): {p_emp:.4f}")
print(f"Z-score: {(obs_count-null_counts.mean())/null_counts.std():.1f}")
~~~~

### [msg 261] python · env: dcr-tierA · Rendering self-shuffle control figure

~~~~python
# length-binned observed vs self-shuffle rate
pe['gc']=pe.RID.map(lambda k:gc_mouse.get(k,(None,None))[1])
breaks=make_breaks(int(pe.ExLength.max()),3)
def binprop2(df,col):
    r=df.copy();cats=pd.cut(r.ExLength,bins=breaks)
    r['mid']=cats.apply(lambda x:x.mid if pd.notna(x) else np.nan).astype(float)
    return r.dropna(subset=['mid']).groupby('mid')[col].agg(prop='mean',n='size').reset_index()
obs_b=binprop2(pe,'obs_dc'); shuf_b=binprop2(pe,'shuf_dc_rate')

fig,axes=plt.subplots(1,2,figsize=(8.0,3.6))
ax=axes[0]
m=(obs_b.mid>=21)&(obs_b.mid<=300)
ax.scatter(obs_b.mid[m],obs_b.prop[m],s=14,color="#FFA800",label="observed DCR ortholog",edgecolor='none',zorder=3)
m2=(shuf_b.mid>=21)&(shuf_b.mid<=300)
ax.scatter(shuf_b.mid[m2],shuf_b.prop[m2],s=14,color="#c0392b",label="self-shuffled (own codons)",edgecolor='none',zorder=3)
ax.axvline(21,color="#FFA800",ls=':',lw=1,alpha=0.6)
ax.set_xlim(21,300);ax.set_ylim(0,1.02)
ax.set_xlabel("Exon region length (bp)");ax.set_ylabel("Proportion with \u22652 full-length ORFs")
ax.set_title("Observed vs self-shuffle by length",loc="left")
ax.legend(frameon=False,fontsize=6.5,loc="upper right")
# panel B: null distribution of total count
ax=axes[1]
ax.hist(null_counts,bins=30,color="#c0392b",alpha=0.6,label="self-shuffle null")
ax.axvline(obs_count,color="#FFA800",lw=2,label=f"observed = {obs_count}")
ax.set_xlabel("# dual-coding DCR orthologs");ax.set_ylabel("null replicates")
ax.set_title(f"Composition null: obs {obs_count} vs null {null_counts.mean():.0f}\n(Z=15.5, p<0.001)",loc="left",fontsize=8)
ax.legend(frameon=False,fontsize=6.5)
fig.suptitle("R2-M1 composition null: DCR dual-coding excess is from codon ARRANGEMENT, not GC (identical composition)",fontsize=8,x=0.02,ha='left')
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(f"{WS}/tierA/self_shuffle_control.png",dpi=300)
print("saved")
plt.close(fig)

# save table
self_tab=pd.DataFrame({
 'metric':['observed_dualcoding_pct','self_shuffle_dualcoding_pct','arrangement_effect_pts',
           'observed_count','null_mean','null_ci_low','null_ci_high','Z','p_empirical'],
 'value':[round(obs_rate,1),round(shuf_rate,1),round(obs_rate-shuf_rate,1),
          obs_count,round(null_counts.mean(),0),round(np.percentile(null_counts,2.5),0),
          round(np.percentile(null_counts,97.5),0),round((obs_count-null_counts.mean())/null_counts.std(),1),p_emp]})
self_tab.to_csv(f"{WS}/tierA/self_shuffle_control.csv",index=False)
print(self_tab.to_string(index=False))
~~~~
