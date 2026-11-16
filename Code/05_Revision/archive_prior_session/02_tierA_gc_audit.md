# 02_tierA_gc_audit — as executed, not refactored

**Tier A GC audit: per-exon GC of DCR vs non-DCR exons (human) and mouse orthologs vs random control.**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[210..229]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 7 (msg 212 … msg 225) |
| Language / environment | bash (env: dcr-tierA); bash (env: python); python (env: dcr-tierA) |
| Saved outputs | `gc_audit.png`, `gc_audit_summary.csv`, `gc_audit_tests.csv` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 7 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

---

### [msg 212] bash · env: python · Inspecting exon bed formats

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
C="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
echo "=== bed formats (human exon classes) ==="
for f in DCR_exons_DCR_genes_UNIQ.bed nonDCR_exon_DCR_gene_UNIQ.bed nonDCR_exon_nonDCR_gene_UNIQ.bed; do
  echo "--- $f ($(wc -l < "$C/$f") lines) ---"; head -2 "$C/$f"
done
echo ""
echo "hg38.fa present: $(ls -lh "$REF/hg38.fa" | awk '{print $5}')"
echo "mouse ortholog fa header: $(head -1 "$C/DCE_overlapping_regions_LO.to.mm39.fa")"
echo "random ctrl fa: $(head -1 "$C/fake_DCE_intervals.fa")"
~~~~

### [msg 214] bash · env: dcr-tierA · Extracting human exon sequences via bedtools

~~~~bash
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
WS="<OLD_WORKSPACE>"
C="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation"
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
GC="$WS/tierA/gc"; mkdir -p "$GC"
# Extract first 3 cols to clean bed (getfasta needs valid bed), then getfasta
for f in DCR_exons_DCR_genes_UNIQ nonDCR_exon_DCR_gene_UNIQ nonDCR_exon_nonDCR_gene_UNIQ; do
  cut -f1-3 "$C/$f.bed" > "$GC/$f.b3.bed"
  bedtools getfasta -fi "$REF/hg38.fa" -bed "$GC/$f.b3.bed" -fo "$GC/$f.fa" 2>/dev/null
  echo "$f: $(grep -c '>' "$GC/$f.fa") seqs"
done
# copy the mouse ortholog + random control FASTAs locally too
cp "$C/DCE_overlapping_regions_LO.to.mm39.fa" "$GC/mouse_DCR_ortholog.fa"
cp "$C/fake_DCE_intervals.fa" "$GC/random_ctrl.fa"
echo "mouse_DCR_ortholog: $(grep -c '>' "$GC/mouse_DCR_ortholog.fa") seqs"
echo "random_ctrl: $(grep -c '>' "$GC/random_ctrl.fa") seqs"
~~~~

### [msg 216] python · env: dcr-tierA · Computing per-class GC content

~~~~python
import numpy as np, pandas as pd
GC=f"{WS}/tierA/gc"
def read_fasta_gc(path):
    lens=[]; gcs=[]; seq=[]
    def flush(s):
        if s:
            s=s.upper(); L=len(s)
            if L>0:
                g=s.count('G')+s.count('C'); n=g+s.count('A')+s.count('T')
                if n>0: lens.append(L); gcs.append(100.0*g/n)
    with open(path) as fh:
        cur=[]
        for line in fh:
            if line.startswith('>'):
                flush(''.join(cur)); cur=[]
            else: cur.append(line.strip())
        flush(''.join(cur))
    return np.array(lens), np.array(gcs)

classes={
 "DCR_exon_human":         f"{GC}/DCR_exons_DCR_genes_UNIQ.fa",
 "nonDCR_exon_DCR_gene":   f"{GC}/nonDCR_exon_DCR_gene_UNIQ.fa",
 "nonDCR_exon_nonDCR_gene":f"{GC}/nonDCR_exon_nonDCR_gene_UNIQ.fa",
 "mouse_DCR_ortholog":     f"{GC}/mouse_DCR_ortholog.fa",
 "random_ctrl_mouse":      f"{GC}/random_ctrl.fa",
}
gcdata={}
rows=[]
for name,path in classes.items():
    L,G=read_fasta_gc(path); gcdata[name]=(L,G)
    rows.append(dict(cls=name,n=len(G),gc_median=round(np.median(G),2),gc_mean=round(G.mean(),2),
                     gc_q25=round(np.percentile(G,25),2),gc_q75=round(np.percentile(G,75),2),
                     len_median=int(np.median(L))))
gc_summary=pd.DataFrame(rows)
print(gc_summary.to_string(index=False))
~~~~

### [msg 218] python · env: dcr-tierA · Computing per-class GC content

~~~~python
import numpy as np, pandas as pd
WS="<OLD_WORKSPACE>"
GC=f"{WS}/tierA/gc"
def read_fasta_gc(path):
    lens=[]; gcs=[]
    def flush(s):
        if s:
            s=s.upper(); L=len(s)
            if L>0:
                g=s.count('G')+s.count('C'); n=g+s.count('A')+s.count('T')
                if n>0: lens.append(L); gcs.append(100.0*g/n)
    with open(path) as fh:
        cur=[]
        for line in fh:
            if line.startswith('>'): flush(''.join(cur)); cur=[]
            else: cur.append(line.strip())
        flush(''.join(cur))
    return np.array(lens), np.array(gcs)

classes={
 "DCR_exon_human":         f"{GC}/DCR_exons_DCR_genes_UNIQ.fa",
 "nonDCR_exon_DCR_gene":   f"{GC}/nonDCR_exon_DCR_gene_UNIQ.fa",
 "nonDCR_exon_nonDCR_gene":f"{GC}/nonDCR_exon_nonDCR_gene_UNIQ.fa",
 "mouse_DCR_ortholog":     f"{GC}/mouse_DCR_ortholog.fa",
 "random_ctrl_mouse":      f"{GC}/random_ctrl.fa",
}
gcdata={}; rows=[]
for name,path in classes.items():
    L,G=read_fasta_gc(path); gcdata[name]=(L,G)
    rows.append(dict(cls=name,n=len(G),gc_median=round(np.median(G),2),gc_mean=round(G.mean(),2),
                     gc_q25=round(np.percentile(G,25),2),gc_q75=round(np.percentile(G,75),2),
                     len_median=int(np.median(L))))
gc_summary=pd.DataFrame(rows)
print(gc_summary.to_string(index=False))
~~~~

### [msg 220] python · env: dcr-tierA · Running GC statistical tests

~~~~python
from scipy import stats
def rank_biserial(a,b):
    # from MWU: r = 1 - 2U/(n1*n2)
    u,p=stats.mannwhitneyu(a,b,alternative='two-sided')
    n1,n2=len(a),len(b)
    r=1-2*u/(n1*n2)
    return u,p,r

comps=[("DCR_exon_human","nonDCR_exon_DCR_gene"),
       ("DCR_exon_human","nonDCR_exon_nonDCR_gene"),
       ("mouse_DCR_ortholog","random_ctrl_mouse")]
res=[]
for a,b in comps:
    Ga=gcdata[a][1]; Gb=gcdata[b][1]
    u,p,r=rank_biserial(Ga,Gb)
    res.append(dict(group_A=a,group_B=b,
                    median_A=round(np.median(Ga),2),median_B=round(np.median(Gb),2),
                    median_diff=round(np.median(Ga)-np.median(Gb),2),
                    MWU_p=f"{p:.2e}",rank_biserial=round(r,3)))
gc_tests=pd.DataFrame(res)
print(gc_tests.to_string(index=False))
gc_summary.to_csv(f"{WS}/tierA/gc_audit_summary.csv",index=False)
gc_tests.to_csv(f"{WS}/tierA/gc_audit_tests.csv",index=False)
print("\nsaved gc_audit_summary.csv + gc_audit_tests.csv")
~~~~

### [msg 222] python · env: dcr-tierA · Rendering GC audit figure

~~~~python
import matplotlib as mpl, matplotlib.pyplot as plt
from figure_style_kernel import apply_figure_style
apply_figure_style()

fig, axes = plt.subplots(1,2,figsize=(7.2,3.4))
# panel A: human classes
human=[("DCR_exon_human","#FFA800","DCR exon"),
       ("nonDCR_exon_DCR_gene","#d9a441","non-DCR exon,\nDCR gene"),
       ("nonDCR_exon_nonDCR_gene","#9a9a9a","non-DCR exon,\nnon-DCR gene")]
ax=axes[0]
bins=np.linspace(20,90,50)
for name,col,lab in human:
    G=gcdata[name][1]
    ax.hist(G,bins=bins,density=True,histtype='step',color=col,lw=1.6,label=lab)
    ax.axvline(np.median(G),color=col,ls=':',lw=1,alpha=0.8)
ax.set_xlabel("GC content (%)"); ax.set_ylabel("density")
ax.set_title("Human exons",loc="left",fontsize=8)
ax.legend(frameon=False,fontsize=6)
# panel B: mouse DCR ortholog vs random control
ax=axes[1]
mouse=[("mouse_DCR_ortholog","#FFA800","DCR mouse ortholog"),
       ("random_ctrl_mouse","#9a9a9a","random control")]
for name,col,lab in mouse:
    G=gcdata[name][1]
    ax.hist(G,bins=bins,density=True,histtype='step',color=col,lw=1.6,label=lab)
    ax.axvline(np.median(G),color=col,ls=':',lw=1,alpha=0.8)
ax.set_xlabel("GC content (%)"); ax.set_ylabel("density")
ax.set_title("Mouse exons",loc="left",fontsize=8)
ax.legend(frameon=False,fontsize=6.5)
fig.suptitle("GC audit (R2-M1): DCR exons modestly GC-enriched vs controls (+2.8–3.0 pts human, +1.2 pts mouse)",
             fontsize=8.5,x=0.02,ha='left')
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig(f"{WS}/tierA/gc_audit.png",dpi=300)
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
ov=[(a.get_text(),b.get_text()) for i,(a,ba) in enumerate(texts) for b,bb in texts[i+1:] if ba.overlaps(bb)]
print("overlaps:",len(ov))
plt.close(fig)
~~~~

### [msg 225] python · env: dcr-tierA · Rendering GC audit figure

~~~~python
import matplotlib as mpl, matplotlib.pyplot as plt
mpl.rcParams.update({'font.size':8,'axes.titlesize':8,'axes.labelsize':8,'xtick.labelsize':7,
 'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':110,
 'font.family':'sans-serif','axes.linewidth':0.8,'legend.fontsize':6.5})

fig, axes = plt.subplots(1,2,figsize=(7.2,3.4))
bins=np.linspace(20,90,50)
human=[("DCR_exon_human","#FFA800","DCR exon"),
       ("nonDCR_exon_DCR_gene","#d9a441","non-DCR exon, DCR gene"),
       ("nonDCR_exon_nonDCR_gene","#9a9a9a","non-DCR exon, non-DCR gene")]
ax=axes[0]
for name,col,lab in human:
    G=gcdata[name][1]
    ax.hist(G,bins=bins,density=True,histtype='step',color=col,lw=1.6,label=lab)
    ax.axvline(np.median(G),color=col,ls=':',lw=1,alpha=0.8)
ax.set_xlabel("GC content (%)"); ax.set_ylabel("density")
ax.set_title("Human exons",loc="left"); ax.legend(frameon=False,fontsize=6)
ax=axes[1]
for name,col,lab in [("mouse_DCR_ortholog","#FFA800","DCR mouse ortholog"),("random_ctrl_mouse","#9a9a9a","random control")]:
    G=gcdata[name][1]
    ax.hist(G,bins=bins,density=True,histtype='step',color=col,lw=1.6,label=lab)
    ax.axvline(np.median(G),color=col,ls=':',lw=1,alpha=0.8)
ax.set_xlabel("GC content (%)"); ax.set_ylabel("density")
ax.set_title("Mouse exons",loc="left"); ax.legend(frameon=False,fontsize=6.5)
fig.suptitle("GC audit (R2-M1): DCR exons modestly GC-enriched (+2.8–3.0 pts human, +1.2 pts mouse)",fontsize=8.5,x=0.02,ha='left')
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig(f"{WS}/tierA/gc_audit.png",dpi=300)
print("saved gc_audit.png")
plt.close(fig)
~~~~
