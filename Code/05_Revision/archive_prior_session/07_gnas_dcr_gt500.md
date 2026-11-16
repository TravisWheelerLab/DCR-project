# 07_gnas_dcr_gt500 — as executed, not refactored

**GNAS: DCRs > 500 nt and the mouse dual-coding ORFs of DCR 201.**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[303..330]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 10 (msg 303 … msg 322) |
| Language / environment | bash (env: python); python (env: dcr-tierA) |
| Saved outputs | `gnas_data.png`, `gnas_dcr_gt500.csv` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 4 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

---

### [msg 303] bash · env: python · Finding GNAS mentions in preprint

~~~~bash
FT="<OLD_WORKSPACE>/articles/10.64898_2026.02.24.707805.txt"
grep -n -i "gnas\|data not shown\|not shown" "$FT" | head -30
~~~~

### [msg 306] bash · env: python · Extracting GNAS context

~~~~bash
FT="<OLD_WORKSPACE>/articles/10.64898_2026.02.24.707805.txt"
grep -o -i -E ".{280}gnas.{280}" "$FT" | head -8
echo "===== data not shown context ====="
grep -o -i -E ".{220}(data not shown|not shown).{160}" "$FT" | head -8
~~~~

### [msg 308] python · env: dcr-tierA · Finding DCRs over 500nt

~~~~python
# DCR interval bed with gene + dceID
bed=pd.read_csv("<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed",
                sep="\t",header=None,dtype=str)
print("bed shape",bed.shape)
print(bed.head(3).to_string())
bed.columns=['chr','start','end','dceID','strand','gene','n'][:bed.shape[1]]
bed['start']=bed.start.astype(int);bed['end']=bed.end.astype(int)
bed['len']=bed.end-bed.start
big=bed[bed.len>500].copy()
print("\nDCR intervals >500nt:",len(big))
print(big[['chr','start','end','gene','len']].to_string(index=False))
print("\nunique genes with DCR>500nt:",sorted(big.gene.unique()))
~~~~

### [msg 310] python · env: dcr-tierA · Loading mouse conservation status

~~~~python
# DCRmiceStatus.tsv gives per-DCR mouse conservation status
mice=pd.read_csv(f"{WS}/tierA/inputs/DCRmiceStatus.tsv",sep="\t",header=None,dtype=str)
print("mice status shape",mice.shape); print(mice.head(3).to_string())
~~~~

### [msg 312] python · env: dcr-tierA · Mapping big DCRs to mouse status

~~~~python
mice=pd.read_csv(f"{WS}/tierA/inputs/DCRmiceStatus.tsv",sep=r"\s+",dtype=str)
mice.columns=['DCE_ID','mice_status','group']
print(mice.mice_status.value_counts().to_dict())
# map big DCRs to their dceID -> mouse status
big2=big.merge(mice,left_on='dceID',right_on='DCE_ID',how='left')
print("\nThe 3 DCR>500nt genes and their mouse conservation status:")
print(big2[['gene','chr','start','end','len','dceID','mice_status']].to_string(index=False))
~~~~

### [msg 314] python · env: dcr-tierA · Mapping big DCRs to mouse status

~~~~python
rows=[]
with open(f"{WS}/tierA/inputs/DCRmiceStatus.tsv") as fh:
    header=fh.readline()
    for line in fh:
        t=line.split()
        if len(t)>=3:
            # DCE_ID, mice_status may be multi-token (e.g. 'not DCR','too short','not mapped')
            did=t[0]; grp=t[-1]; status=" ".join(t[1:-1])
            rows.append((did,status,grp))
mice=pd.DataFrame(rows,columns=['DCE_ID','mice_status','group'])
print(mice.mice_status.value_counts().to_dict())
big2=big.merge(mice,left_on='dceID',right_on='DCE_ID',how='left')
print("\nThe 3 DCR>500nt genes and their mouse conservation status:")
print(big2[['gene','chr','start','end','len','dceID','mice_status']].to_string(index=False))
~~~~

### [msg 316] python · env: dcr-tierA · Loading human-to-mouse liftover bed

~~~~python
# Confirm GNAS mouse ortholog dual-coding directly. Find its ortholog sequence (dceID 201).
# The mouse ortholog FASTA is keyed by coordinate; need the GNAS mouse ortholog coords.
# The DCE_LO.to.mm39.bed maps human DCR -> mouse coords with dceID. Load it.
lo=pd.read_csv("<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation/DCE_overlapping_regions_LO.to.mm39.bed",
               sep="\t",header=None,dtype=str)
print("LO bed shape",lo.shape); print(lo.head(2).to_string())
~~~~

### [msg 318] python · env: dcr-tierA · Confirming GNAS mouse dual-coding ORFs

~~~~python
import re
lo.columns=['chr','start','end','info','score','strand']
lo['dceID']=lo.info.str.extract(r'DCRid=(\d+)')
gnas_lo=lo[lo.dceID=='201']
print("GNAS (dceID 201) mouse ortholog interval(s):")
print(gnas_lo.to_string(index=False))
# get its sequence + run the ORF scan
def read_fasta_seqs(path):
    out={};hdr=None;cur=[]
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if hdr is not None:out[hdr]=''.join(cur).upper()
                h=line[1:].strip().split('(')[0];c,rng=h.split(':');a,b=rng.split('-')
                hdr=f"{c}:{a}:{b}";cur=[]
            else:cur.append(line.strip())
        if hdr is not None:out[hdr]=''.join(cur).upper()
    return out
mouse_seqs=read_fasta_seqs(f"{GC}/mouse_DCR_ortholog.fa")
for _,r in gnas_lo.iterrows():
    key=f"{r.chr}:{r.start}:{r.end}"
    seq=mouse_seqs.get(key)
    if seq is None:
        print("  no seq for",key); continue
    orfs=esl_translate_watson(seq,min_aa=1)
    L=len(seq)
    full=[(s1,e1,fr,len(p)) for s1,e1,fr,p in orfs if (s1-1-3)<=0 and (L-e1)<=3]
    print(f"\n  GNAS mouse ortholog {key} len={L}nt")
    print(f"  total ORFs={len(orfs)}; full-length (span region) ORFs={len(full)} -> dual-coding={len(full)>1}")
    for s1,e1,fr,alen in full: print(f"    frame {fr}: nt {s1}-{e1}, {alen} aa")
~~~~

### [msg 320] python · env: dcr-tierA · Confirming GNAS mouse dual-coding ORFs

~~~~python
lo.columns=['chr','start','end','descr','score','strand']
lo['dceID']=lo['descr'].str.extract(r'DCRid=(\d+)')
gnas_lo=lo[lo.dceID=='201']
print("GNAS (dceID 201) mouse ortholog interval(s):")
print(gnas_lo.to_string(index=False))
def read_fasta_seqs(path):
    out={};hdr=None;cur=[]
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if hdr is not None:out[hdr]=''.join(cur).upper()
                h=line[1:].strip().split('(')[0];c,rng=h.split(':');a,b=rng.split('-')
                hdr=f"{c}:{a}:{b}";cur=[]
            else:cur.append(line.strip())
        if hdr is not None:out[hdr]=''.join(cur).upper()
    return out
mouse_seqs=read_fasta_seqs(f"{GC}/mouse_DCR_ortholog.fa")
gnas_orf_info=[]
for _,r in gnas_lo.iterrows():
    key=f"{r.chr}:{r.start}:{r.end}"; seq=mouse_seqs.get(key)
    if seq is None: print("  no seq for",key); continue
    orfs=esl_translate_watson(seq,min_aa=1); L=len(seq)
    full=[(s1,e1,fr,len(p)) for s1,e1,fr,p in orfs if (s1-1-3)<=0 and (L-e1)<=3]
    print(f"\n  GNAS mouse ortholog {key} len={L}nt: total ORFs={len(orfs)}; full-length ORFs={len(full)} -> dual-coding={len(full)>1}")
    for s1,e1,fr,alen in full: print(f"    frame {fr}: nt {s1}-{e1}, {alen} aa")
    gnas_orf_info.append((key,L,len(orfs),len(full),full))
~~~~

### [msg 322] python · env: dcr-tierA · Rendering GNAS supporting figure

~~~~python
# Build the GNAS supporting figure: (A) table of the 3 DCR>500nt genes + mouse status;
# (B) GNAS mouse ortholog dual-ORF architecture schematic.
fig,axes=plt.subplots(1,2,figsize=(8.4,3.2),gridspec_kw={'width_ratios':[1.05,1]})
# Panel A: table
ax=axes[0];ax.axis('off')
tbl=big2[['gene','len','mice_status']].copy()
tbl.columns=['gene','DCR len (nt)','mouse status']
cellcols=[['#FFF6E0']*3 if s=='DCR' else ['#f2f2f2']*3 for s in big2.mice_status]
t=ax.table(cellText=tbl.values,colLabels=tbl.columns,loc='center',cellLoc='center',cellColours=cellcols)
t.auto_set_font_size(False);t.set_fontsize(8);t.scale(1,1.8)
for (row,col),cell in t.get_celld().items():
    if row==0: cell.set_text_props(weight='bold');cell.set_facecolor('#FFA800')
ax.set_title("3 human genes with a DCR > 500 nt",loc='left',fontsize=8.5,y=0.86)
ax.text(0.5,0.14,"Only GNAS conserved as dual-coding in mouse",ha='center',fontsize=7.5,style='italic',transform=ax.transAxes)
# Panel B: GNAS mouse ortholog ORF architecture
ax=axes[1]
key,L,ntot,nfull,full=gnas_orf_info[0]
ax.add_patch(plt.Rectangle((0,2.4),L,0.5,color='#dddddd',ec='k',lw=0.5))
ax.text(L/2,2.9,f"GNAS mouse ortholog region ({L} nt)",ha='center',va='bottom',fontsize=7)
ycols={1:'#FFA800',3:'#c0392b'}
for i,(s1,e1,fr,alen) in enumerate(full):
    y=1.4-i*0.9
    ax.add_patch(plt.Rectangle((s1-1,y),e1-s1+1,0.5,color=ycols.get(fr,'#2c7fb8'),ec='k',lw=0.5))
    ax.text((s1-1+e1)/2,y+0.25,f"frame {fr} ORF: {alen} aa",ha='center',va='center',fontsize=7,color='white',weight='bold')
ax.set_xlim(-50,L+50);ax.set_ylim(0,3.4);ax.axis('off')
ax.set_title("GNAS: two full-length ORFs in distinct frames\n(conserved dual-coding, R1-m1 'data not shown')",loc='left',fontsize=8)
fig.suptitle("R1-m1: substantiating the GNAS conserved-dual-coding claim",fontsize=9,x=0.02,ha='left')
fig.tight_layout(rect=[0,0,1,0.94]);fig.savefig(f"{WS}/tierA/gnas_data.png",dpi=300)
print("saved");plt.close(fig)
# save table
big2out=big2[['gene','chr','start','end','len','dceID','mice_status']].copy()
big2out['mouse_full_length_ORFs']=[nfull if g=='gnas' else '' for g in big2out.gene]
big2out.to_csv(f"{WS}/tierA/gnas_dcr_gt500.csv",index=False)
print(big2out.to_string(index=False))
~~~~
