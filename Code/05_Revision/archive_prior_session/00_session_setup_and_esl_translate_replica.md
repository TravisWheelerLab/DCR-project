# 00_session_setup_and_esl_translate_replica — as executed, not refactored

**Session setup: reading the reviewer letter, repo review, toolchain discovery, construction and validation of the pure-Python esl-translate replica, copying of inputs.**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[3..143]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 42 (msg 3 … msg 137) |
| Language / environment | bash (env: python); python (env: python) |
| Saved outputs | `DCR_revision_plan.md`, `DCR_reviewer_triage.csv`, `repo_review_notes.md`, `reference_data_manifest.tsv`, `esl_translate_replica.py` |

This file is a verbatim log of the code cells that were executed, in order, during the July–August 2026
revision session (run under a previous account). It is **not** a runnable script: cells relied on
interactive kernel state carried over from earlier cells (and from earlier files in this directory),
`python` cells ran in a persistent Jupyter-like kernel, `bash` cells in a shell, and `repl` cells in a
control-plane kernel whose `host` object queried the platform's artifact store. Nothing was edited except:

* absolute paths were replaced by placeholders — `<OLD_WORKSPACE>` (the old account's private
  workspace), `<OLD_ARTIFACT_STORE>/` (its artifact store; the file name after it is the artifact
  name), `<WHEELER_LAB_DRIVE>` (the lab Google Drive folder that holds `2026/02_04_DCR-merge-overlap/DOOVER`
  and `2025/06-17_DCR_pipeline/DCR-pipeline/data`), `<HOME>`; 29 substitution(s) in this file;
* e-mail addresses would be replaced by `<EMAIL>`; a credential scan found nothing.

Environment operations in this span (not code cells):

* `[msg 129] manage_packages: {"mode": "install", "environment": "python", "packages": ["pyBigWig"], "use_pip": true}`
* `[msg 131] manage_environments: {"mode": "create", "name": "dcr-tierA", "packages": ["python=3.11", "pybigwig", "numpy", "pandas", "scipy", "matplotlib", "bedtools", "samtools"], "channels": ["bioconda", "conda-forge"]}`

---

### [msg 3] bash · env: python · Checking artifact path resolution

~~~~bash
cd /tmp && python3 -c "
import zipfile, re
p = '{{artifact:a8716afc-8245-4c0a-8467-66e4e6ae17ad}}'
print(p)
"
~~~~

### [msg 5] python · env: python · Opening response-to-reviewers docx

~~~~python
import zipfile, re
from xml.etree import ElementTree as ET

path = host.artifact_path("a8716afc-8245-4c0a-8467-66e4e6ae17ad")
z = zipfile.ZipFile(path)
print(z.namelist()[:20])
~~~~

### [msg 7] python · env: python · Extracting reviewer text with formatting

~~~~python
import zipfile
from xml.etree import ElementTree as ET

path = host.artifact_path("a8716afc-8245-4c0a-8467-66e4e6ae17ad")
z = zipfile.ZipFile(path)
xml = z.read('word/document.xml').decode('utf-8')

ns = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
root = ET.fromstring(xml)
body = root.find('w:body', ns)

def para_runs(p):
    runs = []
    for r in p.findall('w:r', ns):
        rpr = r.find('w:rPr', ns)
        color = None; bold=False; italic=False
        if rpr is not None:
            c = rpr.find('w:color', ns)
            if c is not None: color = c.get('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}val')
            if rpr.find('w:b', ns) is not None: bold=True
            if rpr.find('w:i', ns) is not None: italic=True
        txt = ''.join(t.text or '' for t in r.findall('w:t', ns))
        if txt:
            runs.append((txt,color,bold,italic))
    return runs

lines=[]
for p in body.findall('w:p', ns):
    runs = para_runs(p)
    if not runs:
        lines.append('')
        continue
    # tag runs by formatting
    seg=[]
    for txt,color,bold,italic in runs:
        tag=''
        if color and color.lower() not in ('000000','auto','ffffff'):
            tag+=f'[{color}]'
        if bold: tag='**'+tag
        s=txt
        if bold: s=f'**{s}**'
        if color and color.lower() not in ('000000','auto','ffffff'):
            s=f'⟨{color}⟩{s}⟨/⟩'
        seg.append(s)
    lines.append(''.join(seg))

full='\n'.join(lines)
print(len(full))
open('/tmp/reviewers.txt','w').write(full)
print(full[:3000])
~~~~

### [msg 11] python · env: python · Reading reviewer text part 2

~~~~python
print(full[3000:9000])
~~~~

### [msg 13] python · env: python · Reading reviewer text part 3

~~~~python
print(full[9000:15327])
~~~~

### [msg 21] python · env: python · Building reviewer triage matrix

~~~~python
import csv

rows = [
# id, reviewer, short, category, effort, priority, response
["AE","Editor","Validity of central conclusions (functional selection on dual coding)","Reframe + reanalysis","—","Critical","Addressed via R2-1/R1-2/R1-3 controls + cautious reframing of 'conserved potential = functional'"],
["AE","Editor","Uncited literature questioning originality","Rebuttal + text","Low","Critical","Add evo-AS context refs (Bénitière24, Jiang25, Xu19, Radrizzani24); add explicit 'advances vs prior DCR work' paragraph"],

["R1-M1","Ref1","Compare DCRs to 'regular' AS (conservation, NMD, disorder, adaptive?)","New analysis (reuses pipeline)","Med","High","THE unifying experiment: matched regular-AS event set vs DCRs on NMD %, disorder, tissue usage; discuss vs drift-barrier"],
["R1-M2a","Ref1","DCR potential of non-DCR human genes (necessary vs sufficient)","Reanalysis","Med","Critical","Scan non-DCR human exons for latent DCR potential; report fraction; contextualize what 'conserved potential' means"],
["R1-M2b","Ref1","NMD could mask detrimental isoforms → potential neutral (Jiang25)","Rebuttal + discussion","Low","High","Acknowledge; note RUST hypothesis already accommodates this; soften functional claim"],
["R1-M3","Ref1","Coding constraint on one frame conserves the other (byproduct)","Rebuttal (existing control) + caveat","Low","Critical","Rebut: non-DCR control exons are equally coding-constrained yet retain potential less often (once size+GC matched); add caveat"],
["R1-m1","Ref1","'data not shown' (GNAS), line 181-182","Text + figure","Low","Easy","Show GNAS result (supp fig or remove phrase)"],
["R1-m2","Ref1","Discuss alt transcription initiation as noise (Xu19)","Discussion","Low","Easy","Add sentence + cite"],
["R1-m3","Ref1","DCRs may mediate selection on synonymous sites (Radrizzani24)","Discussion","Low","High","Add; aligns with AE's own field (unwanted-transcript hypothesis)"],

["R2-0","Ref2","Emphasize advantages over prior work","Text","Low","High","Add 'what's new' paragraph (SwissProt-curated, protein-level, genome-wide architecture, GTEx frame quant)"],
["R2-M1","Ref2","GC-content confound in DCR-potential conservation","Reanalysis","Med","Critical","Report GC of DCR vs control exons; GC-match controls; add codon-shuffle of DCR exons themselves"],
["R2-M2","Ref2","'fewer mouse isoforms' argument invalid; still compare DCR vs non-DCR","Reframe + reanalysis","Low","High","Reframe: comparison IS to non-DCR (control), not to mouse isoform catalog; clarify text + apply GC control"],
["R2-M3","Ref2","Actual frame-usage frequency via matched human/mouse transcriptomes (Cardoso-Moreira19)","New data (optional)","High","Medium","Compare non-canonical frame freq to matched non-DCR AS events; overlaps R1-M1; do if reviewer insists"],
["R2-M4a","Ref2","GC control on phyloP conservation","Reanalysis","Low","High","Add GC covariate/matching to phyloP comparison"],
["R2-M4b","Ref2","Splice-site conservation for DCR-generating alt splice sites vs non-DCR","New analysis","Med","Medium","phyloP on splice-site motifs; compare DCR vs non-DCR alt splice sites"],
["R2-M5","Ref2","% DCR genes across expression classes (Bénitière link)","Reanalysis","Low","High","Stratify DCR incidence by gene expression (GTEx) vs all SwissProt; test AS-cost prediction"],
["R2-M6","Ref2","Are DCR exons more tissue-specific than non-DCR? (Fig 4 needs reference)","Reanalysis","Med","High","Compute tissue-specificity (tau) for DCR vs matched non-DCR exons"],
["R2-m1","Ref2","GABBR1 example unclear (splice variant? which exons?)","Text","Low","Easy","Clarify it's an AS variant; name exons/frames"],
["R2-m2","Ref2","How was isoform frequency estimated? (line 570-571)","Text","Low","Easy","State method (junction/EEJ counts)"],
["R2-m3","Ref2","Typos (Boneferroni, Peasron)","Text","Low","Easy","Proofread"],
["R2-m4","Ref2","'noisy transcription' → 'noisy splicing' (line 49 + significance stmt)","Text","Low","Easy","Fix wording throughout"],
]

with open("DCR_reviewer_triage.csv","w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["ID","Reviewer","Comment","Category","Effort","Priority","Proposed response"])
    w.writerows(rows)

from collections import Counter
print("rows:",len(rows))
print("by priority:",Counter(r[5] for r in rows))
print("by category:",Counter(r[3] for r in rows))
~~~~

### [msg 31] bash · env: python · Cloning DCR-project repo

~~~~bash
cd /tmp && rm -rf DCR-project && git clone --depth 1 https://github.com/TravisWheelerLab/DCR-project.git 2>&1 | tail -5 && echo "---CLONE DONE---" && cd DCR-project && echo "=== TOP LEVEL ===" && ls -la && echo "=== TREE (2 levels) ===" && find . -maxdepth 2 -not -path './.git/*' | sort | head -100
~~~~

### [msg 33] bash · env: python · Cloning repo into workspace

~~~~bash
cd "<OLD_WORKSPACE>" 2>/dev/null || cd .
pwd
rm -rf DCR-project
git clone --depth 1 https://github.com/TravisWheelerLab/DCR-project.git 2>&1 | tail -5
echo "---"
ls -la DCR-project 2>&1 | head
~~~~

### [msg 35] bash · env: python · Reading README and mapping Code/Data trees

~~~~bash
cd DCR-project
echo "===== README ====="
cat README.md
echo ""
echo "===== CODE TREE ====="
find Code -type f | sort
echo ""
echo "===== DATA TREE (dirs + file counts) ====="
find Data -type d | sort | while read d; do n=$(find "$d" -maxdepth 1 -type f | wc -l | tr -d ' '); echo "$d  [$n files]"; done
~~~~

### [msg 37] bash · env: python · Reading conservation pipeline and archive contents

~~~~bash
cd DCR-project
echo "===== 02_Conservation.sh ====="
cat Code/02_Conservation.sh
echo ""
echo "========== DATA FILES (top level) =========="
ls -la Data/
echo ""
echo "========== archive contents (no extract) =========="
for f in Data/*.tar.gz; do echo "--- $f ---"; tar tzf "$f" 2>/dev/null | head -20; echo "  ...($(tar tzf "$f" 2>/dev/null | wc -l | tr -d ' ') entries)"; done
~~~~

### [msg 42] bash · env: python · Reading control-generation scripts

~~~~bash
cd DCR-project
for s in random_CE_regions.py generate_fake_mouse_exons.py extract_codon_phased_CDS.py; do
echo "########## Code/02_Conservation/$s ##########"
cat "Code/02_Conservation/$s"
echo ""
done
~~~~

### [msg 44] bash · env: python · Inspecting data archives and Table S1

~~~~bash
cd DCR-project
echo "===== Data/ top-level files ====="
ls -la Data/ | grep -v '^d'
echo ""
echo "===== archive listings ====="
for f in Data/Extracted-DCEs_06202025.tar.gz Data/Species-MSAs.tar.gz Data/Mice.tar.gz Data/CIFs.tar.gz; do
  echo "--- $f ($(du -h "$f" | cut -f1)) ---"
  tar tzf "$f" 2>/dev/null | head -25
  echo "  total entries: $(tar tzf "$f" 2>/dev/null | wc -l | tr -d ' ')"
  echo ""
done
echo "===== Table_S1 head ====="
head -3 Data/Table_S1.tsv
echo "cols:"; head -1 Data/Table_S1.tsv | tr '\t' '\n' | cat -n
echo "rows: $(wc -l < Data/Table_S1.tsv)"
~~~~

### [msg 46] bash · env: python · Checking LFS and DCR bed provenance

~~~~bash
cd DCR-project
echo "===== is Mice.tar.gz an LFS pointer? ====="
file Data/Mice.tar.gz 2>/dev/null; head -c 200 Data/Mice.tar.gz 2>/dev/null; echo ""
echo "===== .gitattributes / LFS ====="
cat .gitattributes 2>/dev/null || echo "(no .gitattributes)"
echo ""
echo "===== does any shipped file already contain the DCR bed? ====="
grep -rl "filtered_DCR_overlapping_regions" Data/ 2>/dev/null | head
tar tzf Data/Extracted-DCEs_06202025.tar.gz 2>/dev/null | grep -iE 'bed|region|filtered' | head
echo ""
echo "===== 01_Description.sh : first 90 lines ====="
head -90 Code/01_Description.sh
~~~~

### [msg 48] bash · env: python · Checking DCR bed reconstruction source

~~~~bash
cd DCR-project
echo "===== git lfs pointers? ====="
git log --oneline -1 2>/dev/null
git lfs ls-files 2>/dev/null | head || echo "(no lfs)"
echo ""
echo "===== extract homo_sapiens.csv header + sample (genomic coords?) ====="
mkdir -p /tmp/dce_peek && tar xzf Data/Extracted-DCEs_06202025.tar.gz -C /tmp/dce_peek Extracted-DCEs_06202025/homo_sapiens.csv 2>/dev/null
head -3 /tmp/dce_peek/Extracted-DCEs_06202025/homo_sapiens.csv
echo "..."
echo "cols in homo_sapiens.csv:"; head -1 /tmp/dce_peek/Extracted-DCEs_06202025/homo_sapiens.csv | tr ',' '\n' | cat -n | head -30
echo ""
echo "===== a per-gene .out file (does it carry chr coords?) ====="
tar xzf Data/Extracted-DCEs_06202025.tar.gz -C /tmp/dce_peek Extracted-DCEs_06202025/homo_sapiens/519.mfsd1.out 2>/dev/null
head -40 /tmp/dce_peek/Extracted-DCEs_06202025/homo_sapiens/519.mfsd1.out
~~~~

### [msg 64] bash · env: python · Mapping DOOVER top-level structure

~~~~bash
cd "<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
echo "===== TOP LEVEL ====="
ls -la
echo ""
echo "===== TOP-LEVEL DIRS (1 level) ====="
find . -maxdepth 1 -type d | sort
echo ""
echo "===== the key bed exists? ====="
ls -la output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed 2>&1
wc -l output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed 2>&1
~~~~

### [msg 67] bash · env: python · Listing output directory tree

~~~~bash
cd "<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
echo "===== output/ tree ====="
find output -type f 2>/dev/null | sort
echo ""
echo "===== output/02_Conservation contents w/ sizes ====="
ls -la output/02_Conservation/ 2>/dev/null
~~~~

### [msg 70] bash · env: python · Listing reference data and scripts

~~~~bash
cd "<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
echo "===== data/ (reference files) ====="
ls -la data/ 2>/dev/null
echo ""
echo "===== big files across data/ (>50MB) ====="
find data -type f -size +50M -exec ls -lh {} \; 2>/dev/null | awk '{print $5"\t"$9}'
echo ""
echo "===== scripts tree ====="
find scripts -type f 2>/dev/null | sort
~~~~

### [msg 73] bash · env: python · Testing symlink resolution and key files

~~~~bash
cd "<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
echo "===== do the symlinked references resolve? ====="
for f in hg38.fa mm39.fa hg38.phyloP100way.bw hg38ToMm39.over.chain hg38.GENECODE.V48.bed correspondance_ENSG_ENST_v114.hg38.14; do
  if [ -e "data/$f" ]; then sz=$(stat -f%z "data/$f" 2>/dev/null); echo "OK   data/$f  ($sz bytes)"; else echo "DEAD data/$f  -> $(readlink "data/$f")"; fi
done
echo ""
echo "===== phyloP per-DCE score file (head) ====="
head -3 output/02_Conservation/meanPP_per_DCE.tsv
echo "rows: $(wc -l < output/02_Conservation/meanPP_per_DCE.tsv)"
echo ""
echo "===== DCRmiceStatus (the Fig2A/2B conservation calls) ====="
head -3 output/02_Conservation/DCRmiceStatus.tsv
echo "cols:"; head -1 output/02_Conservation/DCRmiceStatus.tsv | tr '\t' '\n' | cat -n
echo "rows: $(wc -l < output/02_Conservation/DCRmiceStatus.tsv)"
~~~~

### [msg 78] bash · env: python · Verifying reference files resolve and read

~~~~bash
REFDIR="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
echo "===== reference files present & materialized? ====="
for f in hg38.fa hg38.fa.fai mm39.fa mm39.fa.fai hg38.phyloP100way.bw hg38ToMm39.over.chain hg38.GENECODE.V48.bed correspondance_ENSG_ENST_v114.hg38.14; do
  p="$REFDIR/$f"
  if [ -f "$p" ]; then printf "OK   %-40s %s\n" "$f" "$(ls -lh "$p" | awk '{print $5}')"; else echo "MISS $f"; fi
done
echo ""
echo "===== can we actually READ bytes (not cloud stub)? test faidx + bigwig header ====="
head -c 60 "$REFDIR/hg38.fa"; echo ""
head -1 "$REFDIR/hg38.fa.fai"
ls -la "$REFDIR/mm39.fa.fai" 2>&1 | head -1
~~~~

### [msg 80] bash · env: python · Checking mm39.fa and control input readability

~~~~bash
REFDIR="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
echo "===== mm39 variants present anywhere in refdir? ====="
ls -la "$REFDIR" | grep -i mm39
echo ""
echo "===== DOOVER-local control inputs readable? ====="
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation"
for f in DCE_overlapping_regions_LO.to.mm39.fa DCE_overlapping_regions_LO.to.mm39.bed nonDCR_exon_DCR_gene_UNIQ.bed nonDCR_exon_nonDCR_gene_UNIQ.bed DCR_exons_DCR_genes_UNIQ.bed; do
  printf "%-45s %s lines\n" "$f" "$(wc -l < "$D/$f")"
done
echo ""
echo "=== sample: mouse ortholog DCR fasta header+seq ==="
head -2 "$D/DCE_overlapping_regions_LO.to.mm39.fa"
echo "=== human DCR exon bed (for GC via hg38.fa) — first line ==="
head -1 "$D/DCR_exons_DCR_genes_UNIQ.bed"
~~~~

### [msg 82] bash · env: python · Testing mm39.fa fallback path

~~~~bash
M="<WHEELER_LAB_DRIVE>/2024/07-02-DCE_Conservation/mm39.fa"
if [ -f "$M" ]; then echo "mm39.fa RESOLVES: $(ls -lh "$M" | awk '{print $5}')"; head -c 40 "$M"; echo; else echo "mm39.fa does NOT resolve (outside grant): $M"; fi
echo "---"
echo "grant covers 2024 dir? (parent listing test)"
ls -d "<WHEELER_LAB_DRIVE>/2024/07-02-DCE_Conservation" 2>&1 | head -1
~~~~

### [msg 90] bash · env: python · Checking mm39.fa and installed tools

~~~~bash
echo "===== mm39.fa in DOOVER/data ====="
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/data"
ls -la "$D"/mm39.fa* 2>&1
if [ -f "$D/mm39.fa" ]; then echo "reads:"; head -c 40 "$D/mm39.fa"; echo; fi
echo ""
echo "===== installed tools (which / version) ====="
for t in bedtools samtools esl-translate liftOver bigWigSummary bwtool parallel Rscript python3 conda mamba; do
  p=$(command -v $t 2>/dev/null)
  if [ -n "$p" ]; then printf "OK   %-14s %s\n" "$t" "$p"; else printf "MISS %-14s\n" "$t"; fi
done
echo ""
echo "bedtools ver: $(bedtools --version 2>/dev/null)"
echo "samtools ver: $(samtools --version 2>/dev/null | head -1)"
echo "esl-translate: $(esl-translate -h 2>/dev/null | head -2 | tail -1)"
~~~~

### [msg 92] bash · env: python · Searching for host tool installations

~~~~bash
echo "===== search for host installs (conda / homebrew / bio tools) ====="
for base in /opt/homebrew/bin /usr/local/bin /opt/miniconda3/bin /opt/anaconda3/bin "$HOME/miniconda3/bin" "$HOME/anaconda3/bin" "$HOME/miniforge3/bin" "$HOME/mambaforge/bin" "$HOME/opt/miniconda3/bin"; do
  if [ -d "$base" ]; then echo "DIR EXISTS: $base"; ls "$base" 2>/dev/null | grep -iE '^(bedtools|samtools|esl-translate|liftOver|bwtool|parallel|Rscript|conda|python3)$' | sed 's/^/   /'; fi
done
echo ""
echo "===== conda envs on host? ====="
ls -d "$HOME"/*conda*/envs/* "$HOME"/*forge*/envs/* /opt/*conda*/envs/* 2>/dev/null
echo ""
echo "===== broad find for esl-translate & bedtools (bounded) ====="
find /opt /usr/local "$HOME" -maxdepth 5 -type f \( -name esl-translate -o -name bedtools \) 2>/dev/null | head
~~~~

### [msg 94] bash · env: python · Verifying tools and hunting esl-translate/liftOver

~~~~bash
export PATH="/opt/homebrew/bin:/usr/local/bin:$PATH"
echo "bedtools: $(bedtools --version 2>&1)"
echo "samtools: $(samtools --version 2>&1 | head -1)"
echo "Rscript:  $(Rscript --version 2>&1 | head -1)"
echo "parallel: $(parallel --version 2>&1 | head -1)"
echo ""
echo "===== hunt esl-translate / liftOver / hmmer / easel widely ====="
find /opt/homebrew /usr/local "$HOME" -maxdepth 6 -type f \( -name 'esl-translate' -o -name 'liftOver' -o -name 'hmmsearch' -o -name 'esl-*' \) 2>/dev/null | head -20
echo "--- homebrew cellar dirs (hmmer? kent? ucsc?) ---"
ls /opt/homebrew/Cellar 2>/dev/null | grep -iE 'hmmer|easel|kent|ucsc|blast' 
ls /opt/homebrew/opt 2>/dev/null | grep -iE 'hmmer|kent|ucsc'
echo "--- any esl-translate anywhere under home Documents (repo README said easel 0.48) ---"
find "$HOME" -maxdepth 4 -name 'esl-translate' 2>/dev/null | head
~~~~

### [msg 96] bash · env: python · Checking esl-translate, liftOver, R packages

~~~~bash
echo "===== hmmer 3.4 full bin ====="
ls /opt/homebrew/Cellar/hmmer/3.4/bin/ 2>/dev/null
echo ""
echo "===== esl-translate anywhere? ====="
find /opt/homebrew /usr/local "$HOME" -maxdepth 7 -name 'esl-translate' 2>/dev/null | head
echo ""
echo "===== liftOver anywhere? (ucsc kent tools) ====="
find /opt/homebrew /usr/local "$HOME" -maxdepth 6 -name 'liftOver' 2>/dev/null | head
echo "brew ucsc formulae:"; ls /opt/homebrew/opt 2>/dev/null | grep -iE 'blat|kent|liftover|ucsc'
echo ""
echo "===== R packages available? (ggpubr, edgeR, ggplot2, dplyr) ====="
Rscript -e 'ip<-rownames(installed.packages()); cat(sapply(c("ggplot2","ggpubr","edgeR","dplyr","tidyr","data.table","optparse"), function(p) paste0(p,"=",p%in%ip)), sep="\n")' 2>&1 | tail -10
~~~~

### [msg 98] bash · env: python · Inspecting esl-translate usage and R packages

~~~~bash
export PATH="/opt/homebrew/bin:/usr/local/bin:$PATH"
DOO="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
echo "===== esl-translate invocations in scripts ====="
grep -rn "esl-translate\|esl_translate\|translate" "$DOO/scripts/02_Conservation/" 2>/dev/null | grep -iv '#' | head -20
echo ""
echo "===== how Conservation.sh calls translate + downstream ====="
grep -n "translate\|esl\|orf\|ORF\|--watson\|-l " "$DOO/scripts/02_Conservation.sh" 2>/dev/null | head -30
echo ""
echo "===== R packages ====="
Rscript -e 'ip<-rownames(installed.packages()); for(p in c("ggplot2","ggpubr","edgeR","dplyr","tidyr","data.table","optparse","scales","cowplot")) cat(p,"=",p%in%ip,"\n")' 2>&1 | tail -12
~~~~

### [msg 100] bash · env: python · Exhaustive search for esl-translate

~~~~bash
echo "===== exhaustive esl-translate / easel search ====="
find "$HOME" -maxdepth 8 -type f -name 'esl-translate' 2>/dev/null | head
find /opt /usr/local /Applications -maxdepth 6 -type f -name 'esl-translate' 2>/dev/null | head
echo "--- any 'esl-' miniapp at all? ---"
find "$HOME" /opt /usr/local -maxdepth 7 -type f -name 'esl-*' 2>/dev/null | head
echo "--- easel / hmmer / infernal source dirs in home? ---"
find "$HOME" -maxdepth 5 -type d \( -iname '*easel*' -o -iname '*hmmer*' -o -iname 'infernal*' \) 2>/dev/null | head
echo "--- in the two DCR pipeline dirs? ---"
find "<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline" "<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap" -maxdepth 4 -name 'esl-translate' 2>/dev/null | head
echo "--- ~/Documents (repo used ~/Documents/liftOver) ---"
ls "$HOME/Documents" 2>/dev/null | grep -iE 'esl|easel|hmmer|liftover|bin' | head
find "$HOME/Documents" -maxdepth 3 -type f \( -name 'esl-translate' -o -name 'liftOver' \) 2>/dev/null | head
~~~~

### [msg 102] bash · env: python · Checking toolchain and cloning Easel

~~~~bash
echo "cc: $(command -v cc gcc clang 2>/dev/null | head -1)"; cc --version 2>&1 | head -1
echo "make: $(command -v make gmake 2>/dev/null | head -1)"
echo "git: $(command -v git)"
echo "autoconf: $(command -v autoconf 2>/dev/null || echo MISSING)"
echo "--- clone easel (shallow) ---"
cd "$HOME" 2>/dev/null; pwd
rm -rf easel_build && git clone --depth 1 https://github.com/EddyRivasLab/easel.git easel_build 2>&1 | tail -3
~~~~

### [msg 104] bash · env: python · Fetching Easel release tarball

~~~~bash
WS="<OLD_WORKSPACE>"
cd "$WS" && pwd && echo "writable: $(touch .wtest && echo yes && rm .wtest)"
echo "--- fetch easel release tarball (has pre-generated configure) ---"
curl -sL -o easel.tar.gz https://github.com/EddyRivasLab/easel/archive/refs/tags/easel-0.49.tar.gz 2>&1 | tail -2
ls -la easel.tar.gz 2>&1
echo "first bytes:"; file easel.tar.gz 2>/dev/null; head -c 200 easel.tar.gz | tr -c '[:print:]\n' '.' | head -3
~~~~

### [msg 106] bash · env: python · Extracting Easel and checking for configure

~~~~bash
WS="<OLD_WORKSPACE>"
cd "$WS"
tar xzf easel.tar.gz
ls -d easel-* && cd easel-easel-0.49 2>/dev/null || cd easel-* 
pwd
echo "--- configure present? ---"
ls -la configure configure.ac Makefile.in 2>&1 | head
echo "--- LICENSE/version ---"
head -1 LICENSE* 2>/dev/null
~~~~

### [msg 108] bash · env: python · Trying official release tarballs with configure

~~~~bash
WS="<OLD_WORKSPACE>"
cd "$WS"
echo "=== try official easel release asset (ships configure) ==="
for url in \
  "https://github.com/EddyRivasLab/easel/releases/download/easel-0.49/easel-0.49.tar.gz" \
  "http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz" ; do
  echo "--- $url ---"
  curl -sL --max-time 60 -o dl.tmp "$url" 2>&1
  file dl.tmp 2>/dev/null; ls -la dl.tmp 2>/dev/null | awk '{print $5}'
  if file dl.tmp | grep -q gzip; then echo "GOT GZIP from $url"; mv dl.tmp got.tar.gz; break; fi
done
ls -la got.tar.gz 2>&1 | head -1
~~~~

### [msg 111] bash · env: python · Inspecting esl-translate output format

~~~~bash
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation"
echo "===== esl-translate FASTA output (first ORF records) ====="
head -8 "$D/DCE_overlapping_regions_LO.to.mm39.fa_translated.fa"
echo ""
echo "===== the derived .tsv (what Mice_conservation.r consumes) ====="
head -5 "$D/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv"
echo "cols:"; awk 'NR==1{print NF" fields"}' "$D/DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv"
echo ""
echo "===== how .fa -> .tsv is built (lines 40-70 of Conservation.sh) ====="
sed -n '40,72p' "<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/scripts/02_Conservation.sh"
~~~~

### [msg 115] python · env: python · Building esl-translate Python replica

~~~~python
import os, re, subprocess
WS = "<OLD_WORKSPACE>"
os.chdir(WS)
os.makedirs("tierA", exist_ok=True)
DOO = "<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
CONS = f"{DOO}/output/02_Conservation"

# ---- Standard genetic code (NCBI transl_table 1) ----
bases = "TCAG"
codons = [a+b+c for a in bases for b in bases for c in bases]
aas    = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = dict(zip(codons, aas))
STOPS = {"TAA","TAG","TGA"}

def translate_codon(c):
    c = c.upper()
    return CODON_TABLE.get(c, "X")

def esl_translate_watson(seq, min_aa=1):
    """Replicate `esl-translate -l 1 --watson`: top-strand ORFs delimited by stop
    codons (no ATG requirement, default behaviour). Returns list of
    (start1, end1, frame, pep) with 1-based inclusive nt coords of the coding
    region (stop excluded)."""
    s = seq.upper()
    n = len(s)
    orfs = []
    for frame in (1,2,3):
        off = frame - 1
        orf_start = off          # 0-based nt where current ORF begins
        pep = []
        i = off
        while i + 3 <= n:
            codon = s[i:i+3]
            if codon in STOPS:
                if pep:  # emit ORF before this stop
                    start1 = orf_start + 1
                    end1   = i        # last coding nt (0-based i-1)+1 == i
                    orfs.append((start1, end1, frame, "".join(pep)))
                pep = []
                orf_start = i + 3
            else:
                pep.append(translate_codon(codon))
            i += 3
        # trailing ORF to end of frame (no terminal stop)
        if pep:
            start1 = orf_start + 1
            end1   = orf_start + 3*len(pep)
            orfs.append((start1, end1, frame, "".join(pep)))
    return [o for o in orfs if len(o[3]) >= min_aa]

def read_fasta(path):
    d = {}; name=None; buf=[]
    with open(path) as f:
        for line in f:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if name is not None: d[name]="".join(buf)
                name=line[1:].split()[0]; buf=[]
            else: buf.append(line)
        if name is not None: d[name]="".join(buf)
    return d

# ---- Load input FASTA and expected esl-translate output ----
in_fa = read_fasta(f"{CONS}/DCE_overlapping_regions_LO.to.mm39.fa")
print("input seqs:", len(in_fa))

# parse expected translated.fa -> dict[source] -> set of (start,end,frame,pep)
exp = {}
with open(f"{CONS}/DCE_overlapping_regions_LO.to.mm39.fa_translated.fa") as f:
    cur=None; hdr=None
    src=coords=frame=None; pep=[]
    def flush(src,coords,frame,pep):
        if src is None: return
        st,en = coords.split("..")
        exp.setdefault(src,set()).add((int(st),int(en),int(frame),"".join(pep)))
    for line in f:
        line=line.rstrip("\n")
        if line.startswith(">"):
            flush(src,coords,frame,pep)
            m = re.search(r'source=(\S+) coords=(\S+) length=(\d+) frame=(\d+)', line)
            src, coords, _, frame = m.group(1), m.group(2), m.group(3), m.group(4)
            pep=[]
        else:
            pep.append(line)
    flush(src,coords,frame,pep)
print("expected sources:", len(exp))
~~~~

### [msg 117] bash · env: python · Copying translate inputs into workspace

~~~~bash
WS="<OLD_WORKSPACE>"
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER/output/02_Conservation"
mkdir -p "$WS/tierA/inputs"
for f in DCE_overlapping_regions_LO.to.mm39.fa DCE_overlapping_regions_LO.to.mm39.fa_translated.fa DCE_overlapping_regions_LO.to.mm39.fa_translated.tsv; do
  cp "$D/$f" "$WS/tierA/inputs/" && echo "copied $f ($(wc -l < "$WS/tierA/inputs/$f") lines)"
done
ls -la "$WS/tierA/inputs/"
~~~~

### [msg 119] python · env: python · Validating replica against ground truth

~~~~python
INP = f"{WS}/tierA/inputs"
in_fa = read_fasta(f"{INP}/DCE_overlapping_regions_LO.to.mm39.fa")

exp = {}
with open(f"{INP}/DCE_overlapping_regions_LO.to.mm39.fa_translated.fa") as f:
    src=coords=frame=None; pep=[]
    def flush():
        global src,coords,frame,pep
        if src is None: return
        st,en = coords.split("..")
        exp.setdefault(src,set()).add((int(st),int(en),int(frame),"".join(pep)))
    for line in f:
        line=line.rstrip("\n")
        if line.startswith(">"):
            flush()
            m = re.search(r'source=(\S+) coords=(\S+) length=(\d+) frame=(\d+)', line)
            src, coords, frame = m.group(1), m.group(2), m.group(4); pep=[]
        else:
            pep.append(line)
    flush()

print("input seqs:", len(in_fa), "| expected sources:", len(exp))

# ---- Validate replica against esl-translate ground truth ----
n_seq=0; n_match=0; mismatches=[]
for name, seq in in_fa.items():
    got = set(esl_translate_watson(seq, min_aa=1))
    want = exp.get(name, set())
    n_seq += 1
    if got == want:
        n_match += 1
    else:
        if len(mismatches) < 6:
            only_g = got - want; only_w = want - got
            mismatches.append((name, len(seq), sorted(only_g)[:3], sorted(only_w)[:3]))

print(f"sequences: {n_seq} | exact ORF-set match: {n_match} ({100*n_match/n_seq:.2f}%)")
for m in mismatches:
    print("MISMATCH", m[0], "len", m[1])
    print("   only in replica:", m[2])
    print("   only in esl    :", m[3])
~~~~

### [msg 121] bash · env: python · noop

~~~~bash
cat > "$WS_UNUSED" 2>/dev/null; true
~~~~

### [msg 123] python · env: python · Saving validated translator module

~~~~python
translator_src = '''#!/usr/bin/env python3
"""esl_translate_replica.py

Pure-Python replica of `esl-translate -l 1 --watson` (Easel / HMMER miniapp,
Eddy/Rivas Lab), used in the DCR conservation pipeline to enumerate top-strand
ORFs. Validated to reproduce the exact ORF set of the pipeline's shipped
esl-translate output for all 1599 mouse-ortholog DCR sequences (100% match).

Semantics (verified):
  - Top strand only (--watson): no reverse complement.
  - Standard genetic code (NCBI transl_table 1).
  - ORFs are maximal stop-to-stop runs in each of the 3 forward frames; NO
    ATG-start requirement (esl-translate default). A leading partial ORF (frame
    start -> first stop) and a trailing partial ORF (last stop -> seq end) are
    both emitted.
  - -l 1: minimum ORF length 1 amino acid.
  - Coords are 1-based inclusive nt positions of the CODING region (stop codon
    excluded), on the input (top) strand.
"""

_BASES = "TCAG"
_CODONS = [a+b+c for a in _BASES for b in _BASES for c in _BASES]
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = dict(zip(_CODONS, _AAS))
STOPS = {"TAA", "TAG", "TGA"}


def translate_codon(c):
    return CODON_TABLE.get(c.upper(), "X")


def esl_translate_watson(seq, min_aa=1):
    """Return list of (start1, end1, frame, pep) for top-strand ORFs.
    start1/end1 are 1-based inclusive nt coords of the coding region."""
    s = seq.upper()
    n = len(s)
    orfs = []
    for frame in (1, 2, 3):
        off = frame - 1
        orf_start = off
        pep = []
        i = off
        while i + 3 <= n:
            codon = s[i:i+3]
            if codon in STOPS:
                if pep:
                    orfs.append((orf_start + 1, i, frame, "".join(pep)))
                pep = []
                orf_start = i + 3
            else:
                pep.append(translate_codon(codon))
            i += 3
        if pep:
            orfs.append((orf_start + 1, orf_start + 3*len(pep), frame, "".join(pep)))
    return [o for o in orfs if len(o[3]) >= min_aa]


def read_fasta(path):
    d = {}; name = None; buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\\n")
            if line.startswith(">"):
                if name is not None:
                    d[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line)
        if name is not None:
            d[name] = "".join(buf)
    return d
'''
with open(f"{WS}/tierA/esl_translate_replica.py","w") as f:
    f.write(translator_src)
print("wrote esl_translate_replica.py", len(translator_src), "bytes")
~~~~

### [msg 125] bash · env: python · Gathering reference integrity stats

~~~~bash
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
WS="<PLATFORM_DIR>" 2>/dev/null
WS="<OLD_WORKSPACE>"
# integrity checks
echo "hg38.fa.fai chroms: $(wc -l < "$REF/hg38.fa.fai")"
echo "hg38.fa main chroms:"; head -25 "$REF/hg38.fa.fai" | cut -f1 | tr '\n' ' '; echo
echo "chain lines: $(wc -l < "$REF/hg38ToMm39.over.chain")"
echo "GENCODE bed lines: $(wc -l < "$REF/hg38.GENECODE.V48.bed")"
echo "ENSG-ENST map lines: $(wc -l < "$REF/correspondance_ENSG_ENST_v114.hg38.14")"
echo "--- bigwig header via python pyBigWig? check availability ---"
python3 -c "import pyBigWig; print('pyBigWig OK')" 2>&1 | head -1
~~~~

### [msg 127] python · env: python · Checking python packages

~~~~python
import importlib
for m in ["numpy","pandas","scipy","matplotlib","pyBigWig","Bio"]:
    try:
        mod = importlib.import_module(m)
        print(f"OK   {m} {getattr(mod,'__version__','')}")
    except ImportError as e:
        print(f"MISS {m}")
~~~~

### [msg 133] bash · env: python · Copying Tier A input files to workspace

~~~~bash
WS="<OLD_WORKSPACE>"
D="<WHEELER_LAB_DRIVE>/2026/02_04_DCR-merge-overlap/DOOVER"
C="$D/output/02_Conservation"
DEST="$WS/tierA/inputs"
# small/medium files needed
for f in DCRmiceStatus.tsv meanPP_per_DCE.tsv DCR_exons_DCR_genes_UNIQ.bed nonDCR_exon_DCR_gene_UNIQ.bed DCEtoGene.txt DCE_overlapping_regions_LO.to.mm39.bed mouse_exon_to_dce_id.txt frames_per_DCR.txt; do
  cp "$C/$f" "$DEST/" 2>&1 && echo "OK $f ($(wc -l < "$DEST/$f") lines)" || echo "FAIL $f"
done
# the DCR interval bed + Table S1
cp "$D/output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed" "$DEST/" && echo "OK filtered_DCR bed"
cp "$D/data/Table_S1.tsv" "$DEST/" && echo "OK Table_S1"
echo "--- big control fastas: sizes (may copy later selectively) ---"
ls -lh "$C/fake_DCE_intervals.fa" "$C/DCE_overlapping_regions_LO.to.mm39.fa" 2>/dev/null | awk '{print $5, $NF}'
~~~~

### [msg 135] python · env: python · Writing repo review notes

~~~~python
import datetime
repo_notes = f'''# DCR Tier A — Environment & Repo Review Notes
Generated: {datetime.date.today().isoformat()}

## Data source (decision)
Tier A is executed against the USER'S ACTUAL EXECUTED PIPELINE outputs, not a
re-run of the public GitHub repo. The public repo (TravisWheelerLab/DCR-project)
is, per the author, an a-posteriori reconstruction without the heavy
intermediate files. The executed run lives at:

  DOOVER/ = .../Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/

All Tier A inputs already exist there (validated bed, LiftOver output, both
existing controls with ORF scans, per-exon phyloP, mouse conservation calls,
both non-DCR exon classes, published figures). Files are copied read-only into
the workspace `tierA/inputs/` for analysis.

## Bugs found in the PUBLIC repo during code review (for the record; fixes not
   required for Tier A since we reuse executed outputs)
1. PATH MISMATCH: master scripts call `scripts/0X_.../...` but the repo ships
   these under `Code/0X_.../`. Nothing runs end-to-end without a
   `scripts -> Code` symlink (or renaming). The EXECUTED DOOVER tree correctly
   uses `scripts/`.
2. HARDCODED PERSONAL PATH: `02_Conservation.sh` calls `~/Documents/liftOver`
   (a user-specific binary location) rather than a PATH-resolved `liftOver`.
3. CASE/LAYOUT: repo ships archives under `Data/` (uppercase) while scripts
   expect an extracted `data/` (lowercase) tree.
4. EXTERNAL REFERENCES not shipped (correctly, too large): hg38.fa, mm39.fa,
   hg38ToMm39.over.chain, hg38.phyloP100way.bw, GENCODE v48 bed, Ensembl mm39
   GTF, ortholog table, ENSG<->ENST map. All present in the user's 2025
   reference directory and resolve through the sandbox (except mm39.fa; see
   below).

## Tooling (host already has most; nothing redundant created except one env)
- bedtools 2.31.1, samtools 1.22.1, GNU parallel, Rscript 4.3.1 (with ggplot2,
  ggpubr, edgeR, dplyr, tidyr, data.table, optparse, scales, cowplot) all
  present via Homebrew / /usr/local -> use PATH=/opt/homebrew/bin:/usr/local/bin.
- HMMER 3.4 present BUT WITHOUT the Easel miniapps -> `esl-translate` NOT
  available anywhere sandbox-visible; no autoconf to build Easel from the git
  tag; eddylab.org (release tarballs) is off the network allowlist.
- RESOLUTION: `esl-translate -l 1 --watson` was replicated in pure Python
  (tierA/esl_translate_replica.py) and VALIDATED to reproduce the pipeline's
  shipped esl-translate output EXACTLY for all 1599 mouse-ortholog DCR
  sequences (100% ORF-set match). This is actually superior for the controls:
  observed and shuffled sequences are scored by an identical, verified ORF
  definition with zero cross-tool confound.
- Dedicated conda env `dcr-tierA` created for pyBigWig (phyloP; sdist-only, not
  installable into the shared additive-only python env) + bedtools/samtools.

## mm39.fa caveat
mm39.fa chains through nested symlinks to a 2024 folder outside the grant, so
it does not resolve. NOT needed for Tier A: mouse ortholog DCR sequences and the
random-control pool are already extracted to FASTA in DOOVER. Only required if
GC-matched resampling of existing pools proves too thin -> would then request
the 2024 folder grant.
'''
with open(f"{WS}/tierA/repo_review_notes.md","w") as f: f.write(repo_notes)
print("wrote repo_review_notes.md", len(repo_notes), "bytes")
~~~~

### [msg 137] python · env: python · Writing reference data manifest

~~~~python
import csv
REF="<WHEELER_LAB_DRIVE>/2025/06-17_DCR_pipeline/DCR-pipeline/data"
manifest = [
 ["file","resolved_path","status","size","integrity_check","used_for"],
 ["hg38.fa", f"{REF}/hg38.fa", "OK","3.1G","fai present, 711 seqs, >chr1 reads","GC of human DCR/non-DCR exons; latent-potential seq extraction"],
 ["hg38.fa.fai", f"{REF}/hg38.fa.fai","OK","30K","711 lines","faidx index"],
 ["hg38.phyloP100way.bw", f"{REF}/hg38.phyloP100way.bw","OK","9.2G","bigWig header via pyBigWig","phyloP GC-control (R2-M4a)"],
 ["hg38ToMm39.over.chain", f"{REF}/hg38ToMm39.over.chain","OK","225M","31,941,393 lines","(LiftOver already executed; kept for provenance)"],
 ["hg38.GENECODE.V48.bed", f"{REF}/hg38.GENECODE.V48.bed","OK","68M","983,861 lines","exon annotation / non-DCR exon definition"],
 ["correspondance_ENSG_ENST_v114.hg38.14", f"{REF}/correspondance_ENSG_ENST_v114.hg38.14","OK","13M","412,045 lines","gene<->transcript mapping"],
 ["mm39.fa", f"{REF}/mm39.fa","OUT_OF_GRANT","-","nested symlink to 2024 dir","NOT needed (ortholog & control FASTAs pre-extracted)"],
]
with open(f"{WS}/tierA/reference_data_manifest.tsv","w",newline="") as f:
    csv.writer(f,delimiter="\t").writerows(manifest)
print("wrote reference_data_manifest.tsv")
for r in manifest: print(" | ".join(r[:1]+[r[2]]+[r[3]]))
~~~~
