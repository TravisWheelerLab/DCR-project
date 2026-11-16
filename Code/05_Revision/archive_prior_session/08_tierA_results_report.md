# 08_tierA_results_report — as executed, not refactored

**Assembly of the Tier A results document (no new analysis).**

| | |
|---|---|
| Source | exported revision session `31323eb8_2026-09-16T22-37-01.json`, `frames[0].messages[332..342]` (assistant `tool_use` blocks, `input.code` / `input.command`) |
| Tool calls | 4 (msg 332 … msg 338) |
| Language / environment | python (env: dcr-tierA); repl |
| Saved outputs | `TierA_results.md` |

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

### [msg 332] python · env: dcr-tierA · Preparing artifact ids for report

~~~~python
import json
# collect artifact ids for embedding
art={
 'gc_audit':'27592a7f-a6e2-4362-a075-fe737b76d984',
 'gc_matched':'6981ae4f-98c2-4fd6-988d-c0edca298d21',
 'self_shuffle':'d980e8ca-96c6-4273-8f06-ad64b98b2c3a',
 'latent':'497162ff-63a0-4ee7-a8ae-b20de0d3a742',
 'phylop':'7e474197-09da-4d82-b5a2-c47406aacb82',
 'gnas':'c1ac97ab-857f-46c8-a111-c06cecaf29a5',
 'baseline_2B':'2ca2af3a-4440-4e21-a9f9-86b5a460b21f',
 'baseline_2C':'ab564c52-4440-4e21-a9f9-86b5a460b21f'.replace('4440-4e21-a9f9-86b5a460b21f','92d9-...'),  # placeholder guard
}
# fix baseline_2C from summary
art['baseline_2C']='ab564c52-4440'  # will not embed if wrong; use known good below
print("ok")
~~~~

### [msg 334] repl · env: control-plane (stdlib) · Resolving baseline artifact ids

~~~~python
r=host.artifacts(filename="baseline_fig2", limit=20)
for a in r['artifacts']:
    print(a['filename'], a['artifact_id'])
~~~~

### [msg 336] repl · env: control-plane (stdlib) · Resolving baseline artifact ids

~~~~python
r=host.artifacts(filename="baseline_fig2", limit=20)
import json
for a in r['artifacts']:
    print(a.get('filename'), a.get('id'), a.get('latest_version_id'))
~~~~

### [msg 338] python · env: dcr-tierA · Writing Tier A results document

~~~~python
A={
 'baseline_2B':'2ca2af3a-4440-4e21-a9f9-86b5a460b21f',
 'baseline_2C':'ab564c52-47c9-45c6-81f2-2fcbcff7b1ef',
 'gc_audit':'27592a7f-a6e2-4362-a075-fe737b76d984',
 'gc_matched':'6981ae4f-98c2-4fd6-988d-c0edca298d21',
 'self_shuffle':'d980e8ca-96c6-4273-8f06-ad64b98b2c3a',
 'latent':'497162ff-63a0-4ee7-a8ae-b20de0d3a742',
 'phylop':'7e474197-09da-4d82-b5a2-c47406aacb82',
 'gnas':'c1ac97ab-857f-46c8-a111-c06cecaf29a5',
}
def img(k): return "{{artifact:art_%s}}"%A[k]

md=f"""# DCR revision — Tier A results package

**Paper:** Goubert et al., *Alternatively Spliced Dual-Coding Regions Contribute to the Human Gene Regulatory Program* (bioRxiv v2, DOI 10.64898/2026.02.24.707805)
**Scope:** Tier A = confound controls and clarifications answerable by reanalysis of the already-executed pipeline (DOOVER outputs), no new wet-lab or transcriptome work. All analyses re-derive dual-coding calls with the pipeline's exact rule and a pure-Python replica of `esl-translate -l 1 --watson` (validated to 100% ORF-set identity on all 1599 mouse-ortholog sequences).

Reviewer items addressed here: **R2-M1** (GC confound + missing composition null), **R2-M4a** (phyloP GC control), **R1-M2a** (necessary-vs-sufficient), **R1-m1** (GNAS "data not shown").

---

## 0. Baseline reproduction (ground-truth check)

Before modifying anything, the published Figure 2 results were reproduced from the executed pipeline outputs to confirm we operate on the same data.

- **DCR interval set:** 1741 intervals / 1407 unique DCRs — exact match to the paper.
- **Mouse conservation status:** DCR = 981, not-DCR = 250, not-mapped = 112, too-short = 63. The 981 conserved (≈80% of mappable) matches the paper's headline figure.
- **Fig 2B (dual-coding proportions):** DCR mouse-ortholog **82.6%** > random-region control **60.5%** > codon-shuffle control **32.3%**.
- **Fig 2C (phyloP):** DCR exons are *less* nucleotide-conserved than controls — median phyloP D = 2.61 < D\\* = 2.95 < ND = 3.05; all pairwise comparisons significant, same direction as the paper.

![Baseline Fig 2B]({img('baseline_2B')})

![Baseline Fig 2C]({img('baseline_2C')})

The dual-coding rule (from `Mice_conservation.r`) counts full-length ORFs — those spanning the entire region within ±3 bp at each end — and calls a region dual-coding when it has **>1** such ORF. Every control below applies this identical rule.

---

## 1. R2-M1 — GC audit: is dual-coding an artifact of GC content?

Referee 2's central concern is that DCR exons might be GC-enriched, and that GC-rich sequence trivially yields more ORFs (fewer stop codons). We first quantified the GC difference.

**DCR exons ARE GC-enriched, but modestly.**

| Class | n | median GC | vs DCR |
|---|---|---|---|
| Human DCR exon | 3,324 | **54.2%** | — |
| non-DCR exon, DCR gene | 21,031 | 51.5% | −2.8 pts |
| non-DCR exon, non-DCR gene | 286,492 | 51.2% | −3.0 pts |
| Mouse DCR ortholog | 1,599 | 53.5% | — |
| Random-region control (mouse) | 100,000 | 52.3% | −1.2 pts |

The enrichment is highly significant given the sample sizes (MWU p < 1e-24 human, p = 3e-5 mouse) but the **effect size is small**: rank-biserial ≈ 0.11 (human) and ≈ 0.06 (mouse). So the confound is real and must be controlled — but it is a few percentage points, not a dominant driver.

![GC audit]({img('gc_audit')})

---

## 2. R2-M1 — GC-matched Fig 2B: the excess survives GC control

We re-ran the Fig 2B comparison after matching the random-region control pool to the DCR orthologs on **joint length × GC** (two independent methods: 2-D bin matching and nearest-neighbour matching in standardized length–GC space, with bootstrap CI).

- DCR observed: **82.6%** (unchanged).
- Length-only control: 60.5%.
- **GC + length-matched control: 64.8%** (95% CI 62.9–66.8; 2-D bin method agrees at 65.4%).

Matching quality: matched-control median GC 53.6 = DCR 53.6; median length 64 ≈ DCR 64.

**GC explains ~5 of the 22-point gap. A ~18-point DCR excess survives GC + length control.**

![GC-matched Fig 2B]({img('gc_matched')})

---

## 3. R2-M1 — the composition null (DCR self-shuffle)

Referee 2 specifically noted the paper's codon-shuffle control was built from **non-DCR** mouse CDS, so it does not hold the DCR exons' own composition constant. We supplied the missing control: apply the pipeline's exact shuffle procedure (frame-0 codon split, drop trailing partial codon, remove stop codons, shuffle) to **each DCR ortholog's own sequence**. This holds nucleotide composition — and therefore GC — *exactly* constant by construction, randomizing only codon arrangement.

- Observed DCR orthologs dual-coding: **82.4%** (from-sequence recomputation; validates the rule against the 82.6% TSV value).
- Self-shuffled (own codons, N = 20 shuffles/exon): **68.6%**.
- **Arrangement effect: 13.8 points.**

Empirical null over 1,000 whole-set shuffles: observed **1,318** dual-coding orthologs vs null mean **1,096** (95% interval 1,067–1,124), **Z = 15.5, p < 0.001**.

This is the strongest single answer to R2-M1: with GC and codon usage identical, the actual codon *arrangement* of DCR orthologs produces far more dual-coding than chance.

![Self-shuffle composition null]({img('self_shuffle')})

---

## 4. R1-M2a — necessary vs sufficient: latent potential in ordinary exons

Referee 1 asked whether dual-coding potential is *sufficient*, or merely *necessary*, for a region to be used as a DCR. We scored the same dual-coding potential (esl rule) on human exons and compared DCR exons against GC + length-matched ordinary (non-DCR) exons.

| Class | latent dual-coding potential |
|---|---|
| DCR exon (observed) | **48.3%** (CI 46.7–49.9) |
| non-DCR exon, DCR gene (GC+len matched) | 40.9% (CI 39.2–42.4) |
| non-DCR exon, non-DCR gene (GC+len matched) | 41.4% (CI 39.7–43.0) |

Dual-coding potential is **enriched** in DCRs (48% vs 41%) but is **common** in ordinary matched exons — roughly two in five carry latent potential they never realize. **Potential is necessary but not sufficient.** This supports softening the causal language around pillar 1 (conserved potential) and shifting argumentative weight to pillar 2 (tissue-specific *realized* frame usage), which distinguishes true DCRs.

![Latent potential]({img('latent')})

---

## 5. R2-M4a — phyloP conservation with GC control

Referee 2 asked whether the reduced nucleotide conservation of DCR exons is a GC artifact (GC-rich sequence tends to score lower on phyloP). Two complementary controls:

**(a) GC + length-matched resampling** (controls matched to DCR exons on length and GC):

| | mean phyloP | median |
|---|---|---|
| DCR exon | 2.60 | 2.61 |
| matched non-DCR, DCR gene | 2.73 | 2.82 |
| matched non-DCR, non-DCR gene | 2.81 | 3.00 |

DCR exons remain significantly *less* conserved (MWU p = 1.3e-3 and 2.8e-7).

**(b) Linear model** `meanPP ~ is_DCR + GC + length`: the DCR coefficient is **−0.194 phyloP units** (95% CI −0.253 to −0.134, **p = 2.1e-10**), i.e. DCR exons are less conserved even after adjusting for GC and length. GC itself has a strong negative association (−0.039 per %, p ≈ 0), confirming the confound exists — but the DCR effect is independent of it.

**Reduced DCR nucleotide conservation is not a GC artifact; the direction matches the paper.**

![phyloP with GC control]({img('phylop')})

---

## 6. R1-m1 — the GNAS "data not shown"

The paper states only 3 human genes have an exon overlapping a DCR > 500 nt, and that GNAS is the only one whose dual-coding potential is conserved in mouse ("data not shown"). We reproduced the supporting data:

| gene | DCR length (nt) | mouse status |
|---|---|---|
| GBGT1 | 523 | not DCR |
| ZNF586 | 594 | not mapped |
| **GNAS** | **1,878** | **DCR (conserved)** |

Exactly 3 genes, as claimed, and only GNAS is conserved as dual-coding. Direct ORF confirmation: the GNAS mouse ortholog region (2,154 nt) carries **two full-length ORFs** — frame 1 (718 aa) and frame 3 (717 aa) — both spanning the region. This figure/table can replace "data not shown".

![GNAS data]({img('gnas')})

---

## Summary for the rebuttal

| Reviewer item | Result | Verdict |
|---|---|---|
| R2-M1 GC confound | DCR exons +2.8–3.0 pts GC (small effect size) | Confound real but modest |
| R2-M1 GC-matched Fig 2B | 82.6% vs 64.8% matched control | ~18-pt excess survives |
| R2-M1 composition null | 82.4% → 68.6% self-shuffle; Z=15.5, p<0.001 | Arrangement effect, GC held exactly constant |
| R1-M2a necessary/sufficient | 48% DCR vs ~41% matched ordinary exons | Necessary, not sufficient → lean on pillar 2 |
| R2-M4a phyloP GC control | is_DCR coef −0.19, p=2e-10 (GC-adjusted) | Not a GC artifact |
| R1-m1 GNAS | 3 genes >500nt; only GNAS conserved (2 full ORFs) | Claim substantiated |

**Bottom line:** the GC confound Referee 2 identified is genuine but accounts for only a small fraction of the dual-coding signal. Every DCR result holds after GC control, and the self-composition null — the specific control R2 said was missing — shows the effect is driven by codon arrangement with GC held exactly constant. The one honest concession is R1-M2a: dual-coding *potential* is necessary but not sufficient, which motivates foregrounding the tissue-specific realized-usage evidence (pillar 2) as the distinguishing feature of true DCRs.
"""
with open(f"{WS}/tierA/TierA_results.md","w") as f: f.write(md)
print("wrote TierA_results.md,",len(md),"chars")
~~~~
