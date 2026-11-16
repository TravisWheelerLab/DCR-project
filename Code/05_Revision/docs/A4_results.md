# A4 — expression class and tissue specificity (R2-M5, R2-M6)

GTEx v8 (GENCODE v26, 54 tissue sites) median TPM via the pinned connector (in this repository: `controls/gtex_fetch_medians.py`, public GTEx Portal API v2; identical values).
Background = all genes with mapped SwissProt isoforms (17,472); **16,532
resolved (94.6%), including 1,230 of 1,296 DCR genes (94.9%)**, 0 failed
batches. Tissue specificity is Yanai τ = Σ(1 − x_i/x_max)/(n − 1).

Overall DCR rate in the resolved background: **7.44%**.

## R2-M5 — DCR incidence is independent of expression level

| expression quintile | median TPM | genes | DCR genes | % DCR [95% CI] |
|---|---|---|---|---|
| Q1 (lowest) | 0.05 | 3,307 | 203 | 6.14 [5.32, 6.96] |
| Q2 | 1.50 | 3,306 | 246 | 7.44 [6.55, 8.34] |
| Q3 | 7.17 | 3,306 | 271 | 8.20 [7.26, 9.13] |
| Q4 | 16.88 | 3,306 | 269 | 8.14 [7.21, 9.07] |
| Q5 (highest) | 46.57 | 3,307 | 241 | 7.29 [6.40, 8.17] |

Logistic regression on log₁₀(median TPM + 1): OR = 1.039, **p = 0.41**. Flat.

**This answers R2-M5 favourably.** The referee's link to Bénitière et al. implies
that if DCRs were a by-product of transcriptional/splicing output, their
incidence should track expression. It does not — the span across quintiles is
~2 percentage points with overlapping intervals and no trend. DCR incidence
cannot be explained as a simple function of how much a gene is expressed.

## R2-M6 — DCR genes are not more tissue-specific

This is the reference class the referee asked for, and the answer is negative.

| | τ median | n |
|---|---|---|
| DCR genes | 0.7765 | 1,229 |
| non-DCR genes | 0.7952 | 15,253 |

Mann–Whitney p = 0.0135, **rank-biserial 0.042** — negligible, significant only
at this n, and in the *opposite* direction to the hypothesis.

τ is dominated by expression level, so it must be stratified:
`OLS τ ~ is_DCR + log₁₀(TPM+1)` gives **is_DCR = −0.0087 (p = 0.018)** against
**logexp = −0.183 (p < 1e-300)**. Within every expression quintile the
difference is **not significant** (p = 0.26, 0.58, 0.76, 0.41, 0.35).

**Brain dominance (gene level):** 15.95% of DCR genes have their maximum median
TPM in a brain tissue, versus 17.99% of non-DCR genes — OR 0.865, p = 0.075.
**No gene-level brain enrichment.**

## How to write this up

1. **Concede R2-M6 at the gene level.** DCR-containing genes are not unusually
   tissue-specific. State it plainly; it costs nothing, because Fig 4's claim is
   a different object — *which frame dominates* within a tissue, not whether the
   gene is tissue-restricted. A gene can be broadly expressed while the frame
   ratio flips between tissues. Make that distinction explicit so the reference
   class does not read as a contradiction of Fig 4.
2. **Claim R2-M5.** Incidence independent of expression is a positive result
   against the noise interpretation, and it is the cheapest win in the revision.
3. **Check the brain language.** Any sentence implying DCR genes are
   brain-enriched as a class must be softened — the gene-level test says they are
   not. The paper's existing framing (a small number of ancient DCRs with
   brain-specific transcripts) is consistent with this result and should be the
   only form the claim takes. Tress et al.'s protein-level confirmation of
   ATP2B1, ATP2B4, MADD and MARK4 supports the specific-genes version, not a
   class-wide one.

   Two further details from their text. First, and this is the one that matters:
   they report that a large majority of their genes with **ancient** dual coding
   regions are highly expressed in brain (their Figure 4), with only **a small
   number of the 20** not enriched in nervous tissues. So at gene level their
   ancient set *is* brain-enriched.

   That is not in conflict with our null — the two tests cover different sets.
   Ours asks whether DCR genes **as a class** (1,230 genes) are brain-dominant,
   and they are not (15.95% vs 17.99%). Theirs asks about a small, *age-selected*
   subset. Taken together they support exactly the claim the paper already makes
   — that the brain signal comes from a subset of ancient DCRs rather than from
   DCR genes generally — which means the existing hedge is the correct framing
   and is now backed from two directions rather than needing to be softened
   further. **Do not state or imply the reverse** (that most of their 20 were
   unenriched); the sentence is split across a page break in the PDF and reads
   that way only if the running footer is not stripped.

   They also note that a few genes which are *not* brain-enriched nonetheless
   carry brain-specific protein isoforms (ADGRL2, ATP2B4, LMO7) or transcripts
   (CCSER2, MAGI1, PBX1) — a minority pattern, but the one that shows gene-level
   expression breadth and isoform-level tissue restriction are separable, which
   is what licenses conceding R2-M6 at gene level without touching Fig 4.

   Second — relevant to how we present the selection test — they state that in
   all but a handful of their regions one frame (the principal frame) is under
   much clearer purifying selection than the other, citing Szklarczyk et al.
   2007 as consistent. That is precisely why our estimator conditions on
   canonical-synonymous sites, and it is worth citing at that point in the
   Methods.

   Note for the response letter: do not merge their gene lists. The
   protein-level confirmation set is ATP2B1/ATP2B4/MADD/MARK4; the
   ancient-DCR-with-brain-specific-transcript set is larger and splits by age
   (ATP2B1, ATP2B4, PBX1, MARK4 pre-vertebrate; CCSER2, MADD, MAGI1, PLCB4,
   PTBP2 across jawed vertebrates; GABBR1 across euteleostomi). Eight of those
   ten were in their analysis — PTBP2 and GABBR1 were not.

## Follow-up: does conservation recover brain enrichment? No.

Because Tress's ancient set *is* brain-enriched at gene level, the obvious
question is whether our own conserved DCRs are too. Stratifying our DCR genes by
mouse status (`DCRmiceStatus.tsv`, parsed manually — multi-word statuses break
delimiter readers):

| group | n | brain-dominant | τ median | vs non-DCR |
|---|---|---|---|---|
| DCR, mouse-conserved | 884 | 15.84% | 0.771 | OR 0.858, p = 0.11 |
| DCR, not conserved | 208 | 16.83% | 0.775 | OR 0.922, p = 0.72 |
| DCR, unscored | 137 | 15.33% | 0.815 | — |
| non-DCR genes | 15,253 | **17.99%** | 0.795 | — |

**No DCR stratum is brain-enriched at gene level** — not even the mouse-conserved
one. Three strata, three nulls, all slightly below the non-DCR baseline.

This does not contradict Tress. Their "ancient" means pre-vertebrate to
jawed-vertebrate (~450 Mya) on a hand-curated set of 20; mouse conservation
(~90 Mya) is a much shallower proxy and we have no deeper alignment to test
theirs. The conclusion is about what **our** data can support:

> The brain claim cannot be made at the level of DCR-containing genes, in any
> conservation stratum. It can only be made about **specific ancient DCRs
> producing brain-specific transcripts or isoforms** — which is the
> frame/isoform-level claim, carried by the 108 significant tissue×frame
> interactions and corroborated at protein level by Tress for ATP2B1, ATP2B4,
> MADD and MARK4.

Write the brain sentence that way and it is defensible from both studies. Write
it as a property of DCR genes and it is contradicted by our own data.

**Limitation to state:** gene-level median TPM. This is not exon-level or
frame-level tissue specificity; an exon-level answer would require a matched
control run of Tallyman and should be offered only if the referee insists.
