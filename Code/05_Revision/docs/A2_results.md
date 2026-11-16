# A2: DCRs compared with ordinary (non-dual-coding) alternative splicing (Referee 1, R1-M1)

## Deviations from the brief (declared up front)
1. **Event-level NMD rule.** The pipeline's DCR call (`analyze_DCR`) has two parts: a core rule (one reading-frame arm of isoforms is 100% NMD-predicted while another arm is not) and a refinement that requires the NMD arm to read to the C-terminus ('C-TERM'). The refinement depends on frame termini, which regular AS events do not have. So both classes were compared on the core rule, with inclusion vs exclusion isoforms as the two arms of a regular-AS event. Under the core rule 479/1,407 DCRs (34.0%) are NMD-susceptible; the official pipeline figure, reproduced exactly here, is 466/1,407 (33.1%). Both numbers are reported.
2. **Isoform-level NMD calls** come from the pipeline's `analyze_isoforms()`, copied verbatim and run on all genes rather than DCR genes only. Two harness changes were needed: the gene list was chunked, and the strand string was serialised in a way that does not depend on the pandas version. In one chunk, the first gene was moved so that the chunk does not start with a gene whose merged model is a single exon. Such genes make the verbatim code reuse `EJCindex` from the previous gene, so their call depends on gene order. All 1,249 single-merged-exon genes were therefore left out of the NMD comparisons. This removed 31 regular-AS events from non-DCR genes and 0 DCRs.
3. **Reading frames.** Each base's codon position was computed from each SwissProt isoform model, counting cumulative CDS length from the first codon. I could not decode the convention of the BED phase column, so it was not used. No GENCODE CDS was used. Validation: this recovers 100% of the bp of the 1,741 catalogue DCR intervals as dual-coding. 957 of 39,149 isoforms (2.4%) have an internal stop when read this way. They were kept, and excluding every event that involves one leaves all results unchanged (NMD 34.0 / 19.2 / 10.9%, median phyloP 2.95 / 2.26 / 2.39).
4. **Data reductions.** Regular-AS segments shorter than 6 nt (the shortest DCR interval) were dropped: 3,313 of 24,619. 772 "grey-zone" genes were excluded from every regular-AS arm: these have dual coding in the raw MIRAGE output or in the protein-model derivation but no catalogue DCR. The constitutive non-DCR context pool (123,702 segments) was subsampled to 20,000 for nearest-neighbour matching only. Under the length/GC calipers, 1,402–1,407 of 1,407 DCRs found matches; under the inclusion-matched design, 448–476 of 478. 5 genes with conflicting exon strands are skipped by the pipeline code itself, so NMD calls exist for 39,120 of 39,149 isoforms. Two fs-AS segments with <90% phyloP coverage were dropped.
5. **Not addressed here:** protein structure/disorder (out of scope; the paper already has AlphaFold3 on 56 isoform pairs) and "adaptiveness". The latter would need polymorphism or divergence-based selection tests, which this analysis does not do.

## Definition of the reference class
- **Source.** The same SwissProt isoform models (MIRAGE2 mappings) that the DCRs come from: `output/01_Description/UniprotKB_isoforms_mapping.bed`, with the pipeline's 70 ambiguous-chromosome genes removed. That leaves 39,149 isoforms in 17,403 genes. All isoform models start at their first codon (96.8% start with ATG) and exclude the stop codon.
- **Atomic segments.** For each gene, every isoform exon boundary was used to cut the coding territory into atomic intervals. For each interval I recorded the set of isoforms that code it and each isoform's codon position. Adjacent intervals with the same coverage set and class were merged. This gives 206,446 segments:
  - **dual-coding**: ≥2 codon positions among the coding isoforms. 3,393 segments; these contain 100% of catalogue DCR bp.
  - **constitutive**: coded by every isoform in a single frame. 134,338.
  - **variable-terminal**: absent only from isoforms that start or end before the segment (alternative termini or truncations). 44,096; not used.
  - **regular internal AS**: single frame, coded by some isoforms, and skipped by at least one isoform that has coding exons on both sides. 24,619.
- **Arms** (segments ≥6 nt):
  - **D**: the 1,407 catalogue DCRs (1,741 intervals, 1,296 genes); intervals pooled per DCR.
  - **fp-D\***: frame-preserving regular AS (length a multiple of 3) in DCR genes. 1,467 events, 700 genes.
  - **fp-ND**: frame-preserving regular AS in non-DCR genes. 9,121 events, 5,075 genes. Non-DCR genes means no catalogue DCR, no raw MIRAGE DCR call and no derived dual-coding segment (15,343 genes).
  - **fs-D\*** and **fs-ND**: internal single-frame AS whose length is not a multiple of 3. These are frameshifting, but the shifted frame is never shared with another isoform, so no DCR forms. 2,574 events (1,242 genes) and 5,637 events (3,100 genes). They are kept as a second reference because they separate the effect of frameshifting from the effect of dual coding.
  - **c-D\*** and **c-ND**: constitutive segments, used as context for conservation. 6,757 and 123,702.
  - fp + fs together form the brief's "regular AS in DCR genes" (4,041) and "regular AS in non-DCR genes" (14,758).
- **Inclusion fraction** is the share of the gene's isoforms that code the segment. 929/1,407 DCRs (66%) are coded by every isoform of their gene. Regular AS cannot be in that state, by definition.

## (c) Length and frame preservation (descriptive)
| arm | n | genes | median length (IQR), nt | % length multiple of 3 |
|---|---|---|---|---|
| D | 1,407 | 1,296 | 70 (40.5–119.5) | 36.0 |
| fp-D\* | 1,467 | 700 | 81 (39–132) | 100 (definition) |
| fp-ND | 9,121 | 5,075 | 81 (39–129) | 100 (definition) |
| fs-D\* | 2,574 | 1,242 | 86 (49–137) | 0 (definition) |
| fs-ND | 5,637 | 3,100 | 80 (37–136) | 0 (definition) |
| c-D\* / c-ND | 6,757 / 123,702 | 1,171 / 15,112 | 120 / 124 | 38.8 / 40.8 |

- **Length.** DCRs are about as long as regular AS. Cliff's δ vs fp-D\* = −0.03 (p = 0.16) and vs fp-ND = −0.02 (p = 0.33). The only difference is vs fs-D\*: δ = −0.09, p = 1e-6, negligible.
- **Frame preservation, all internal regular AS.** 36.3% of regular AS in DCR genes has a length that is a multiple of 3, the same as DCRs (36.0%; OR 0.99, p = 0.85). In non-DCR genes the figure is 61.8% (DCR vs non-DCR-gene AS: OR 0.35, p = 1e-77). DCR genes are therefore enriched for frameshifting AS in general, not only at the DCR itself.

## (a) NMD susceptibility (pipeline 50-nt rule, identical code in all arms)
- **Validation of the rerun.** It reproduces the pipeline's `NMDtable.tsv` exactly: prediction and distance are identical for 5,053/5,053 isoforms. Over all genes, 3,777/39,120 isoforms are NMD-predicted. The official DCR result is also reproduced: 466/1,407 `in_NMD_zone`.

| arm | n events | NMD-susceptible | % | NMD arm = inclusion / exclusion |
|---|---|---|---|---|
| D | 1,407 | 479 | 34.0 | n/a (frame arms) |
| fp-D\* | 1,467 | 279 | 19.0 | 182 / 97 |
| fp-ND | 9,090 | 1,012 | 11.1 | 869 / 143 |
| fs-D\* | 2,574 | 716 | 27.8 | 281 / 435 |
| fs-ND | 5,637 | 971 | 17.2 | 909 / 62 |

DCR vs each reference. The adjusted model is a GEE logistic regression clustered by gene, adjusted for log length, log₂(isoforms per gene) and log₂(smallest arm size). The within-gene model is conditional logistic regression stratified by gene, adjusted for log length and arm size.
| reference | risk difference | Cohen's h | OR raw | OR adjusted (95% CI), p | OR within gene (95% CI), p |
|---|---|---|---|---|---|
| fp-D\* | +15.0 pts | 0.34 | 2.20 | 2.11 (1.78–2.50), 1.2e-17 | 3.75 (2.75–5.11), 7.4e-17 (246 informative genes) |
| fp-ND | +22.9 pts | 0.57 | 4.12 | 3.92 (3.41–4.52), 1.6e-80 | n/a |
| fs-D\* | +6.2 pts | 0.13 | 1.34 | 1.32 (1.22–1.44), 9.0e-11 | 2.34 (1.79–3.07), 7.0e-10 (241 genes) |
| fs-ND | +16.8 pts | 0.39 | 2.48 | 2.11 (1.83–2.43), 4.0e-25 | n/a |

- **Inclusion control.** Restricted to DCRs with inclusion fraction < 1 (n = 478, 38.3% NMD) and adjusting additionally for inclusion fraction, the excess holds: OR 2.48 (1.90–3.24) vs fp-D\*, 9.35 vs fp-ND, 1.25 (1.01–1.56, p = 0.043) vs fs-D\*, and 5.82 vs fs-ND.
- **By stratum.** DCR > fp-D\* > fp-ND in every inclusion stratum (`A2_table_NMD_by_inclusion.csv`).

## (b) Sequence conservation (mean phyloP 100-way per segment)
Medians: D 2.95; fp-D\* 2.27; fp-ND 2.40; fs-D\* 2.32; fs-ND 1.75; c-D\* 3.39; c-ND 3.65.

DCR vs each reference. Matching is 5 nearest neighbours on log₁₀ length and GC, with calipers 0.1 and 0.03. After matching, |SMD| < 0.011 for both covariates. The OLS model adjusts for length, GC, GC² and log₂(isoforms per gene), with gene-clustered standard errors.
| reference | raw Cliff's δ | matched median Δ (DCR − match) | matched rank-biserial | p (Wilcoxon) | OLS β (95% CI) |
|---|---|---|---|---|---|
| fp-D\* | +0.14 | +0.47 | +0.25 | 1.4e-16 | +0.42 (0.28–0.57) |
| fp-ND | +0.11 | +0.48 | +0.21 | 4.5e-12 | +0.39 (0.26–0.52) |
| fs-D\* | +0.13 | +0.60 | +0.30 | 6.3e-22 | +0.43 (0.29–0.56) |
| fs-ND | +0.21 | +0.78 | +0.41 | 1.0e-40 | +0.70 (0.57–0.84) |
| c-D\* | −0.14 | −0.17 | −0.07 | 0.019 | −0.38 (−0.49 to −0.27) |
| c-ND | −0.20 | −0.25 | −0.17 | 5.5e-8 | −0.69 (−0.79 to −0.58) |

**Inclusion level explains the raw advantage, and the sign reverses.** Conservation rises steeply with inclusion fraction in every class (Fig. E). 66% of DCRs sit at inclusion fraction 1, which regular AS cannot reach. Comparing like with like:
- **DCRs with inclusion < 1 (n = 478) vs regular AS matched on length, GC and inclusion fraction** (caliper 0.1). DCRs are *less* conserved:
  - fp-D\*: Δ −0.36, r_rb −0.21, p = 1.3e-4
  - fp-ND: Δ −0.27, r_rb −0.15, p = 0.005
  - fs-D\*: Δ −0.27, r_rb −0.13, p = 0.014
  - fs-ND: Δ −0.20, r_rb −0.10, p = 0.056
- **DCRs coded by every isoform (n = 929) vs constitutive segments matched on length and GC:**
  - c-D\* (same gene class): indistinguishable, Δ −0.05, r_rb 0.00, p = 0.93
  - c-ND: lower, Δ −0.19, r_rb −0.11, p = 0.003
- **OLS with inclusion fraction added** (extrapolating, because the two classes barely overlap in inclusion): DCR β becomes −0.57 vs fp-D\* (p = 5e-8) and −0.70 vs fs-D\* (p = 1e-13).
- **Within-gene paired comparison** (no inclusion control): DCR minus regular-AS gene means +0.39 (700 genes, p = 4e-17).

## Interpretation (draft paragraph for the response letter)
We built a reference set of ordinary alternative splicing from the same curated SwissProt isoform models. It contains 24,619 internal, single-frame coding segments that some isoforms include and others skip. From these we compared, with identical code, 1,407 DCRs against frame-preserving and frameshifting regular AS events, in DCR genes and in genes with no dual coding.

DCRs are clearly more NMD-prone than regular AS. With the pipeline's own 50-nt rule, 34.0% of DCRs are NMD-susceptible, against 19.0% of frame-preserving AS in the same genes and 11.1% in non-DCR genes. The odds ratio is 2.1–3.9 after adjusting for length, isoform number and arm size, and 3.7 within genes. DCRs also exceed frameshifting AS events that do not create a second shared frame (27.8% and 17.2%). This is consistent with a substantial share of DCR isoforms being NMD targets, as Benitiere et al. describe for low-abundance splice variants.

For conservation, DCRs look more conserved than regular AS at face value: +0.47 phyloP after length/GC matching. That difference comes entirely from inclusion level, because two-thirds of DCRs are coded by every isoform of their gene. At equal inclusion, DCRs are modestly *less* conserved than regular AS (Δ −0.2 to −0.4; rank-biserial −0.10 to −0.21). DCRs coded by all isoforms are as conserved as constitutive coding segments of the same genes (Δ ≈ 0).

DCR and regular-AS segments have similar lengths, and the proportion that preserves frame (36%) is identical to that of all regular AS in DCR genes. Overall, at the population level DCRs do not behave as a more constrained class than ordinary AS. Their conservation is at or slightly below what their inclusion level predicts, and their NMD susceptibility is higher. This agrees with our earlier finding that DCRs are less conserved than non-DCR exons (Fig. 2C), and with a view in which many DCRs are weakly constrained, while the conserved, mouse-shared subset (Fig. 2A, D) is the candidate functional fraction.

Protein structure was addressed separately (AlphaFold3, 56 pairs). We did not test adaptiveness, which would require polymorphism- or divergence-based selection tests.

## Files
- **Scripts:**
  - `A2_build_events.py`: segments and frames
  - `nmd_allgenes.py`: verbatim `analyze_isoforms`
  - `A2_event_tables.py`: arms, phyloP, GC
  - `A2_stats.py`: all tests
  - `A2_figure.R`: figure
- **Tables:**
  - `A2_events_annotated.csv`: one row per event, all arms
  - `A2_table_*.csv`
  - `A2_matched_pairs.csv`
  - `A2_NMDtable_allgenes.tsv`
