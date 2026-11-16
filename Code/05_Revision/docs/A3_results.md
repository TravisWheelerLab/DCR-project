# A3 — GABBR1 / PTBP2 ascertainment check

Answers R2-m1 and tests whether Tress et al.'s non-recovery of these two regions
is explained by their ascertainment. Source: the executed pipeline
(`filtered_DCR_overlapping_regions_with_gene_and_dceID.bed`, `New_TableS1.tsv`)
joined to our own selection-test exclusion table by genomic coordinate.

## Result

| gene | interval | len | frames | our selection test | contig | GC |
|---|---|---|---|---|---|---|
| **GABBR1** | chr6_GL000255v2_alt:874219-874245 | 26 nt | 2 | **too_short — excluded** | **ALT haplotype** | 46.2% |
| **PTBP2** | chr1:96806865-96806889 | 24 nt | 2 | **too_short — excluded** | primary | 50.0% |

Both regions fall below the 30-nt floor of our alternate-frame selection test, so
**neither is supported by our own strongest evidence.** This has to be stated
plainly; it is not a detail a careful referee would miss.

**Evidence tier on the alternate isoform — CORRECTED.** An earlier version of
this note listed only isoforms -2/-3/-4 and said every alternative product
carries only `ECO:0000303`/`ECO:0000305`. That omitted the -5 isoforms and was
wrong. From `New_TableS1.tsv`:

| gene | isoform | frame role | isoform evidence | ECO codes |
|---|---|---|---|---|
| GABBR1 | Q9UBS5-1 | canonical | Protein-level-evidence(1) | …, **0000269**, … |
| GABBR1 | Q9UBS5-2/-3 | canonical | No-evidence-indicator | 0000303 |
| GABBR1 | Q9UBS5-4 | canonical | No-evidence-indicator | 0000305 |
| GABBR1 | **Q9UBS5-5** | **alternative (+2)** | No-evidence-indicator | **0000269**, 0000305 |
| PTBP2 | Q9UKA9-1 | canonical | Protein-level-evidence(1) | …, 0000269, … |
| PTBP2 | Q9UKA9-2/-3/-4/-6 | — | No-evidence-indicator | 0000303 |
| PTBP2 | Q9UKA9-5 | — | No-evidence-indicator | **0000269**, 0000303 |

`ECO:0000269` is *experimental evidence used in manual assertion*. So the
GABBR1 alternative-frame isoform has curated experimental support, even though
UniProt's isoform-level indicator is not "protein level". For R2-m1 this is
**better** than the earlier version: the alternative frame is experimentally
documented, not just predicted. It also means the peptide-route argument about
Tress's ascertainment can't be built on SwissProt tags. An experimental
assertion is not the same as a PeptideAtlas peptide, and only a direct
PeptideAtlas query would settle it. **Drop the ascertainment argument for
GABBR1.** Describe it on its own evidence, which is what the user wanted
anyway. Check the current UniProt entry for Q9UBS5-5 before the evidence
statement goes into the legend, since the table is a snapshot.

**Conclusion for the response letter.** Their non-recovery is consistent with
their design: no protein-level evidence on the alternate isoform closes the
peptide route, and these are short internal regions rather than the C-terminal
frame changes their screen targets. But the honest framing is **a region that
evades both methods** — their peptide-ascertained curation and our length-gated
selection test — not a region we can show is under selection.

**Caveat to disclose:** GABBR1's DCR sits on an ALT haplotype contig in the MHC
region, not primary chr6. Its GC (46.2%) is below the DCR median, so it is not a
GC artifact, but the coordinate placement should be stated rather than left for
a referee to find.

## Showcase choice

The ancient brain-specific DCRs Tress confirmed split sharply by length, and the
split runs straight through our 30-nt floor:

| fate | genes |
|---|---|
| **kept** (in our selection test) | CCSER2 92 nt · PLCB4 66 nt · MARK4 62 + 45 nt · MAGI1 54 nt · PBX1 44 nt |
| **too_short** (excluded) | GABBR1 26 · PTBP2 24 · ATP2B4 23 · ATP2B1 23 · MADD 21 nt |

**Use CCSER2 as the pillar-1 worked example.** 92 nt, primary chr10, GC 54.3%,
`kept` in our selection set — and Tress independently report CCSER2 as a
triple-coding region under selection in **all three** frames. It is the one
region where our population-level test and their per-region curation both apply
and agree, which makes it both the strongest example and a natural place to
credit them.

Keep GABBR1 as the detailed R2-m1 answer (isoforms, exon, frames, the 26-nt
interval, the ALT-contig caveat); do not ask it to carry the selection argument.

## Incidental finding worth one Methods line

The 24 interval rows on non-primary contigs (20 DCRs, 20 genes: aadacl2, alox5,
atp6v1g3, c17orf78, cars1, ccl4l1, ddr1, gabbr1, gp6, gpx5, ltb, mboat7, nrm,
phrf1, rpp21, scarf1, slc2a11, tbc1d3, trim31, vwa7) were **not** excluded from
the selection analysis — 8 of them are in the `kept` set. phyloP and multiz
coverage on ALT haplotypes is unreliable, so this is worth either excluding and
re-reporting (the effect on n = 948 kept intervals is negligible) or disclosing.
Separately, 24 intervals (all primary) have no row in the exclusion table at
all; unexplained, and also negligible, but it should not be discovered by a
referee.
