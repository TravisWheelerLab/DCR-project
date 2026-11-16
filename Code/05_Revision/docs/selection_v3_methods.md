# Alternate-frame selection test, v3: reimplementation, v2 reproduction, non-primary-contig exclusion

Scripts: `selection/selection_v3.py` (functions + command-line driver; in this repository the former `run_selection_v3.py` driver is merged into it). Run: `python Code/05_Revision/selection/selection_v3.py --inputs-dir $WORK/selection_inputs --outdir $WORK/selection` (or `run_revision.sh selection`). *(Original note: the analysis was run interactively in the `dcr-tierA` environment by exec'ing the two files.)*
Runtime is about 1 min. Bootstrap: B = 2000 interval resamples per arm (1000 for length bins and the length-matched test), seed 20261006.

## 1. Reproducing v2 (hard requirement 1): reproduced within tolerance

The v2 code was lost, so this is a reimplementation. The full check is in `selection_v2_reproduction_check.csv`; every row passes except the length-matched null value (see deviation D4).

| quantity | v2 | v3 reimplementation |
|---|---|---|
| DCR n | 1,228 | **1,228** (661 ss13 / 567 ss12, exact) |
| DCR dNN/dNS | 0.2437 [0.2233, 0.2641] | **0.2437** [0.2245, 0.2656] |
| DCR dSN/dSS | 0.4465 [0.3757, 0.5292] | **0.4465** [0.3719, 0.5282] |
| DCR ω_ref | 0.2604 | **0.2604** |
| null n | 21,340 | 21,409 (+0.32%) |
| null dNN/dNS | 0.7981 | 0.7943 |
| null dSN/dSS | 0.9839 [0.951, 1.017] | 0.9831 [0.9511, 1.0168] |
| null ω_ref | 0.1097 | 0.1100 |
| ss13: DCR / null dSN/dSS | 0.4427 (661) / 0.8269 (11,559) | 0.4427 (661) / 0.8284 (11,611) |
| ss12: DCR / null dSN/dSS | 0.4511 (567) / 1.1955 (9,781) | 0.4511 (567) / 1.1921 (9,798) |

The DCR arm reproduces exactly, to every printed digit, including the frame-relation split. The null arm is within tolerance on every n and every point estimate. Its bootstrap CI matches v2 to the third decimal.

### Implementation as reproduced
- **Site classification.** OLGenie classes NN/NS/SN/SS (first letter: reference frame; second letter: alternate frame). Only positions whose codon is complete in both frames are scored. Each of the three alternative bases adds 1/3 site. Human–mouse differences are classified as single-site changes **in the human codon context**. (Averaging over human and mouse contexts gives DCR dSN/dSS = 0.461 rather than 0.4465, so v2 used human context.) A change to or from a stop counts as nonsynonymous. Counts are pooled over intervals, converted to p = diffs/sites, Jukes–Cantor corrected, and then ratioed.
- **DCR arm.** Human sequences come from `dcr_regions.fa` and mouse sequences from the UCSC liftOver set `D_ucsc.fa`. Frames come from `dcr_protein_frames.csv`; `prot_offs` are direct 0-based in-interval offsets (99.3% of protein frames are stop-free, against 39% of non-protein frames). **Reference = lowest offset, alternate = second offset**: this is the only assignment that gives 661/567. For a DCR the choice is arbitrary, so a swapped-role sensitivity is reported below.
- **Null arm.** Inputs are the v2 intermediates exported from the previous account (`09-16_DCR-revisions/`): `N_human.fa` (36,177 latent intervals), `N_ucsc.fa`, `N_cds.ov`, `cds_exons.tsv`, and `nondcr_pool.fa`. The reference frame is the GENCODE v48 frame from fully containing CDS blocks: cumulative preceding CDS length → codon position at the 5′ base → offset = (3 − cpos) % 3. Only intervals with exactly one annotated frame are used.
- **Latent-frame rule (reverse-engineered).** The rule is not documented anywhere, so it was recovered by matching v2's per-frame-relation counts. The reference frame must be stop-free. If exactly one other frame is stop-free, it is the alternate. If both other frames are stop-free, the alternate is the lowest offset when the reference offset is 0 or 1, and **the interval is dropped when the reference offset is 2** (1,558 intervals before other filters).
  - Every category-level assignment that meets the 1% n-tolerance requires that last drop. This rule also reproduces v2's per-length-bin null values far better than the alternatives: (0,30] 1.419 vs 1.434, against 1.024 under the ORF-pair rule without the drop.
  - It has no biological rationale and is almost certainly an artefact of the lost code. It is kept in the main v3 numbers so that v3 differs from v2 only by the contig exclusion. A principled alternative is reported as a sensitivity.
- **Filters.** Interval length < 30 nt, failed lift, and human/mouse length mismatch (indel). The indel rate among lifted null intervals is 14.98%, which matches v2's quoted 15.0%.

## 2. v3: non-primary contigs excluded from both arms

Excluded contigs: `_alt`, `_random`, `_fix`, `chrUn*`. Among intervals that v2 scored, this removes **13 DCR intervals** (all `_alt`: chr6 MHC haplotypes GL000251/252/253, chr10_GL383546, chr17_KI270861, chr11_KI270831, chr19_GL949749, chr1_GL383520) and **1,461 null intervals** (1,453 `_alt`, 7 `_random`, 1 `chrUn`).

| | DCR (n = 1,215) | latent null (n = 19,948) | null − DCR [95% CI] | p |
|---|---|---|---|---|
| ω_alt = dNN/dNS | **0.2427** [0.2238, 0.2630] | 0.7932 [0.7817, 0.8051] | 0.5505 [0.5268, 0.5736] | < 0.0005 |
| ω_alt = dSN/dSS | **0.4430** [0.3718, 0.5216] | 0.9881 [0.9550, 1.0236] | 0.5451 [0.4577, 0.6244] | < 0.0005 |
| ω_ref = dNN/dSN (positive control) | 0.2605 [0.2383, 0.2837] | 0.1050 [0.1029, 0.1074] | | |

p is the two-sided bootstrap p-value: none of the 2000 replicates of the difference crossed 0, so p < 1/B.

**By frame relation** (`selection_by_framerel_v3.csv`):

| | DCR dSN/dSS | null dSN/dSS | diff [95% CI] |
|---|---|---|---|
| ss13 (+1) | 0.4383 [0.3506, 0.5452], n = 653 | 0.8371 [0.7982, 0.8782], n = 10,845 | 0.399 [0.282, 0.497] |
| ss12 (+2) | 0.4486 [0.3487, 0.5827], n = 562 | 1.1903 [1.1314, 1.2489], n = 9,103 | 0.742 [0.596, 0.859] |

dNN/dNS is in the same file: ss13 0.243 vs 0.689; ss12 0.243 vs 0.937. Both differences have p < 0.0005.

**Length-matched** (`selection_length_matched_v3.csv`). The null is reweighted to the exact DCR length distribution. 1,211 DCRs have at least one null interval of identical length.
- dSN/dSS: DCR 0.4344 [0.3645, 0.5166] vs null 0.9899. A bootstrap that redraws 10 length-identical null intervals per DCR in each replicate gives a null CI of [0.9345, 1.0466] and a difference CI of [0.4585, 0.6463], p ≤ 0.001.
- dNN/dNS: 0.2405 vs 0.7801.

**Per-length-bin dSN/dSS** (`omega_by_length_bin_v3.csv`; CIs from 1000 bootstrap replicates):

| bin (nt) | DCR | null |
|---|---|---|
| (0,30], i.e. exactly 30 nt | 0.165 [0.045, 0.708], n = 21 | 1.123 [0.422, 3.300], n = 113 |
| (30,45] | 0.379 [0.202, 0.628], n = 243 | 1.140 [0.934, 1.392], n = 1,729 |
| (45,60] | 0.449 [0.276, 0.715], n = 211 | 1.018 [0.888, 1.162], n = 2,746 |
| (60,90] | 0.374 [0.259, 0.536], n = 357 | 0.975 [0.904, 1.048], n = 6,047 |
| (90,120] | 0.688 [0.473, 1.032], n = 181 | 1.044 [0.969, 1.123], n = 4,248 |
| (120,180] | 0.348 [0.240, 0.507], n = 150 | 0.926 [0.860, 0.989], n = 3,559 |
| (180,300] | 0.504 [0.301, 0.858], n = 48 | 0.988 [0.897, 1.087], n = 1,185 |
| > 300 | 0.835 [0.253, 2.293], n = 4 | 1.022 [0.899, 1.158], n = 321 |

The DCR point estimate is below the null in every bin. However, **the per-bin 95% CIs overlap in (90,120]**, as well as in the sparse (0,30] and > 300 bins. v2's wording that DCR and null are "separated in every length bin" should become "the DCR estimate is below the null in every bin". The per-bin CIs are new in v3; v2 reported point values only.

**Sensitivities** (rows in `selection_omega_v3.csv`, v3 scope):
- *Principled null* (`sens_null_orfpair_latent`). The alternate frame is the frame of the ORF pair whose intersection defines the latent interval in its source exon, and the reference-offset-2 drop is not applied (n = 20,503). Null dSN/dSS = 0.982 [0.951, 1.016]; the difference is 0.539 [0.452, 0.620]. The v2 latent-frame quirk does not affect the conclusion.
- *DCR frame roles swapped* (`sens_DCR_roles_swapped`; reference = higher offset). DCR dSN/dSS = 0.4755 [0.4015, 0.5589] and dNN/dNS = 0.2605. The difference from the null is 0.513 [0.421, 0.598]. The result does not depend on which DCR frame is called "canonical".

**Effect of the exclusion:** negligible. DCR dSN/dSS moves from 0.4465 to 0.4430 and null from 0.9831 to 0.9881. The difference moves from 0.537 to 0.545. The null positive control moves from 0.110 to 0.105.

## 3. Exclusion fates (`dcr_exclusion_fates_v3.csv`, all 1,741 intervals)

- **v3 cascade** (non-primary → too short → lift failed → indel): kept 1,215; too_short 271; indel 185; lift_failed 46; nonprimary_contig 24.
- **v2-scope cascade:** kept 1,228; too_short 275; indel 190; lift_failed 48.
- Every individual flag is a separate column, so any other ordering can be rebuilt.

**Why the counts differ from the old table** (245 / 201 / 49): the old table tested lift and indel *before* length. In the old table, 23 sub-30-nt intervals are labelled `indel` (15) or `lift_failed` (8). Here they are `too_short`. The scored set is identical.

**The 24 primary-contig intervals with no row in the old fates table.** The old table (`dcr_exclusion_fates.csv`, 1,717 rows) was built in v1 from `D_cds.ov`, the overlap of DCR intervals with GENCODE v48 CDS blocks. Its 1,717 IDs are exactly the 1,717 intervals present in `D_cds.ov`. The 24 missing intervals have no qualifying GENCODE CDS block:
- 21 overlap no GENCODE v48 CDS at all.
- 3 overlap only partially: dce2|495 (36/43 nt), dce2|1106 (301/456 nt), dce2|856 (13/60 nt).

All 24 are present in the protein-isoform overlap table (`D_iso.ov`), so they are DCRs by protein evidence. They are the extreme case of the GENCODE-vs-protein authority error fixed in v2. v2 scored them silently (1,228 = 956 + 266 + **6** of these 24) but never added them to the fates table. Their fates now: kept 6, too_short 7, lift_failed 7, indel 4. None is on a non-primary contig. Column `old_table_missing_reason` gives the reason per interval.

The per-interval site-count table `selection_site_counts_v3.csv` (37,918 rows: 1,741 D + 36,177 N) contains:
- coordinates, frames, flags and fates;
- the 8 OLGenie counts (sites and differences per class);
- swapped-role counts (D) and orfpair-rule counts (N);
- `scored_v2` / `scored_v3` masks and a length bin.

Summing the count columns over a mask reproduces every number above.

## Deviations from v2 (explicit)

- **D1.** The code is a reimplementation. The null latent-frame rule was **reverse-engineered** by fitting v2's counts, not taken from documentation. That includes the dropping of 1,558 three-open intervals with reference offset 2. The null therefore has 21,409 rather than 21,340 intervals (+69; ss13 +52, ss12 +17). This is within tolerance but not identical, so the null's exact membership is a reconstruction.
- **D2.** The null inputs are the v2 intermediates recovered from the previous account's export folder, not regenerated from DOOVER. `cds_exons.tsv` (93 MB) and `nondcr_pool.fa` were read from that folder and are not re-saved as artifacts.
- **D3.** Bootstrap CIs differ from v2 in the third decimal because of the reseeded bootstrap (B = 2000). v2's seed is unknown.
- **D4.** The **length-matched null does not reproduce**: 0.9906 here (exact-length reweighting) vs v2's 1.034. The v2 matching procedure is unknown and was probably a single random matched draw. v2's value lies inside this analysis's matched-bootstrap null CI [0.938, 1.049]. The DCR side does reproduce (0.438 vs 0.439).
- **D5.** Per-bin null values match v2 within 0.02 in 5 of 7 bins. The exceptions are (0,30] (1.419 vs 1.434; 0.015, within tolerance but on n = 127) and **(30,45], 1.142 vs 1.204**, which is inside that bin's CI. Bin n's differ from v2 by ≤ 1%, except (90,120] (+0.8%) and (180,300] (+1.1%).
- **D6.** New in v3: bootstrap CIs for the null and for per-bin values, a bootstrap length-matched test, and two sensitivity analyses.
- **D7.** The exclusion list is applied exactly as specified. **chrM is not excluded**, so 2 chrM null intervals remain in both scopes even though they are scored with the nuclear genetic code. Their effect on a 19,948-interval pool is negligible, but they could be removed in a later pass. `_fix` contigs contributed 0 scored intervals.
- **D8.** Fate cascade order changed: length is tested first. See section 3.

No subsampling was done; all intervals in both arms were used.
