# Code/05_Revision — analyses added for the GBE revision (2026)

Code behind the new and revised analyses of the revised manuscript (Goubert et al., *Alternatively
spliced dual-coding regions…*, bioRxiv 10.64898/2026.02.24.707805, revision for GBE). Every analysis
starts from the outputs of the executed pipeline of `Code/01–04` (called **DOOVER** below:
`output/01_Description`, `output/02_Conservation`, `output/New_TableS1.tsv`, `data/`). DCR calls,
frames and isoform models always come from the protein-level (UniProtKB/Swiss-Prot isoform) mappings
of that pipeline, never from GENCODE CDS annotation. Nothing is written inside DOOVER; all outputs go
to `$WORK`.

## Layout

| Path | Content |
|---|---|
| `config.sh` | all paths (DOOVER, reference files, `$WORK`, interpreters); override any by environment variable |
| `run_revision.sh` | driver; runs the steps below in dependency order (`./run_revision.sh [step ...]`) |
| `environment.yml`, `requirements.txt`, `R_packages.txt` | software, pinned to the versions used |
| `selection/` | alternate-frame (OLGenie-style) selection test v3 + input staging/md5 manifest |
| `splicing_vs_AS/` | DCRs vs ordinary alternative splicing (length, frame preservation, NMD, conservation) |
| `controls/` | splice-site phyloP; GTEx expression class / tissue specificity (tau) / brain dominance; GABBR1–PTBP2 ascertainment |
| `figures/` | Figure 2B and Supplementary Figures S6–S9 (R, house style `dcr_style.R`) + their data-prep scripts |
| `utils/` | pure-Python replica of `esl-translate -l 1 --watson` |
| `results/` | small result tables (< 1 MB each) read by the figures; see `results/README.md` |
| `docs/` | methods/results notes written during the revision and the draft figure legends |
| `archive_prior_session/` | **as executed, not refactored**: verbatim code cells of the July–August 2026 session (Tier A controls, v1/v2 selection test) |

## Quick start

```bash
conda env create -f Code/05_Revision/environment.yml && conda activate dcr-revision   # + R packages
export DOOVER=/path/to/DOOVER                       # executed pipeline directory
export REFDIR=/path/to/reference/data               # hg38.fa(+.fai), hg38.phyloP100way.bw, ...
export PHYLOP_BW=/local/disk/hg38.phyloP100way.bw   # local copy! (see below)
export PRIOR_EXPORT=/path/to/exported/revision/artifacts   # selection-test intermediates
export WORK=$PWD/revision_work
Code/05_Revision/run_revision.sh                    # everything
Code/05_Revision/run_revision.sh selection figure_inputs figures   # or a subset
```

Every script also runs on its own with explicit arguments (`python <script> --help`); defaults are
derived from `$DOOVER`, `$REFDIR` and `$WORK`.

**phyloP bigWig.** `hg38.phyloP100way.bw` is 9.9 GB. Copy it to a local disk before running
`splice_sites` or `splicing_vs_AS`: random reads over a cloud-synced drive ran at about 1.6 queries/s
(≈ 11 h for the 64,064 splice sites), against ~10 s from a local SSD. The same applies to `hg38.fa`.

## Steps (dependency order)

| Step | Script(s) | Inputs | Main outputs (`$WORK/…`) | Manuscript | Runtime* |
|---|---|---|---|---|---|
| `stage_selection` | `selection/stage_selection_inputs.py` | `$PRIOR_EXPORT` | `selection_inputs/` (11 files, md5-checked) | — | < 1 min |
| `selection` | `selection/selection_v3.py` | `selection_inputs/` | `selection/selection_omega_v3.csv`, `…_by_framerel_v3.csv`, `…_length_matched_v3.csv`, `…_by_lenbin_v3.csv`, `selection_v2_reproduction_check.csv`, `selection_site_counts_v3.csv` | **Figure 2B(iii)**; Methods/Results, selection test | ~1 min |
| `ascertainment` | `controls/gabbr1_ascertainment.py` | `selection_inputs/`, DOOVER DCR bed | `controls/selection_interval_fates_with_coords.csv`, `controls/A3_showcase_candidates.csv` | **Fig S6** (GABBR1) legend; Results | < 1 min |
| `splice_sites` | `controls/splice_site_phylop.py` | DOOVER splice-site beds, phyloP bigWig | `controls/splice_site_phylop.tsv`, `controls/splice_site_conservation.csv` | **Fig S9 A** | ~10 s (local bigWig) |
| `gtex_fetch` | `controls/gtex_fetch_medians.py` | DOOVER DCR + isoform beds; GTEx Portal API v2 (network) | `controls/gtex_v8_gene_medians.json`, `controls/a4_genes.json` | **Fig S9 B–C** | ~3 min |
| `expression_tau` | `controls/expression_tau.py` | the JSON above; DOOVER `DCRmiceStatus.tsv`, `DCEtoGene.txt` | `controls/A4_dcr_by_expression.csv`, `A4_tau_by_expression.csv`, `A4_gene_expression_tau.csv`, `A4_brain_by_conservation.csv`, `A4_tests.txt` | **Fig S9 B–C**; Results | < 1 min |
| `splicing_vs_AS` | `splicing_vs_AS/A2_build_events.py` → `nmd_allgenes.py` (×`$NMD_CHUNKS`) → `A2_event_tables.py` → `A2_stats.py` | DOOVER `output/01_Description/*`, hg38.fa, phyloP | `splicing_vs_AS/A2_events_annotated.csv`, `A2_table_*.csv` | **Fig S10**; Results/Methods (DCR vs ordinary AS) | ~5 min (8 NMD chunks) |
| `figure_data` | `figures/make_gc_tables.py`, `figures/prep_gene_models.py` | DOOVER beds/FASTAs, hg38.fa | `figure_data/gc_per_*.tsv`, `gm_*.csv`, `gnas_mouse_*.csv` | Fig S6, S7, S8 | ~10 s (local hg38) |
| `latent_check` (optional) | `figures/verify_latent_exon_fraction.py` | DOOVER exon beds, hg38.fa | `figure_data/latent_exon_fraction.csv` | check of Fig 2B(i) | ~20 s |
| `figure_inputs` | (driver) | `results/*.csv` + fresh `$WORK` tables | `figure_inputs/` | — | < 1 s |
| `figures` | `figures/*.R`, `splicing_vs_AS/A2_figure.R` | `figure_inputs/`, `figure_data/`, `splicing_vs_AS/` | `figures/*.pdf/png` + `*_plotted.csv` / `*_check.csv` | Fig 2B, S6–S10 | ~10 s |

\*Measured on the analysis laptop (Apple silicon, 10 cores) with hg38.fa and the phyloP bigWig on local disk.

### Figure ↔ script map

| Figure | Script | Output stem |
|---|---|---|
| Figure 2B(i)–(iii) (dual-coding potential vs. selection) | `figures/fig2B_potential_vs_selection.R` | `Figure2B_dual_coding_potential_vs_selection` |
| Supplementary Figure S6 (GABBR1 DCR gene model) | `figures/supp_GABBR1_gene_model.R` (+ `gene_model.R`) | `Supp_GABBR1_DCR_gene_model` |
| Supplementary Figure S7 (GNAS DCR, mouse ORFs) | `figures/supp_GNAS_gene_model.R` (+ `gene_model.R`) | `Supp_GNAS_DCR_gene_model_mouse_ORFs` |
| Supplementary Figure S8 (GC-content controls) | `figures/supp_gc_content_controls.R` | `Supp_GC_content_controls` |
| Supplementary Figure S9 (splice sites, expression class, tissue specificity) | `figures/supp_splice_expression_tissue.R` | `Supp_splice_site_expression_tissue_specificity` |
| Supplementary Figure S10 (DCRs vs ordinary alternative splicing) | `splicing_vs_AS/A2_figure.R` | `A2_FigureS_DCR_vs_regularAS` |

Legends: `docs/figure_legends.md` (draft; Supplementary figure numbers assigned in the manuscript).

## The analyses

**Alternate-frame selection test (v3), `selection/`.** OLGenie-style site classification
(Nelson, Ardern & Wei 2020) of human–mouse differences in regions where two frames are both
complete: dNN/dNS and dSN/dSS for the alternate frame, dNN/dSN for the reference frame (positive
control), pooled over intervals, Jukes–Cantor corrected, 2,000 bootstrap replicates over intervals.
DCR arm: catalogue DCR intervals with both frames taken from the Swiss-Prot isoform models; null arm:
latent dual-coding intervals in non-DCR coding exons (GENCODE v48 reference frame, open alternate
frame). v3 reimplements the lost v2 code (reproduces v2: DCR arm exactly, null n 21,409 vs 21,340;
`selection_v2_reproduction_check.csv`) and removes non-primary contigs from both arms. Headline
(primary contigs): DCR n = 1,215, dSN/dSS = 0.443 (95% CI 0.372–0.522); null n = 19,948,
dSN/dSS = 0.988. Details: `docs/selection_v3_methods.md`, `docs/OLGenie_methods_explained.md`.
The 11 inputs are intermediates built by the archived session (`archive_prior_session/10_…`, `12_…`);
`stage_selection_inputs.py` finds them in an exported-artifact folder and checks their md5.

**DCRs vs ordinary alternative splicing, `splicing_vs_AS/`.** Atomic coding segments are derived
from the Swiss-Prot isoform models (frame from the protein model itself) and classified as
dual-coding, ordinary internal AS (frame-preserving or frame-shifting), constitutive or variable
terminal; DCRs are compared with ordinary AS segments of DCR genes and of non-DCR genes for length,
frame preservation, NMD (the pipeline's verbatim `analyze_isoforms()` rerun on all genes;
validated against the pipeline's `NMDtable.tsv`), and phyloP conservation (length/GC- and
inclusion-matched). Details: `docs/A2_results.md`.

**Splice-site conservation, `controls/splice_site_phylop.py`.** Mean phyloP100way at the donor and
acceptor sites flanking DCR exons vs control exons (the pipeline's splice-site beds), per junction ×
site, Mann–Whitney U with rank-biserial effect size, Bonferroni-adjusted.

**Expression class and tissue specificity, `controls/gtex_fetch_medians.py` + `expression_tau.py`.**
GTEx v8 median TPM in 54 tissue sites per gene (GENCODE v26 ids resolved from symbols; DCR genes vs
all genes with mapped Swiss-Prot isoforms), DCR incidence per expression quintile (+ logistic
regression), tau within quintiles, and the fraction of brain-dominant genes in DCR genes split by
mouse conservation. During the revision the GTEx data came through an in-platform connector;
`gtex_fetch_medians.py` uses the public GTEx Portal API v2 instead (`--tie-break lowest`
reproduces the connector's gene set exactly). Details: `docs/A4_results.md`.

**GABBR1/PTBP2 ascertainment, `controls/gabbr1_ascertainment.py`.** Joins the selection-test
exclusion fates to interval coordinates and lists the candidate showcase genes (GABBR1's 26-nt DCR
falls below the 30-nt filter). Details: `docs/A3_results.md`.

**Figures, `figures/`.** House style in `dcr_style.R` (Helvetica, vector PDF via the base `pdf()`
device + 300-dpi PNG via ragg; fonts embedded with Ghostscript when available); `gene_model.R`
draws Figure 2E-style gene models. `make_gc_tables.py` and `prep_gene_models.py` derive the per-exon
GC tables and the GNAS/GABBR1 model tables from DOOVER; every figure writes the values it plotted
(`*_plotted.csv` / `*_check.csv`) and asserts consistency with the summary tables.

**Archived prior-session code, `archive_prior_session/`.** The Tier A GC audit, GC+length-matched
controls, DCR self-shuffle, phyloP with GC control, latent potential, GNAS, definition
reconciliation, X/X′ conservation filter, annotation-depth check and the v1/v2 selection test ran as
interactive cells in an earlier session. Their code is archived verbatim (paths anonymised), labelled
"as executed, not refactored"; their result tables are in `results/`.

## Reproduction status (smoke test, October 2026)

Run from this layout against DOOVER, compared with the tables and figures used in the manuscript:

| Step | Result of the re-run vs. the saved manuscript outputs |
|---|---|
| `stage_selection` | 11/11 inputs found in the prior-session export folder, md5 match |
| `selection` | 5/5 tables byte-identical (`selection_omega_v3`, `_by_framerel_v3`, `_length_matched_v3`, `selection_v2_reproduction_check`, `selection_site_counts_v3`). Primary contigs: DCR n = 1,215, dSN/dSS = 0.443; null n = 19,948, dSN/dSS = 0.988 |
| `ascertainment` | 2/2 tables byte-identical |
| `splice_sites` | 2/2 tables byte-identical (64,064 sites, 64,032 with phyloP) |
| `gtex_fetch` + `expression_tau` | public API run: 16,532 genes (1,230 DCR genes), same gene set and per-gene median/max/tau/top tissue as the connector data; 3/4 tables byte-identical, `A4_gene_expression_tau.csv` identical up to row order. Fed the connector JSON, all 4 are byte-identical |
| `splicing_vs_AS` | 10/11 `A2_table_*.csv` and `A2_events.tsv`, `A2_events_annotated.csv`, `A2_matched_pairs.csv`, `segments.tsv`, `isoform_frames.tsv` byte-identical. `NMDtable_allgenes.tsv` differs for 42 isoforms, all in single-merged-exon genes whose call depends on gene order (documented in `docs/A2_results.md`; these genes are excluded from every comparison), so the 11th table, `A2_table_NMD_validation.csv`, reports `rerun_NMD` = 3,778 instead of 3,777 |
| `figure_data` | 11/11 tables with a saved copy byte-identical (`gc_per_region_mouse.tsv` had none) |
| `latent_check` | `latent_exon_fraction.csv` byte-identical |
| `figures` | all 6 figure PNGs (Fig 2B, S6–S10) byte-identical; plotted-value tables identical (Fig 2B table differs only in the recorded input path) |

## Known limitations

* The archive is a log, not a pipeline: it cannot be re-run as is.
* The selection-test inputs are prior-session intermediates (not rebuilt from DOOVER by a script in
  this directory); they are identified by md5 in `selection/selection_inputs_manifest.tsv`.
* `gtex_fetch_medians.py` depends on the live GTEx Portal API; the result is identical to the
  connector-fetched values used for the manuscript as of October 2026 (`A4_gene_expression_tau.csv`
  differs only in row order).
* Pre-existing issues of `Code/01–04` are not addressed here (see the PR description).
