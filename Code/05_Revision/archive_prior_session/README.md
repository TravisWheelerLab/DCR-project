# archive_prior_session/ — code as executed, not refactored

These files hold, verbatim and in execution order, the code cells of the July–August 2026 revision
session that ran under a previous account (exported session `31323eb8_2026-09-16T22-37-01.json`,
`frames[0].messages`, assistant `tool_use` blocks). They document how the Tier A controls and the
v1/v2 selection test were computed. **They are not runnable scripts**: cells depend on kernel state
from earlier cells, on files in that session's private workspace, and occasionally on the platform's
artifact store. Only absolute paths were replaced by placeholders (`<OLD_WORKSPACE>`,
`<OLD_ARTIFACT_STORE>/`, `<WHEELER_LAB_DRIVE>`, `<HOME>`); a credential scan found nothing to remove.

Cells that only read the reviewer letter, browsed directories or checked tool availability are kept so
that the record is complete. Planning documents written with a file editor and save/plan bookkeeping
calls are not included (no code).

| File | Messages | Analysis | Outputs used downstream |
|---|---|---|---|
| [`00_session_setup_and_esl_translate_replica.md`](00_session_setup_and_esl_translate_replica.md) | 3–143 | Session setup: reading the reviewer letter, repo review, toolchain discovery, construction and validation of the pure-Python esl-translate replica, copying of inputs | `esl_translate_replica.py` (now `utils/`) |
| [`01_baseline_figure2_reproduction.md`](01_baseline_figure2_reproduction.md) | 145–205 | Reproduction of the published Figure 2A-2C (baseline before any control); defines helper functions reused by later steps | baseline Fig 2 numbers |
| [`02_tierA_gc_audit.md`](02_tierA_gc_audit.md) | 210–229 | Tier A GC audit: per-exon GC of DCR vs non-DCR exons (human) and mouse orthologs vs random control | `gc_audit_summary.csv`, `gc_audit_tests.csv` (Fig S8 A-B) |
| [`03_gc_length_matched_controls.md`](03_gc_length_matched_controls.md) | 231–247 | GC+length-matched random-region controls and GC-matched Figure 2B curves | `gc_matched_fig2b.csv` (Fig S8 C) |
| [`04_dcr_self_shuffle.md`](04_dcr_self_shuffle.md) | 249–265 | DCR self-shuffle null: codon-shuffle of each DCR ortholog's own sequence | `self_shuffle_control.csv` (Fig S8 D) |
| [`05_latent_potential.md`](05_latent_potential.md) | 267–278 | Latent dual-coding potential in GC+length-matched non-DCR human exons | latent-potential fractions (Fig 2B(i) text) |
| [`06_phyloP_gc_control.md`](06_phyloP_gc_control.md) | 280–301 | phyloP conservation with GC control (matched resampling + linear model) | `phyloP_gc_control.csv` (Fig S8 E) |
| [`07_gnas_dcr_gt500.md`](07_gnas_dcr_gt500.md) | 303–330 | GNAS: DCRs > 500 nt and the mouse dual-coding ORFs of DCR 201 | `gnas_dcr_gt500.csv` (Fig S7) |
| [`08_tierA_results_report.md`](08_tierA_results_report.md) | 332–342 | Assembly of the Tier A results document (no new analysis) | Tier A results document |
| [`09_definition_reconciliation.md`](09_definition_reconciliation.md) | 344–392 | Reconciliation of the 82% vs 48% dual-coding definitions (strict vs sub-interval rule) | `definition_reconciliation.csv` (Fig 2B(i)) |
| [`10_xprime_conservation_filter.md`](10_xprime_conservation_filter.md) | 393–449 | Conservation as a specificity filter, X vs X' (pyliftover, then UCSC liftOver re-lift); builds dcr_regions.fa, N/D human and mouse FASTAs | `xprime_summary.csv` (Fig 2B(i)-(ii)); `dcr_regions.fa`, `D_ucsc.fa`, `N_ucsc.fa`, `nondcr_pool.fa` (selection inputs) |
| [`11_annotation_depth_check.md`](11_annotation_depth_check.md) | 451–481 | Annotation-depth (transcript count) check of DCR discovery and the X' gap | `annotation_depth_check.csv`, `dcr_discovery_vs_annotation.csv` (Results text) |
| [`12_selection_test_v1_v2.md`](12_selection_test_v1_v2.md) | 484–625 | Alternate-frame selection test v1 (GENCODE frames) and v2 (protein-derived DCR frames): CDS tables, intersections, OLGenie scoring, ascertainment checks | `N_human.fa`, `cds_exons.tsv`, `*_cds*.ov`, `dcr_exclusion_fates.csv`, `dcr_protein_frames.csv` (selection inputs); v2 `selection_omega.csv` |

The selection-test intermediates built here are verified by md5 against
`../selection/selection_inputs_manifest.tsv` when they are staged with
`../selection/stage_selection_inputs.py`.
