#!/usr/bin/env bash
# run_revision.sh -- run the 2026 revision analyses (Code/05_Revision) in dependency order.
#
#   ./run_revision.sh                 # all steps
#   ./run_revision.sh selection figures   # selected steps (names below), in the order given
#   CONFIG=my_config.sh ./run_revision.sh ...
#
# Steps (dependency order; approximate runtime on a 2024 laptop, local reference files):
#   stage_selection  md5-checked staging of the 11 selection intermediates from $PRIOR_EXPORT   <1 min
#   selection        alternate-frame selection test v3                         -> Fig 2B(iii)     ~1 min
#   ascertainment    GABBR1/PTBP2 showcase ascertainment table                 -> Fig S6 text     <1 min
#   splice_sites     splice-site phyloP, DCR vs control                         -> Fig S9A        ~10 s (local bigWig)
#   gtex_fetch       GTEx v8 median TPM per gene (public API, network)          -> Fig S9B-C      ~3 min
#   expression_tau   expression class / tau / brain by mouse conservation       -> Fig S9B-C      <1 min
#   splicing_vs_AS   DCR vs ordinary alternative splicing (A2), incl. NMD rerun -> Fig S10        ~5 min
#   figure_data      GC tables + GNAS/GABBR1 gene models from DOOVER            -> Fig S6-S8      ~10 s
#   latent_check     exon-level check of Fig 2B(i) fractions (optional)          -> (check only)   ~20 s
#   figure_inputs    assemble $WORK/figure_inputs from results/ + fresh outputs                  <1 s
#   figures          render Figure 2B and Supplementary Figures S6-S10                           ~1 min
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${CONFIG:-$HERE/config.sh}"
mkdir -p "$WORK"
log() { printf '\n=== [%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

step_stage_selection() {
  if [[ -z "$PRIOR_EXPORT" ]]; then
    echo "PRIOR_EXPORT not set: skipping staging (selection expects $WORK/selection_inputs)"; return 0
  fi
  "$PY" "$HERE/selection/stage_selection_inputs.py" --search "$PRIOR_EXPORT" --dest "$WORK/selection_inputs"
}
step_selection()      { "$PY" "$HERE/selection/selection_v3.py" --inputs-dir "$WORK/selection_inputs" --outdir "$WORK/selection"; }
step_ascertainment()  { "$PY" "$HERE/controls/gabbr1_ascertainment.py" --outdir "$WORK/controls"; }
step_splice_sites()   { "$PY" "$HERE/controls/splice_site_phylop.py" --phylop "$PHYLOP_BW" --outdir "$WORK/controls"; }
step_gtex_fetch()     { "$PY" "$HERE/controls/gtex_fetch_medians.py" --out "$WORK/controls/gtex_v8_gene_medians.json"; }
step_expression_tau() { "$PY" "$HERE/controls/expression_tau.py" --expr "$WORK/controls/gtex_v8_gene_medians.json" --outdir "$WORK/controls"; }

step_splicing_vs_AS() {
  local A="$HERE/splicing_vs_AS" O="$WORK/splicing_vs_AS"
  mkdir -p "$O/nmd_chunks"
  "$PY" "$A/A2_build_events.py" --hg38 "$HG38" --outdir "$O"
  local pids=() k
  for ((k = 0; k < NMD_CHUNKS; k++)); do
    "$PY" "$A/nmd_allgenes.py" "$k" "$NMD_CHUNKS" --outdir "$O" > "$O/nmd_chunks/chunk.$k.log" 2>&1 & pids+=($!)
  done
  for k in "${pids[@]}"; do wait "$k"; done
  : > "$O/NMDtable_allgenes.tsv"
  for ((k = 0; k < NMD_CHUNKS; k++)); do cat "$O/nmd_chunks/NMDtable_all.$k.tsv" >> "$O/NMDtable_allgenes.tsv"; done
  "$PY" "$A/A2_event_tables.py" --hg38 "$HG38" --phylop "$PHYLOP_BW" --outdir "$O"
  "$PY" "$A/A2_stats.py" --outdir "$O"
}

step_figure_data() {
  "$PY" "$HERE/figures/make_gc_tables.py"   --hg38 "$HG38" --datadir "$WORK/figure_data"
  "$PY" "$HERE/figures/prep_gene_models.py" --hg38 "$HG38" --datadir "$WORK/figure_data"
}
step_latent_check() { "$PY" "$HERE/figures/verify_latent_exon_fraction.py" --hg38 "$HG38" --datadir "$WORK/figure_data"; }

step_figure_inputs() {
  # archived summary tables first, then overwrite with anything recomputed in $WORK
  local F="$WORK/figure_inputs" f
  mkdir -p "$F"
  cp "$HERE"/results/*.csv "$F"/
  for f in selection/selection_omega_v3.csv controls/splice_site_phylop.tsv controls/splice_site_conservation.csv \
           controls/A4_dcr_by_expression.csv controls/A4_tau_by_expression.csv \
           controls/A4_gene_expression_tau.csv controls/A4_brain_by_conservation.csv; do
    if [[ -f "$WORK/$f" ]]; then cp "$WORK/$f" "$F"/; else echo "note: $WORK/$f absent, using results/ copy if any"; fi
  done
  for f in splice_site_phylop.tsv A4_gene_expression_tau.csv; do
    [[ -f "$F/$f" ]] || echo "WARNING: $f missing (run splice_sites / gtex_fetch + expression_tau) -> Fig S9 will fail"
  done
}

step_figures() {
  local G="$HERE/figures" IN="$WORK/figure_inputs" D="$WORK/figure_data" OUT="$WORK/figures"
  mkdir -p "$OUT"
  "$RSCRIPT" "$G/fig2B_potential_vs_selection.R"   "$IN" "$OUT"          # Figure 2B(i)-(iii)
  "$RSCRIPT" "$G/supp_GABBR1_gene_model.R"         "$D" "$OUT"           # Fig S6
  "$RSCRIPT" "$G/supp_GNAS_gene_model.R"           "$D" "$IN" "$OUT"     # Fig S7
  "$RSCRIPT" "$G/supp_gc_content_controls.R"       "$IN" "$D" "$OUT"     # Fig S8
  "$RSCRIPT" "$G/supp_splice_expression_tissue.R"  "$IN" "$OUT"          # Fig S9
  if [[ -f "$WORK/splicing_vs_AS/A2_events_annotated.csv" ]]; then         # Fig S10
    "$RSCRIPT" "$HERE/splicing_vs_AS/A2_figure.R" "$WORK/splicing_vs_AS" "$OUT"
  else
    echo "note: $WORK/splicing_vs_AS/A2_events_annotated.csv absent -> run step splicing_vs_AS for Fig S10"
  fi
}

ALL=(stage_selection selection ascertainment splice_sites gtex_fetch expression_tau splicing_vs_AS
     figure_data figure_inputs figures)
STEPS=("$@"); [[ ${#STEPS[@]} -eq 0 ]] && STEPS=("${ALL[@]}")
for s in "${STEPS[@]}"; do
  declare -F "step_$s" > /dev/null || { echo "unknown step: $s (known: ${ALL[*]} latent_check)"; exit 2; }
  log "$s"; SECONDS=0; "step_$s"; echo "--- $s done in ${SECONDS}s"
done
