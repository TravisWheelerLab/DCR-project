# config.sh -- paths for Code/05_Revision. Sourced by run_revision.sh; every value can be overridden
# from the environment (export DOOVER=...; ./run_revision.sh) or by pointing CONFIG=/path/to/my_config.sh.
# All scripts also accept the same paths as CLI arguments (see --help of each script).

# Executed-pipeline directory (the source of truth for every DCR call; the layout written by
# Code/01-04: output/01_Description, output/02_Conservation, output/New_TableS1.tsv, data/)
: "${DOOVER:=$PWD/DOOVER}"

# Reference data. DOOVER/data holds symlinks to these files; point REFDIR at the real copies if the
# links are dangling on your machine.
: "${REFDIR:=$DOOVER/data}"
: "${HG38:=$REFDIR/hg38.fa}"                       # + hg38.fa.fai
: "${PHYLOP_BW:=$REFDIR/hg38.phyloP100way.bw}"     # 9.9 GB; COPY IT TO A LOCAL DISK (see README)

# Prior-session intermediates for the selection test (searched recursively; md5-checked).
# Typically the folder holding the exported artifacts of the July-August 2026 revision session.
: "${PRIOR_EXPORT:=}"

# Where all revision outputs go (created). Nothing is written inside DOOVER.
: "${WORK:=$PWD/revision_work}"

# Interpreters / tools
: "${PY:=python3}"           # python 3.11 with the packages in environment.yml
: "${RSCRIPT:=Rscript}"      # R 4.3 with the packages in environment.yml
: "${NMD_CHUNKS:=8}"          # parallel chunks for splicing_vs_AS/nmd_allgenes.py

export DOOVER REFDIR HG38 PHYLOP_BW PRIOR_EXPORT WORK PY RSCRIPT NMD_CHUNKS
