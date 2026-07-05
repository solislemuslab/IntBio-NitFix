#!/usr/bin/env bash
set -euo pipefail

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2}"
SCRIPT="$OUT/Rscripts_v2/12_make_nifH_clean5_tree_ready_v2.py"
LOGDIR="$OUT/nifH_clean5_tree_ready_v2"

mkdir -p "$LOGDIR"

echo "Step 12 V2: make tree-ready nifH clean5 FASTA files"
echo "OUT: $OUT"
echo "Clean refs: ${CLEAN_REFS:-ref52,ref56,ref60,ref62,ref63}"
echo

python3 "$SCRIPT"

echo
echo "Counts:"
grep -c '^>' "$OUT/nifH_clean5_tree_ready_v2/nifH_clean5_strict_single_dominant.fasta"
grep -c '^>' "$OUT/nifH_clean5_tree_ready_v2/nifH_clean5_iupac_all_pass.fasta"
