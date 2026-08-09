#!/usr/bin/env bash
set -euo pipefail

# nifK: align strict and mixed-IUPAC consensus FASTA files with MAFFT.

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
GENE="nifK"
MIN_PCT="${MIN_PCT:-80}"
MIN_MEAN_DEPTH="${MIN_MEAN_DEPTH:-10}"
MAX_N_PERCENT="${MAX_N_PERCENT:-20}"
THREADS="${THREADS:-8}"
RUN_DIR="${OUT}/gene_trees/nifK/01_consensus"
ALIGN_DIR="${OUT}/gene_trees/nifK/02_alignment"
mkdir -p "$ALIGN_DIR"

if command -v module >/dev/null 2>&1; then
  module load SolisLemus-BioPhylo/2026.04.20 || true
fi

for SET in strict_single_dominant iupac_all_pass; do
  if [[ "$SET" == "strict_single_dominant" ]]; then
    IN="$RUN_DIR/nifK_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.fasta"
    OUTFA="$ALIGN_DIR/nifK_consensus50_strict_single_dominant.mafft.fasta"
  else
    IN="$RUN_DIR/nifK_consensus50_iupac_all_pass_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.fasta"
    OUTFA="$ALIGN_DIR/nifK_consensus50_iupac_all_pass.mafft.fasta"
  fi

  if [[ ! -s "$IN" ]]; then
    echo "ERROR: input FASTA missing or empty: $IN" >&2
    exit 1
  fi

  echo "Aligning nifK $SET ($(grep -c '^>' "$IN") sequences)"
  mafft --thread "$THREADS" --auto "$IN" > "$OUTFA" 2> "$OUTFA.log"
  awk '/^>/ {if (seq!="") print length(seq); seq=""; next} {seq=seq $0} END {if (seq!="") print length(seq)}' "$OUTFA"     | sort | uniq -c > "$OUTFA.lengths.txt"
  cat "$OUTFA.lengths.txt"
done
