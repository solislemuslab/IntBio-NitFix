#!/usr/bin/env bash
set -euo pipefail

# Build a strict single-dominant nifH tree using the consensus50 reference
# with a relaxed sample-selection threshold:
#   percent covered >= 60
#   mean depth >= 10
#   N percent <= 40
#
# This script intentionally writes to new pct60/depth10/Nle40-specific folders
# so it does not overwrite the previous pct80/depth10/Nle20 or
# pct60/depth10/Nle20 consensus, alignment, or tree.

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"

MIN_PCT="${MIN_PCT:-60}"
MIN_MEAN_DEPTH="${MIN_MEAN_DEPTH:-10}"
MAX_N_PERCENT="${MAX_N_PERCENT:-40}"
JOBS="${JOBS:-2}"

RUN60="$OUT/nifH_consensus_tree_pct60_depth10_Nle40_v2"
ALIGN60="$OUT/nifH_consensus_mafft_pct60_depth10_Nle40_v2"
TREE60="$OUT/nifH_consensus_iqtree_pct60_depth10_Nle40_v2"

STRICT_FASTA="$RUN60/nifH_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.fasta"
STRICT_ALIGNED="$ALIGN60/nifH_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.mafft.fasta"
TREE_PREFIX="$TREE60/nifH_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}"

mkdir -p "$RUN60" "$ALIGN60" "$TREE60"

echo "Step 17: strict nifH tree from consensus50 reference at pct${MIN_PCT}/depth${MIN_MEAN_DEPTH}"
echo "OUT:              $OUT"
echo "Consensus folder: $RUN60"
echo "Alignment folder: $ALIGN60"
echo "Tree folder:      $TREE60"
echo "MIN_PCT:          $MIN_PCT"
echo "MIN_MEAN_DEPTH:   $MIN_MEAN_DEPTH"
echo "MAX_N_PERCENT:    $MAX_N_PERCENT"
echo "JOBS:             $JOBS"
echo

echo "1. Extracting nifH consensus sequences..."
MIN_PCT="$MIN_PCT" \
MIN_MEAN_DEPTH="$MIN_MEAN_DEPTH" \
MAX_N_PERCENT="$MAX_N_PERCENT" \
JOBS="$JOBS" \
RUN_DIR="$RUN60" \
bash "$OUT/Rscripts_v2/13_run_nifH_consensus_consensus50.sh" \
  | tee "$RUN60/nifH_consensus_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_run.log"

echo
echo "2. Checking strict FASTA..."
if [[ ! -s "$STRICT_FASTA" ]]; then
  echo "ERROR: strict FASTA was not created or is empty: $STRICT_FASTA"
  exit 1
fi
grep -c '^>' "$STRICT_FASTA" | awk '{print "Strict sequence count:", $1}'

echo
echo "3. Aligning strict FASTA with MAFFT..."
if ! command -v mafft >/dev/null 2>&1; then
  echo "ERROR: mafft is not available in PATH."
  echo "Load the module/environment that contains MAFFT, then run again."
  exit 1
fi

mafft --auto "$STRICT_FASTA" \
  > "$STRICT_ALIGNED" \
  2> "$ALIGN60/nifH_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.mafft.log"

grep -c '^>' "$STRICT_ALIGNED" | awk '{print "Aligned sequence count:", $1}'

echo
echo "4. Building strict tree with IQ-TREE..."
ALIGNMENT="$STRICT_ALIGNED" \
PREFIX="$TREE_PREFIX" \
MODEL=MFP \
BOOTSTRAPS=1000 \
ALRT=1000 \
THREADS=AUTO \
bash "$OUT/Rscripts_v2/14_build_nifH_consensus50_strict_tree_iqtree_v2.sh" \
  | tee "$TREE60/nifH_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}_iqtree_run.log"

echo
echo "Done."
echo "QC table:       $RUN60/nifH_consensus_qc.tsv"
echo "Strict FASTA:   $STRICT_FASTA"
echo "Alignment:      $STRICT_ALIGNED"
echo "Tree file:      ${TREE_PREFIX}.treefile"
echo "IQ-TREE report: ${TREE_PREFIX}.iqtree"
echo "IQ-TREE log:    ${TREE_PREFIX}.log"
