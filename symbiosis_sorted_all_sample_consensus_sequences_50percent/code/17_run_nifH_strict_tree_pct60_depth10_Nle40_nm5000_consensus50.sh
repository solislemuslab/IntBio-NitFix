#!/usr/bin/env bash
set -euo pipefail

# Build a strict single-dominant nifH tree from the consensus50 reference
# using a relaxed recovery threshold and extended UFBoot iterations.
#
# Selection thresholds:
#   percent_covered_from_pileup >= 60
#   mean_depth_from_pileup      >= 10
#   N_percent                   <= 40
#
# Tree settings:
#   IQ-TREE ModelFinder (-m MFP)
#   1000 ultrafast bootstraps (-B 1000)
#   1000 SH-aLRT replicates (-alrt 1000)
#   maximum UFBoot iterations (-nm 5000)
#
# This script writes to nm5000-specific folders, so it does not overwrite
# previous pct80/depth10/Nle20 or earlier pct60/depth10/Nle40 results.

module load SolisLemus-BioPhylo/2026.04.20 || true

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"

MIN_PCT="${MIN_PCT:-60}"
MIN_MEAN_DEPTH="${MIN_MEAN_DEPTH:-10}"
MAX_N_PERCENT="${MAX_N_PERCENT:-40}"
JOBS="${JOBS:-2}"
FORCE="${FORCE:-0}"

RUN_DIR="$OUT/nifH_consensus_tree_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}_nm5000_v2"
ALIGN_DIR="$OUT/nifH_consensus_mafft_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}_nm5000_v2"
TREE_DIR="$OUT/nifH_consensus_iqtree_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}_nm5000_v2"

STRICT_FASTA="$RUN_DIR/nifH_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.fasta"
STRICT_ALIGNED="$ALIGN_DIR/nifH_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.mafft.fasta"
PREFIX="$TREE_DIR/nifH_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}_nm5000"

mkdir -p "$RUN_DIR" "$ALIGN_DIR" "$TREE_DIR"

if [[ -s "$PREFIX.treefile" && "$FORCE" != "1" ]]; then
  echo "ERROR: output tree already exists: $PREFIX.treefile" >&2
  echo "Set FORCE=1 only if you intentionally want to rerun and overwrite it." >&2
  exit 1
fi

echo "Step 17 nm5000: strict nifH tree from consensus50 reference"
echo "OUT:              $OUT"
echo "Consensus folder: $RUN_DIR"
echo "Alignment folder: $ALIGN_DIR"
echo "Tree folder:      $TREE_DIR"
echo "MIN_PCT:          $MIN_PCT"
echo "MIN_MEAN_DEPTH:   $MIN_MEAN_DEPTH"
echo "MAX_N_PERCENT:    $MAX_N_PERCENT"
echo "JOBS:             $JOBS"
echo "IQ-TREE -nm:      5000"
echo

echo "1. Extracting strict single-dominant nifH consensus sequences..."
MIN_PCT="$MIN_PCT" \
MIN_MEAN_DEPTH="$MIN_MEAN_DEPTH" \
MAX_N_PERCENT="$MAX_N_PERCENT" \
JOBS="$JOBS" \
RUN_DIR="$RUN_DIR" \
bash "$OUT/Rscripts_v2/13_run_nifH_consensus_consensus50.sh" \
  | tee "$RUN_DIR/nifH_consensus_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}_run.log"

echo
echo "2. Checking strict FASTA..."
if [[ ! -s "$STRICT_FASTA" ]]; then
  echo "ERROR: strict FASTA was not created or is empty: $STRICT_FASTA" >&2
  exit 1
fi
STRICT_COUNT=$(grep -c '^>' "$STRICT_FASTA")
echo "Strict sequence count: $STRICT_COUNT"

echo
echo "3. Aligning strict FASTA with MAFFT..."
if ! command -v mafft >/dev/null 2>&1; then
  echo "ERROR: mafft is not available in PATH." >&2
  echo "Load the module/environment that contains MAFFT, then run again." >&2
  exit 1
fi

mafft --auto "$STRICT_FASTA" \
  > "$STRICT_ALIGNED" \
  2> "$ALIGN_DIR/nifH_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.mafft.log"

ALIGNED_COUNT=$(grep -c '^>' "$STRICT_ALIGNED")
echo "Aligned sequence count: $ALIGNED_COUNT"

echo
echo "Alignment length distribution:"
awk '
  /^>/ {
    if (seq != "") print length(seq)
    seq=""
    next
  }
  { seq=seq $0 }
  END { if (seq != "") print length(seq) }
' "$STRICT_ALIGNED" | sort | uniq -c | tee "$ALIGN_DIR/alignment_length_distribution.txt"

echo
echo "4. Building strict tree with IQ-TREE..."
if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE="iqtree"
else
  echo "ERROR: neither iqtree2 nor iqtree is available in PATH." >&2
  exit 1
fi
echo "Using: $IQTREE"

REDO_ARG=()
if [[ "$FORCE" == "1" ]]; then
  REDO_ARG=(-redo)
fi

"$IQTREE" \
  -s "$STRICT_ALIGNED" \
  -st DNA \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -nm 5000 \
  -T AUTO \
  --prefix "$PREFIX" \
  "${REDO_ARG[@]}" \
  | tee "$TREE_DIR/nifH_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}_nm5000_iqtree_run.log"

echo
echo "Tree finished. Checking key outputs..."
ls -lh "$PREFIX.treefile" "$PREFIX.contree" "$PREFIX.iqtree" "$PREFIX.log"

echo
echo "Best-fit model:"
grep -i "Best-fit model" "$PREFIX.iqtree" || true

echo
echo "Convergence/warning messages:"
grep -i "correlation coefficient\|converge\|WARNING\|ERROR" "$PREFIX.log" || true

echo
echo "Tree tip count:"
python3 - <<PY
from pathlib import Path

tree = Path("$PREFIX.treefile").read_text().strip()
labels = []
token = []
in_quote = False
quote_char = ""

def flush():
    global token
    if token:
        s = "".join(token).strip()
        if s:
            labels.append(s)
        token = []

i = 0
while i < len(tree):
    ch = tree[i]
    if in_quote:
        if ch == quote_char:
            in_quote = False
        else:
            token.append(ch)
    else:
        if ch in "'\"":
            in_quote = True
            quote_char = ch
        elif ch in ",();":
            flush()
        elif ch == ":":
            flush()
            i += 1
            while i < len(tree) and tree[i] not in ",();":
                i += 1
            continue
        else:
            token.append(ch)
    i += 1
flush()

tips = [x for x in labels if "|nifH|consensus50|" in x]
print(len(tips))
PY

echo
echo "Relaxed strict nifH nm5000 tree complete."
echo "QC table:       $RUN_DIR/nifH_consensus_qc.tsv"
echo "Strict FASTA:   $STRICT_FASTA"
echo "Alignment:      $STRICT_ALIGNED"
echo "Main tree:      $PREFIX.treefile"
echo "Consensus tree: $PREFIX.contree"
echo "IQ-TREE report: $PREFIX.iqtree"
echo "Run log:        $PREFIX.log"
