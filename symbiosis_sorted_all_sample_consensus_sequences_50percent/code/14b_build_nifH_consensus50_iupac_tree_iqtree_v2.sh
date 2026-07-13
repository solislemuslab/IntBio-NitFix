#!/usr/bin/env bash
set -euo pipefail

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
ALIGNDIR="$OUT/nifH_consensus_mafft_v2"
TREEDIR="$OUT/nifH_consensus_iqtree_v2"

ALIGNMENT="${ALIGNMENT:-$ALIGNDIR/nifH_consensus50_iupac_all_pass.mafft.fasta}"
PREFIX="${PREFIX:-$TREEDIR/nifH_consensus50_iupac_all_pass}"
THREADS="${THREADS:-AUTO}"
BOOTSTRAPS="${BOOTSTRAPS:-1000}"
ALRT="${ALRT:-1000}"
MODEL="${MODEL:-MFP}"
SEQTYPE="${SEQTYPE:-DNA}"

mkdir -p "$TREEDIR"

echo "Step 14b V2 consensus50: build mixed/IUPAC nifH tree with IQ-TREE"
echo "Alignment:  $ALIGNMENT"
echo "Prefix:     $PREFIX"
echo "Seq type:   $SEQTYPE"
echo "Model:      $MODEL"
echo "Bootstrap:  $BOOTSTRAPS ultrafast bootstraps"
echo "SH-aLRT:    $ALRT replicates"
echo "Threads:    $THREADS"
echo

if [[ ! -s "$ALIGNMENT" ]]; then
  echo "ERROR: alignment does not exist or is empty: $ALIGNMENT"
  exit 1
fi

if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE="iqtree"
else
  echo "ERROR: neither iqtree2 nor iqtree is available in PATH."
  echo "Load the IQ-TREE module/environment, then run again."
  exit 1
fi

echo "Using: $IQTREE"
echo

"$IQTREE" \
  -s "$ALIGNMENT" \
  -st "$SEQTYPE" \
  -m "$MODEL" \
  -B "$BOOTSTRAPS" \
  -alrt "$ALRT" \
  -T "$THREADS" \
  --prefix "$PREFIX"

echo
echo "IQ-TREE complete."
echo "Main tree:       ${PREFIX}.treefile"
echo "Model/log file:  ${PREFIX}.iqtree"
echo "Run log:         ${PREFIX}.log"
echo "Bootstrap tree:  ${PREFIX}.ufboot"
