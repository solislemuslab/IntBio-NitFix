#!/usr/bin/env bash
set -euo pipefail

# noeB: build strict single-dominant IQ-TREE with ModelFinder, UFBoot, SH-aLRT, nm5000.
# Uses -safe only for nolF because nolF previously triggered IQ-TREE numerical underflow.

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
GENE="noeB"
ALIGN="${OUT}/gene_trees/noeB/02_alignment/noeB_consensus50_strict_single_dominant.mafft.fasta"
TREE_DIR="${OUT}/gene_trees/noeB/03_tree_strict_nm5000"
PREFIX="$TREE_DIR/noeB_consensus50_strict_single_dominant_nm5000"
FORCE="${FORCE:-0}"
mkdir -p "$TREE_DIR"

if [[ ! -s "$ALIGN" ]]; then echo "ERROR: alignment missing: $ALIGN" >&2; exit 1; fi
if [[ -s "$PREFIX.treefile" && "$FORCE" != "1" ]]; then echo "ERROR: output tree exists: $PREFIX.treefile; set FORCE=1 to overwrite" >&2; exit 1; fi

if command -v module >/dev/null 2>&1; then module load SolisLemus-BioPhylo/2026.04.20 || true; fi
if command -v iqtree2 >/dev/null 2>&1; then IQTREE=iqtree2; elif command -v iqtree >/dev/null 2>&1; then IQTREE=iqtree; else echo "ERROR: iqtree not found" >&2; exit 1; fi

"$IQTREE"   -s "$ALIGN"   -st DNA   -m MFP   -B 1000   -alrt 1000   -nm 5000   -T AUTO   --prefix "$PREFIX"   | tee "$PREFIX.run.log"

echo "Tree tip count:"
grep -o "[^,():;][^,():;]*" "$PREFIX.treefile" | grep -v "^[0-9.eE+-]*$" | wc -l || true
echo "Strict tree complete: $PREFIX.treefile"
