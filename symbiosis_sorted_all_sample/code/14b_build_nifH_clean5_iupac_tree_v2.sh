#!/usr/bin/env bash
set -euo pipefail

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2}"

echo "Build mixed-aware IUPAC nifH clean5 tree"
echo "OUT: $OUT"
echo

mkdir -p "$OUT/nifH_clean5_mafft_alignment_v2"
mkdir -p "$OUT/nifH_clean5_iqtree_v2"

INPUT_FASTA="$OUT/nifH_clean5_tree_ready_v2/nifH_clean5_iupac_all_pass.fasta" \
LABEL="iupac_all_pass" \
bash "$OUT/Rscripts_v2/13_align_nifH_clean5_with_mafft_v2.sh" \
  | tee "$OUT/nifH_clean5_mafft_alignment_v2/mafft_iupac_all_pass_run.log"

ALIGNMENT="$OUT/nifH_clean5_mafft_alignment_v2/nifH_clean5_iupac_all_pass.mafft.fasta" \
PREFIX="$OUT/nifH_clean5_iqtree_v2/nifH_clean5_iupac_all_pass" \
bash "$OUT/Rscripts_v2/14_build_nifH_clean5_tree_iqtree_v2.sh" \
  | tee "$OUT/nifH_clean5_iqtree_v2/iqtree_iupac_all_pass_run.log"
