#!/usr/bin/env bash
set -euo pipefail

# Run one selected pipeline step for multiple per-gene folders.
# Example:
#   STEP=02 JOBS=2 bash run_selected_step_for_genes.sh
#   STEP=04 JOBS=2 bash run_selected_step_for_genes.sh

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
GENES="${GENES:-nifH,nifD,nifK,nifJ,nodL,nolG,nolF,noeA,noeB,nodX}"
STEP="${STEP:-}"
JOBS="${JOBS:-2}"
SCRIPT_ROOT="$OUT/Rscripts_v2/gene_tree_pipelines"

if [[ -z "$STEP" ]]; then
  echo "ERROR: set STEP to one of: 02, 03, 04, 05, 06, 07" >&2
  exit 1
fi

script_for_gene () {
  gene="$1"
  case "$STEP" in
    02) echo "$SCRIPT_ROOT/$gene/02_run_${gene}_consensus_consensus50.sh" ;;
    03) echo "$SCRIPT_ROOT/$gene/03_align_${gene}_consensus50_with_mafft.sh" ;;
    04) echo "$SCRIPT_ROOT/$gene/04_build_${gene}_strict_tree_nm5000_consensus50.sh" ;;
    05) echo "$SCRIPT_ROOT/$gene/05_build_${gene}_iupac_tree_nm5000_consensus50.sh" ;;
    06) echo "$SCRIPT_ROOT/$gene/06_blast_${gene}_consensus_to_taxon_reference.sh" ;;
    07) echo "$SCRIPT_ROOT/$gene/07_check_${gene}_tree_outputs.sh" ;;
    *) echo "ERROR: unknown STEP=$STEP" >&2; return 1 ;;
  esac
}
export -f script_for_gene
export STEP SCRIPT_ROOT OUT

IFS=',' read -r -a GENE_ARRAY <<< "$GENES"
printf "%s
" "${GENE_ARRAY[@]}" | parallel -j "$JOBS" 's=$(script_for_gene {}); echo "==== {} STEP '$STEP' ===="; bash "$s"'
