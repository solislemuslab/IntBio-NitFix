#!/usr/bin/env bash
set -euo pipefail

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
GENE="nifH"
ROOT="$OUT/gene_trees/nifH"
echo "==== nifH consensus ===="
if [[ -s "$ROOT/01_consensus/nifH_consensus_qc.tsv" ]]; then
  tail -n +2 "$ROOT/01_consensus/nifH_consensus_qc.tsv" | cut -f18 | sort | uniq -c
fi
for f in "$ROOT"/01_consensus/*.fasta "$ROOT"/02_alignment/*.fasta "$ROOT"/03_tree_strict_nm5000/*.treefile "$ROOT"/04_tree_iupac_nm5000/*.treefile "$ROOT"/05_blast_taxon_annotation/*best_reference_hit_with_taxon.tsv; do
  [[ -e "$f" ]] || continue
  echo "---- $f"
  ls -lh "$f"
  case "$f" in
    *.fasta) grep -c '^>' "$f" ;;
    *.treefile) grep -o "[^,():;][^,():;]*" "$f" | grep -v "^[0-9.eE+-]*$" | wc -l ;;
    *.tsv) tail -n +2 "$f" | wc -l ;;
  esac
done

echo "==== IQ-TREE warnings ===="
grep -ih "Best-fit model\|WARNING\|correlation coefficient\|did not converge\|ERROR" "$ROOT"/03_tree_strict_nm5000/*.log "$ROOT"/04_tree_iupac_nm5000/*.log 2>/dev/null || true
