#!/usr/bin/env bash
set -euo pipefail

# Step 20: Build a nifH reference-taxon tree
# Purpose:
#   Build a phylogenetic tree whose tips are known/reference bacterial nifH taxa,
#   not sample-derived consensus sequences.
#
# Input:
#   nifH.fa.aln = aligned reference nifH sequences, expected headers like:
#   >nifH|Mesorhizobium_sp._AaZ16|XIK06011.1
#
# Output:
#   A cleaned-label reference alignment and an IQ-TREE maximum-likelihood tree.

module load SolisLemus-BioPhylo/2026.04.20 || true
module load blast+/2.17.0 || true

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
REF_ALN="${REF_ALN:-$OUT/reference_nifH_tree/nifH.fa.aln}"
TREE_DIR="${TREE_DIR:-$OUT/reference_nifH_tree/taxon_tree_nm5000}"
FORCE="${FORCE:-0}"

CLEAN_ALN="$TREE_DIR/nifH_reference_746_taxa_labels.aln.fasta"
TREE_PREFIX="$TREE_DIR/nifH_reference_746_taxa_tree_nm5000"

mkdir -p "$TREE_DIR"

if [[ ! -s "$REF_ALN" ]]; then
  echo "ERROR: reference alignment not found or empty: $REF_ALN" >&2
  exit 1
fi

if [[ -s "${TREE_PREFIX}.treefile" && "$FORCE" != "1" ]]; then
  echo "ERROR: tree already exists: ${TREE_PREFIX}.treefile" >&2
  echo "Set FORCE=1 if you intentionally want to overwrite/re-run this output." >&2
  exit 1
fi

if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE="iqtree"
else
  echo "ERROR: neither iqtree2 nor iqtree is available in PATH." >&2
  exit 1
fi

echo "Building reference-taxon nifH tree"
echo "OUT:        $OUT"
echo "REF_ALN:    $REF_ALN"
echo "TREE_DIR:   $TREE_DIR"
echo "IQ-TREE:    $IQTREE"
echo

echo "Checking input reference alignment..."
seq_count=$(grep -c '^>' "$REF_ALN")
echo "Reference sequences: $seq_count"

awk '
  /^>/ {
    if (seq != "") print length(seq)
    seq=""
    next
  }
  { seq=seq $0 }
  END { if (seq != "") print length(seq) }
' "$REF_ALN" | sort | uniq -c > "$TREE_DIR/reference_alignment_lengths.txt"

echo "Alignment length distribution:"
cat "$TREE_DIR/reference_alignment_lengths.txt"
echo

echo "Creating cleaned taxon labels for tree tips..."
awk '
  /^>/ {
    header=substr($0,2)
    split(header,a,"|")
    taxon=a[2]
    accession=a[3]

    if (taxon == "" || accession == "") {
      label=header
    } else {
      label=taxon "|" accession
    }

    gsub(/[^A-Za-z0-9_.|+-]/, "_", label)
    print ">" label
    next
  }
  { print }
' "$REF_ALN" > "$CLEAN_ALN"

echo "Cleaned alignment: $CLEAN_ALN"
echo "Cleaned sequence count: $(grep -c '^>' "$CLEAN_ALN")"
echo

echo "Running IQ-TREE reference-taxon tree..."
"$IQTREE" \
  -s "$CLEAN_ALN" \
  -st DNA \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -nm 5000 \
  -T AUTO \
  --prefix "$TREE_PREFIX" \
  ${FORCE:+-redo} \
  | tee "$TREE_DIR/iqtree_reference_taxon_tree_run.log"

echo
echo "Checking IQ-TREE warnings/convergence messages..."
grep -i "converge\|correlation coefficient\|WARNING\|ERROR" "${TREE_PREFIX}.log" || true

echo
echo "Reference-taxon nifH tree complete."
echo "Main tree:       ${TREE_PREFIX}.treefile"
echo "Consensus tree:  ${TREE_PREFIX}.contree"
echo "IQ-TREE report:  ${TREE_PREFIX}.iqtree"
echo "Run log:         ${TREE_PREFIX}.log"
