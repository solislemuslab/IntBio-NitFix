#!/usr/bin/env bash
set -euo pipefail

# Step 24: Build the mixed-IUPAC nifH sample tree with extended UFBoot iterations.
#
# Purpose:
#   Re-run the mixed-IUPAC nifH tree with -nm 5000 because the earlier mixed tree
#   did not converge under UFBoot. This script writes to a new output folder and
#   does not overwrite the previous mixed-tree results unless FORCE=1 is set.
#
# Input:
#   $OUT/nifH_consensus_mafft_v2/nifH_consensus50_iupac_all_pass.mafft.fasta
#
# Output:
#   $OUT/nifH_consensus_iqtree_iupac_nm5000_v2/
#     nifH_consensus50_iupac_all_pass_nm5000.treefile
#     nifH_consensus50_iupac_all_pass_nm5000.contree
#     nifH_consensus50_iupac_all_pass_nm5000.iqtree
#     nifH_consensus50_iupac_all_pass_nm5000.log

module load SolisLemus-BioPhylo/2026.04.20 || true

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
ALIGN="${ALIGN:-$OUT/nifH_consensus_mafft_v2/nifH_consensus50_iupac_all_pass.mafft.fasta}"
TREE_DIR="${TREE_DIR:-$OUT/nifH_consensus_iqtree_iupac_nm5000_v2}"
PREFIX="$TREE_DIR/nifH_consensus50_iupac_all_pass_nm5000"
FORCE="${FORCE:-0}"

mkdir -p "$TREE_DIR"

if [[ ! -s "$ALIGN" ]]; then
  echo "ERROR: mixed-IUPAC alignment not found or empty: $ALIGN" >&2
  exit 1
fi

if [[ -s "$PREFIX.treefile" && "$FORCE" != "1" ]]; then
  echo "ERROR: output tree already exists: $PREFIX.treefile" >&2
  echo "Set FORCE=1 if you intentionally want to rerun and overwrite this output." >&2
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

echo "Building mixed-IUPAC nifH tree with extended UFBoot iterations"
echo "OUT:      $OUT"
echo "ALIGN:    $ALIGN"
echo "TREE_DIR: $TREE_DIR"
echo "PREFIX:   $PREFIX"
echo "IQ-TREE:  $IQTREE"
echo

echo "Input alignment check:"
echo "Sequences: $(grep -c '^>' "$ALIGN")"
awk '
  /^>/ {
    if (seq != "") print length(seq)
    seq=""
    next
  }
  { seq=seq $0 }
  END { if (seq != "") print length(seq) }
' "$ALIGN" | sort | uniq -c > "$TREE_DIR/iupac_alignment_length_distribution.txt"
cat "$TREE_DIR/iupac_alignment_length_distribution.txt"
echo

REDO_ARG=()
if [[ "$FORCE" == "1" ]]; then
  REDO_ARG=(-redo)
fi

"$IQTREE" \
  -s "$ALIGN" \
  -st DNA \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -nm 5000 \
  -T AUTO \
  --prefix "$PREFIX" \
  "${REDO_ARG[@]}" \
  | tee "$TREE_DIR/nifH_consensus50_iupac_all_pass_nm5000_iqtree_run.log"

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

# Internal node support values can still appear in some trees, so keep labels
# that match the sample-header pattern used in this project.
tips = [x for x in labels if "|nifH|consensus50|" in x]
print(len(tips))
PY

echo
echo "Mixed-IUPAC nm5000 tree complete."
echo "Main tree:      $PREFIX.treefile"
echo "Consensus tree: $PREFIX.contree"
echo "IQ-TREE report: $PREFIX.iqtree"
echo "Run log:        $PREFIX.log"
