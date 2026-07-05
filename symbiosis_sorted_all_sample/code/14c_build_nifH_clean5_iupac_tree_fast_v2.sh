#!/usr/bin/env bash
set -euo pipefail

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2}"
ALIGNDIR="$OUT/nifH_clean5_mafft_alignment_v2"
TREEDIR="$OUT/nifH_clean5_iqtree_v2"
INPUT_FASTA="$OUT/nifH_clean5_tree_ready_v2/nifH_clean5_iupac_all_pass.fasta"
ALIGNMENT="$ALIGNDIR/nifH_clean5_iupac_all_pass.mafft.fasta"
PREFIX="$TREEDIR/nifH_clean5_iupac_all_pass_fast_noboot"

mkdir -p "$ALIGNDIR" "$TREEDIR"

echo "Fast mixed-IUPAC nifH clean5 tree"
echo "This is a quick feedback tree: fixed model GTR+F+R9, no bootstrap."
echo "Input FASTA: $INPUT_FASTA"
echo "Alignment:   $ALIGNMENT"
echo "Tree prefix: $PREFIX"
echo

if [[ ! -s "$INPUT_FASTA" ]]; then
  echo "ERROR: input FASTA does not exist or is empty: $INPUT_FASTA"
  exit 1
fi

if [[ ! -s "$ALIGNMENT" ]]; then
  if ! command -v mafft >/dev/null 2>&1; then
    echo "ERROR: mafft is not available in PATH."
    exit 1
  fi
  echo "MAFFT alignment does not exist yet. Running MAFFT..."
  mafft --auto "$INPUT_FASTA" > "$ALIGNMENT"
else
  echo "Existing MAFFT alignment found. Reusing it."
fi

python3 - "$ALIGNMENT" "$ALIGNDIR/nifH_clean5_iupac_all_pass.fast_alignment_summary.tsv" <<'PY'
from pathlib import Path
import sys

fasta = Path(sys.argv[1])
summary = Path(sys.argv[2])

count = 0
lengths = []
gaps = 0
ns = 0
name = None
chunks = []

def finish(seq):
    global gaps, ns
    lengths.append(len(seq))
    gaps += seq.count("-")
    ns += seq.upper().count("N")

with open(fasta) as handle:
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                finish("".join(chunks))
            name = line[1:]
            count += 1
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        finish("".join(chunks))

with open(summary, "w") as out:
    out.write("metric\tvalue\n")
    out.write(f"sequence_count\t{count}\n")
    out.write(f"alignment_length\t{lengths[0] if lengths else 0}\n")
    out.write(f"all_sequences_same_length\t{len(set(lengths)) == 1}\n")
    out.write(f"total_gap_characters\t{gaps}\n")
    out.write(f"total_N_characters\t{ns}\n")

print(f"Sequence count: {count}")
print(f"Alignment length: {lengths[0] if lengths else 0}")
print(f"Total gap characters: {gaps}")
print(f"Total N characters: {ns}")
PY

if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE="iqtree"
else
  echo "ERROR: neither iqtree2 nor iqtree is available in PATH."
  exit 1
fi

echo
echo "Running fast IQ-TREE without bootstrap..."
"$IQTREE" \
  -s "$ALIGNMENT" \
  -m GTR+F+R9 \
  -T AUTO \
  --prefix "$PREFIX"

echo
echo "Fast mixed-IUPAC tree complete."
echo "Tree: ${PREFIX}.treefile"
echo "Log:  ${PREFIX}.log"
