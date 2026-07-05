#!/usr/bin/env bash
set -euo pipefail

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2}"
INDIR="$OUT/nifH_clean5_tree_ready_v2"
ALIGNDIR="$OUT/nifH_clean5_mafft_alignment_v2"

INPUT_FASTA="${INPUT_FASTA:-$INDIR/nifH_clean5_strict_single_dominant.fasta}"
LABEL="${LABEL:-strict_single_dominant}"

mkdir -p "$ALIGNDIR"

ALIGNED="$ALIGNDIR/nifH_clean5_${LABEL}.mafft.fasta"
SUMMARY="$ALIGNDIR/nifH_clean5_${LABEL}.alignment_summary.tsv"

echo "Step 13 V2: align nifH clean5 sequences with MAFFT"
echo "Input FASTA: $INPUT_FASTA"
echo "Output alignment: $ALIGNED"
echo

if ! command -v mafft >/dev/null 2>&1; then
  echo "ERROR: mafft is not available in PATH."
  echo "Load the module/environment that contains MAFFT, then run again."
  exit 1
fi

if [[ ! -s "$INPUT_FASTA" ]]; then
  echo "ERROR: input FASTA does not exist or is empty: $INPUT_FASTA"
  exit 1
fi

mafft --auto "$INPUT_FASTA" > "$ALIGNED"

python3 - "$ALIGNED" "$SUMMARY" <<'PY'
from pathlib import Path
import sys

fasta = Path(sys.argv[1])
summary = Path(sys.argv[2])

seq_count = 0
lengths = []
gap_counts = []
n_counts = []
name = None
chunks = []

def finish(chunks):
    seq = "".join(chunks).strip()
    if not seq:
        return
    lengths.append(len(seq))
    gap_counts.append(seq.count("-"))
    n_counts.append(seq.upper().count("N"))

with open(fasta) as handle:
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                finish(chunks)
            name = line[1:]
            seq_count += 1
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        finish(chunks)

alignment_length = lengths[0] if lengths else 0
all_same_length = len(set(lengths)) == 1
total_gaps = sum(gap_counts)
seqs_with_gaps = sum(1 for x in gap_counts if x > 0)
total_N = sum(n_counts)

with open(summary, "w") as out:
    out.write("metric\tvalue\n")
    out.write(f"sequence_count\t{seq_count}\n")
    out.write(f"alignment_length\t{alignment_length}\n")
    out.write(f"all_sequences_same_length\t{all_same_length}\n")
    out.write(f"total_gap_characters\t{total_gaps}\n")
    out.write(f"sequences_with_gaps\t{seqs_with_gaps}\n")
    out.write(f"total_N_characters\t{total_N}\n")

print(f"Sequence count: {seq_count}")
print(f"Alignment length: {alignment_length}")
print(f"All sequences same length: {all_same_length}")
print(f"Total gap characters: {total_gaps}")
print(f"Sequences with gaps: {seqs_with_gaps}")
print(f"Total N characters: {total_N}")
PY

echo
echo "Alignment complete."
echo "Aligned FASTA: $ALIGNED"
echo "Summary:       $SUMMARY"
