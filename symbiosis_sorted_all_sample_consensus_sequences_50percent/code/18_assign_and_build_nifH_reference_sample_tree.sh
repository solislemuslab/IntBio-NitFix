#!/usr/bin/env bash
set -euo pipefail

# Step 18: assign sample-derived nifH consensus sequences to closest nifH references
# and build a combined reference + sample nifH tree.
#
# This script assumes nifH.fa.aln has already been copied to:
#   $OUT/reference_nifH_tree/nifH.fa.aln
#
# Main validation before/with tree construction:
#   1) BLAST sample consensus sequences against ungapped nifH reference sequences.
#   2) Save top 10 hits and best hit per sample with percent identity and alignment length.
#   3) Add sample sequences to the existing nifH reference alignment using MAFFT --add --keeplength.
#   4) Build a combined reference + sample tree with IQ-TREE.

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
REF_ALN="${REF_ALN:-$OUT/reference_nifH_tree/nifH.fa.aln}"
SAMPLES="${SAMPLES:-$OUT/nifH_consensus_tree_v2/nifH_consensus50_strict_single_dominant_pct80_depth10_Nle20.fasta}"
WORK="${WORK:-$OUT/nifH_reference_assignment_and_tree_v2}"
THREADS="${THREADS:-8}"
BUILD_TREE="${BUILD_TREE:-yes}"
TREE_BOOTSTRAP="${TREE_BOOTSTRAP:-1000}"
ALRT="${ALRT:-1000}"
TREE_NM="${TREE_NM:-5000}"

mkdir -p "$WORK" "$WORK/blast" "$WORK/alignment" "$WORK/tree" "$WORK/logs"

printf "Step 18: nifH sample-to-reference assignment and combined tree\n"
printf "OUT:       %s\n" "$OUT"
printf "REF_ALN:   %s\n" "$REF_ALN"
printf "SAMPLES:   %s\n" "$SAMPLES"
printf "WORK:      %s\n" "$WORK"
printf "THREADS:   %s\n" "$THREADS"
printf "TREE:      %s\n" "$BUILD_TREE"
printf "\n"

if [[ ! -s "$REF_ALN" ]]; then
  echo "ERROR: reference alignment not found or empty: $REF_ALN" >&2
  exit 1
fi
if [[ ! -s "$SAMPLES" ]]; then
  echo "ERROR: sample consensus FASTA not found or empty: $SAMPLES" >&2
  exit 1
fi

REF_COUNT=$(grep -c '^>' "$REF_ALN" || true)
SAMPLE_COUNT=$(grep -c '^>' "$SAMPLES" || true)
printf "Reference sequences: %s\n" "$REF_COUNT"
printf "Sample consensus sequences: %s\n\n" "$SAMPLE_COUNT"

# Create an ungapped reference FASTA for BLAST. BLAST performs local pairwise alignment;
# this is used as a closest-reference validation, separate from phylogenetic placement.
REF_UNGAPPED="$WORK/blast/nifH_reference_ungapped.fasta"
awk '
  /^>/ {
    if (seq != "") {
      gsub("-", "", seq)
      print seq
    }
    print
    seq=""
    next
  }
  {seq=seq $0}
  END {
    if (seq != "") {
      gsub("-", "", seq)
      print seq
    }
  }
' "$REF_ALN" > "$REF_UNGAPPED"

printf "Wrote ungapped BLAST reference: %s\n" "$REF_UNGAPPED"

if command -v makeblastdb >/dev/null 2>&1 && command -v blastn >/dev/null 2>&1; then
  makeblastdb \
    -in "$REF_UNGAPPED" \
    -dbtype nucl \
    -out "$WORK/blast/nifH_reference_db" \
    > "$WORK/logs/makeblastdb.log" 2>&1

  blastn \
    -query "$SAMPLES" \
    -db "$WORK/blast/nifH_reference_db" \
    -out "$WORK/blast/sample_nifH_vs_reference_top10.tsv" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -max_target_seqs 10 \
    -num_threads "$THREADS" \
    > "$WORK/logs/blastn.log" 2>&1

  sort -k1,1 -k12,12nr "$WORK/blast/sample_nifH_vs_reference_top10.tsv" \
    | awk '!seen[$1]++' \
    > "$WORK/blast/sample_nifH_best_reference_hit.tsv"

  awk -F'\t' 'BEGIN {
      OFS="\t";
      print "sample_id","best_reference_id","best_reference_taxon","percent_identity","alignment_length","query_coverage_percent","reference_coverage_percent","evalue","bitscore"
    }
    {
      tax=$2
      sub(/^nifH\|/, "", tax)
      sub(/\|.*/, "", tax)
      qcov=($14>0 ? 100*$4/$14 : "NA")
      scov=($15>0 ? 100*$4/$15 : "NA")
      print $1,$2,tax,$3,$4,qcov,scov,$11,$12
    }' "$WORK/blast/sample_nifH_best_reference_hit.tsv" \
    > "$WORK/blast/sample_nifH_best_reference_hit_with_taxon.tsv"

  printf "BLAST top10:    %s\n" "$WORK/blast/sample_nifH_vs_reference_top10.tsv"
  printf "BLAST best hit: %s\n" "$WORK/blast/sample_nifH_best_reference_hit_with_taxon.tsv"
else
  echo "WARNING: makeblastdb/blastn not found in PATH. Skipping BLAST validation." >&2
fi

# Add sample sequences to the existing reference alignment.
# --keeplength preserves reference alignment coordinates; insertions relative to the
# reference alignment are not retained. This is appropriate for placing samples in
# the reference coordinate system.
if ! command -v mafft >/dev/null 2>&1; then
  echo "ERROR: mafft is not available in PATH." >&2
  exit 1
fi

COMBINED_ALN="$WORK/alignment/nifH_reference_plus_strict_samples.mafft_add_keeplength.fasta"
mafft --add "$SAMPLES" --keeplength --thread "$THREADS" "$REF_ALN" \
  > "$COMBINED_ALN" \
  2> "$WORK/logs/mafft_add_keeplength.log"

printf "Combined alignment: %s\n" "$COMBINED_ALN"
printf "Combined sequences: %s\n" "$(grep -c '^>' "$COMBINED_ALN" || true)"

awk 'BEGIN{n=0} /^>/{if(n>0) print len; n++; len=0; next} {len+=length($0)} END{if(n>0) print len}' "$COMBINED_ALN" \
  | sort | uniq -c \
  > "$WORK/alignment/combined_alignment_lengths.txt"
printf "Alignment length check:\n"
cat "$WORK/alignment/combined_alignment_lengths.txt"

if [[ "$BUILD_TREE" == "yes" ]]; then
  if command -v iqtree >/dev/null 2>&1; then
    IQTREE="iqtree"
  elif command -v iqtree2 >/dev/null 2>&1; then
    IQTREE="iqtree2"
  else
    echo "ERROR: neither iqtree nor iqtree2 is available in PATH." >&2
    exit 1
  fi

  "$IQTREE" \
    -s "$COMBINED_ALN" \
    -st DNA \
    -m MFP \
    -B "$TREE_BOOTSTRAP" \
    -alrt "$ALRT" \
    -nm "$TREE_NM" \
    -T AUTO \
    --prefix "$WORK/tree/nifH_reference_plus_strict_samples" \
    2>&1 | tee "$WORK/logs/iqtree_reference_plus_samples_run.log"

  printf "\nTree file: %s\n" "$WORK/tree/nifH_reference_plus_strict_samples.treefile"
  printf "IQ-TREE report: %s\n" "$WORK/tree/nifH_reference_plus_strict_samples.iqtree"
fi

printf "\nDone.\n"
