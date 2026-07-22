#!/usr/bin/env bash
set -euo pipefail

# Step 19: BLAST strict sample-derived nifH consensus sequences against
# the 746-sequence nifH reference set.
#
# Purpose:
#   The strict nifH tree tips are sample-derived consensus DNA sequences.
#   This script assigns each strict sample consensus sequence to its closest
#   known/reference nifH sequence so the tree can be interpreted in terms of
#   likely bacterial lineages/taxa.
#
# Primary sample set:
#   nifH consensus50 strict single-dominant sequences
#   threshold: percent covered >= 80, mean depth >= 10, N <= 20%
#   expected sequence count: 387
#
# Validation columns produced:
#   percent_identity
#   alignment_length
#   query_coverage_percent      = how much of the sample consensus aligned
#   reference_coverage_percent  = how much of the reference sequence aligned
#   evalue
#   bitscore

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
WORK="${WORK:-$OUT/nifH_reference_assignment_and_tree_v2}"

REF_ALN="${REF_ALN:-$OUT/reference_nifH_tree/nifH.fa.aln}"
REF_UNGAPPED="${REF_UNGAPPED:-$WORK/blast/nifH_reference_ungapped.fasta}"
SAMPLES="${SAMPLES:-$OUT/nifH_consensus_tree_v2/nifH_consensus50_strict_single_dominant_pct80_depth10_Nle20.fasta}"

THREADS="${THREADS:-8}"

mkdir -p "$WORK/blast" "$WORK/logs"

echo "Step 19: BLAST strict sample nifH consensus sequences to reference nifH taxa"
echo "OUT:           $OUT"
echo "WORK:          $WORK"
echo "Reference aln: $REF_ALN"
echo "Reference DB:  $REF_UNGAPPED"
echo "Samples:       $SAMPLES"
echo "Threads:       $THREADS"
echo

if [[ ! -s "$REF_ALN" ]]; then
  echo "ERROR: reference alignment not found: $REF_ALN" >&2
  exit 1
fi

if [[ ! -s "$SAMPLES" ]]; then
  echo "ERROR: sample FASTA not found: $SAMPLES" >&2
  exit 1
fi

if ! command -v blastn >/dev/null 2>&1; then
  echo "ERROR: blastn is not available in PATH." >&2
  echo "Try: module load blast+/2.17.0" >&2
  exit 1
fi

if ! command -v makeblastdb >/dev/null 2>&1; then
  echo "ERROR: makeblastdb is not available in PATH." >&2
  echo "Try: module load blast+/2.17.0" >&2
  exit 1
fi

echo "BLAST version:"
blastn -version
echo

echo "Reference sequences in alignment:"
grep -c '^>' "$REF_ALN"
echo "Strict sample consensus sequences:"
grep -c '^>' "$SAMPLES"
echo

# Remove alignment gaps from reference sequences for BLAST.
# BLAST performs local pairwise alignment, so it should use ungapped sequences.
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

echo "Wrote ungapped reference FASTA:"
echo "$REF_UNGAPPED"
echo

makeblastdb \
  -in "$REF_UNGAPPED" \
  -dbtype nucl \
  -out "$WORK/blast/nifH_reference_db" \
  > "$WORK/logs/step19_makeblastdb.log" 2>&1

blastn \
  -query "$SAMPLES" \
  -db "$WORK/blast/nifH_reference_db" \
  -out "$WORK/blast/sample_nifH_vs_reference_top10.tsv" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  -max_target_seqs 10 \
  -num_threads "$THREADS" \
  > "$WORK/logs/step19_blastn.log" 2>&1

# Best hit per sample by highest bitscore.
sort -k1,1 -k12,12nr "$WORK/blast/sample_nifH_vs_reference_top10.tsv" \
  | awk '!seen[$1]++' \
  > "$WORK/blast/sample_nifH_best_reference_hit.tsv"

# Add parsed taxon and coverage-style summaries.
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

echo "Done."
echo
echo "Top 10 BLAST hits:"
echo "$WORK/blast/sample_nifH_vs_reference_top10.tsv"
echo "Best reference hit per sample:"
echo "$WORK/blast/sample_nifH_best_reference_hit.tsv"
echo "Best reference hit with parsed taxon:"
echo "$WORK/blast/sample_nifH_best_reference_hit_with_taxon.tsv"
echo
echo "Line counts:"
wc -l "$WORK/blast/sample_nifH_vs_reference_top10.tsv"
wc -l "$WORK/blast/sample_nifH_best_reference_hit_with_taxon.tsv"
echo
echo "First rows:"
head -10 "$WORK/blast/sample_nifH_best_reference_hit_with_taxon.tsv"
