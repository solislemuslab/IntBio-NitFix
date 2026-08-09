#!/usr/bin/env bash
set -euo pipefail

# nolF: extract strict and mixed-IUPAC consensus sequences from good-coverage BAMs.
# Good coverage is selected from consensus50_gene_coverage_all_samples.tsv:
#   gene == nolF, percent_covered >= 80, mean_depth >= 10.
# Consensus QC then requires N_percent <= 20. Per-position consensus calls require >=1 high-quality base; the sample-level depth filter is mean_depth >=10.

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
GENE="nolF"
REF="${OUT}/reference/consensus_sequences_50percent.fasta"
COV="${OUT}/consensus_gene_coverage/consensus50_gene_coverage_all_samples.tsv"
MAP="${OUT}/consensus_mapping_full"
RUN_DIR="${OUT}/gene_trees/nolF/01_consensus"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/01_extract_nolF_consensus_consensus50.py"

MIN_PCT="${MIN_PCT:-80}"
MIN_MEAN_DEPTH="${MIN_MEAN_DEPTH:-10}"
MIN_BASE_DEPTH="${MIN_BASE_DEPTH:-1}"
MAX_N_PERCENT="${MAX_N_PERCENT:-20}"
MIN_BASE_QUAL="${MIN_BASE_QUAL:-20}"
MIN_MAP_QUAL="${MIN_MAP_QUAL:-10}"
MINOR_COUNT_THRESHOLD="${MINOR_COUNT_THRESHOLD:-5}"
MINOR_FRACTION_THRESHOLD="${MINOR_FRACTION_THRESHOLD:-0.20}"
MIXED_SITE_DEPTH_THRESHOLD="${MIXED_SITE_DEPTH_THRESHOLD:-20}"
MIXED_POSITIONS_THRESHOLD="${MIXED_POSITIONS_THRESHOLD:-10}"
JOBS="${JOBS:-2}"
MAX_SAMPLES="${MAX_SAMPLES:-0}"

mkdir -p "$RUN_DIR/tmp" "$RUN_DIR/logs"

STRICT_FASTA="$RUN_DIR/nolF_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.fasta"
IUPAC_FASTA="$RUN_DIR/nolF_consensus50_iupac_all_pass_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.fasta"
QC_TABLE="$RUN_DIR/nolF_consensus_qc.tsv"
SELECTED="$RUN_DIR/nolF_selected_samples.tsv"

: > "$STRICT_FASTA"
: > "$IUPAC_FASTA"

awk -F'	' -v gene="$GENE" -v minpct="$MIN_PCT" -v mindepth="$MIN_MEAN_DEPTH" '
  NR>1 && $3==gene && $6>=minpct && $7>=mindepth {print $1}
' "$COV" | sort -u > "$SELECTED"

if [[ "$MAX_SAMPLES" != "0" ]]; then
  head -n "$MAX_SAMPLES" "$SELECTED" > "$SELECTED.tmp"
  mv "$SELECTED.tmp" "$SELECTED"
fi

N=$(wc -l < "$SELECTED" | tr -d ' ')
echo "nolF: selected samples: $N"
echo "nolF: output: $RUN_DIR"

printf "sample	target	length	covered_bases	percent_covered_from_pileup	mean_depth_from_pileup	median_depth_from_pileup	N_count	N_percent	mixed_positions	max_minor_allele_fraction	min_base_qual	min_map_qual	minor_count_threshold	minor_fraction_threshold	mixed_site_depth_threshold	mixed_positions_threshold	status
" > "$QC_TABLE"

run_one () {
  sample="$1"
  bam="$MAP/${sample}.bam"
  if [[ ! -s "$bam" ]]; then
    printf "%s	%s.fa	0	0	0	0	0	0	100	0	0	%s	%s	%s	%s	%s	%s	fail_missing_bam
"       "$sample" "$GENE" "$MIN_BASE_QUAL" "$MIN_MAP_QUAL" "$MINOR_COUNT_THRESHOLD" "$MINOR_FRACTION_THRESHOLD" "$MIXED_SITE_DEPTH_THRESHOLD" "$MIXED_POSITIONS_THRESHOLD"
    return 0
  fi

  tmp_strict="$RUN_DIR/tmp/${sample}.strict.fasta"
  tmp_iupac="$RUN_DIR/tmp/${sample}.iupac.fasta"
  python3 "$PYTHON_SCRIPT" "$GENE" "$sample" "$bam" "$REF" "$tmp_strict" "$tmp_iupac"     "$MIN_BASE_DEPTH" "$MAX_N_PERCENT" "$MIN_BASE_QUAL" "$MIN_MAP_QUAL"     "$MINOR_COUNT_THRESHOLD" "$MINOR_FRACTION_THRESHOLD" "$MIXED_SITE_DEPTH_THRESHOLD" "$MIXED_POSITIONS_THRESHOLD"
}
export -f run_one
export OUT GENE REF COV MAP RUN_DIR PYTHON_SCRIPT STRICT_FASTA IUPAC_FASTA QC_TABLE
export MIN_BASE_DEPTH MAX_N_PERCENT MIN_BASE_QUAL MIN_MAP_QUAL MINOR_COUNT_THRESHOLD MINOR_FRACTION_THRESHOLD MIXED_SITE_DEPTH_THRESHOLD MIXED_POSITIONS_THRESHOLD

parallel -j "$JOBS" run_one :::: "$SELECTED" >> "$QC_TABLE"

find "$RUN_DIR/tmp" -name "*.strict.fasta" -size +0 -print0 | sort -z | xargs -0 cat >> "$STRICT_FASTA" || true
find "$RUN_DIR/tmp" -name "*.iupac.fasta" -size +0 -print0 | sort -z | xargs -0 cat >> "$IUPAC_FASTA" || true

rm -rf "$RUN_DIR/tmp"

echo "Status counts for nolF:"
tail -n +2 "$QC_TABLE" | cut -f18 | sort | uniq -c

echo "strict $(grep -c '^>' "$STRICT_FASTA" || true)"
echo "iupac  $(grep -c '^>' "$IUPAC_FASTA" || true)"
