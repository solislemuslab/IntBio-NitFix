#!/usr/bin/env bash
set -euo pipefail

OUT=${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}
REF=${REF:-$OUT/reference/consensus_sequences_50percent.fasta}
BAMDIR=${BAMDIR:-$OUT/consensus_mapping_full}
COV=${COV:-$OUT/consensus_gene_coverage/consensus50_gene_coverage_all_samples.tsv}
SCRIPT=${SCRIPT:-$OUT/Rscripts_v2/13_extract_nifH_consensus_consensus50.py}
RUN_DIR=${RUN_DIR:-$OUT/nifH_consensus_tree_v1}

MIN_PCT=${MIN_PCT:-80}
MIN_MEAN_DEPTH=${MIN_MEAN_DEPTH:-10}
MIN_BASE_DEPTH=${MIN_BASE_DEPTH:-1}
MAX_N_PERCENT=${MAX_N_PERCENT:-20}

# Consensus confidence filters used by samtools mpileup.
# Q20 = base-call error probability about 1%; q10 keeps weak but nonzero mapping-confidence reads.
MIN_BASE_QUAL=${MIN_BASE_QUAL:-20}
MIN_MAP_QUAL=${MIN_MAP_QUAL:-10}

# Mixed-site rule: a secondary allele must pass all three thresholds.
MINOR_COUNT_THRESHOLD=${MINOR_COUNT_THRESHOLD:-5}
MINOR_FRACTION_THRESHOLD=${MINOR_FRACTION_THRESHOLD:-0.20}
MIXED_SITE_DEPTH_THRESHOLD=${MIXED_SITE_DEPTH_THRESHOLD:-20}
MIXED_POSITIONS_THRESHOLD=${MIXED_POSITIONS_THRESHOLD:-10}

JOBS=${JOBS:-2}
MAX_SAMPLES=${MAX_SAMPLES:-0}

PER_STRICT="$RUN_DIR/per_sample_strict_single_dominant"
PER_IUPAC="$RUN_DIR/per_sample_iupac_all_pass"
LOGDIR="$RUN_DIR/logs"
mkdir -p "$PER_STRICT" "$PER_IUPAC" "$LOGDIR"

SAMPLE_LIST="$RUN_DIR/nifH_good_samples_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}.txt"
QC="$RUN_DIR/nifH_consensus_qc.tsv"
STRICT_FASTA="$RUN_DIR/nifH_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.fasta"
IUPAC_FASTA="$RUN_DIR/nifH_consensus50_iupac_all_pass_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.fasta"

awk -F'\t' -v minpct="$MIN_PCT" -v mindepth="$MIN_MEAN_DEPTH" '
NR==1 {next}
$3=="nifH" && $6>=minpct && $7>=mindepth {print $1}
' "$COV" | sort -u > "$SAMPLE_LIST"

if [[ "$MAX_SAMPLES" != "0" ]]; then
  head -n "$MAX_SAMPLES" "$SAMPLE_LIST" > "$SAMPLE_LIST.tmp"
  mv "$SAMPLE_LIST.tmp" "$SAMPLE_LIST"
fi

echo "Extracting nifH consensus sequences from consensus_sequences_50percent mapping"
echo "OUT: $OUT"
echo "Reference: $REF"
echo "Coverage table: $COV"
echo "Selected samples: $(wc -l < "$SAMPLE_LIST")"
echo "Sample threshold: percent_covered >= $MIN_PCT; mean_depth >= $MIN_MEAN_DEPTH; N <= $MAX_N_PERCENT%"
echo "Pileup filters: base quality >= $MIN_BASE_QUAL; mapping quality >= $MIN_MAP_QUAL"
echo "Mixed rule: depth >= $MIXED_SITE_DEPTH_THRESHOLD; minor count >= $MINOR_COUNT_THRESHOLD; minor fraction >= $MINOR_FRACTION_THRESHOLD; mixed positions > $MIXED_POSITIONS_THRESHOLD"
echo "Jobs: $JOBS"

echo -e "sample\ttarget\tlength\tcovered_bases\tpercent_covered_from_pileup\tmean_depth_from_pileup\tmedian_depth_from_pileup\tN_count\tN_percent\tmixed_positions\tmax_minor_allele_fraction\tmin_base_qual\tmin_map_qual\tminor_count_threshold\tminor_fraction_threshold\tmixed_site_depth_threshold\tmixed_positions_threshold\tstatus" > "$QC"

export REF BAMDIR PER_STRICT PER_IUPAC SCRIPT MIN_BASE_DEPTH MAX_N_PERCENT MIN_BASE_QUAL MIN_MAP_QUAL
export MINOR_COUNT_THRESHOLD MINOR_FRACTION_THRESHOLD MIXED_SITE_DEPTH_THRESHOLD MIXED_POSITIONS_THRESHOLD

cat "$SAMPLE_LIST" | parallel -j "$JOBS" --joblog "$LOGDIR/nifH_consensus_parallel.joblog" '
  sample={}
  bam="$BAMDIR/${sample}.bam"
  strict="$PER_STRICT/${sample}.nifH.strict_single_dominant.fasta"
  iupac="$PER_IUPAC/${sample}.nifH.iupac_all_pass.fasta"
  if [[ ! -s "$bam" ]]; then
    echo -e "$sample\tnifH.fa\t0\t0\t0\t0\t0\t0\t100\t0\t0\t$MIN_BASE_QUAL\t$MIN_MAP_QUAL\t$MINOR_COUNT_THRESHOLD\t$MINOR_FRACTION_THRESHOLD\t$MIXED_SITE_DEPTH_THRESHOLD\t$MIXED_POSITIONS_THRESHOLD\tfail_missing_bam"
  else
    python3 "$SCRIPT" "$sample" "$bam" "$REF" "$strict" "$iupac" "$MIN_BASE_DEPTH" "$MAX_N_PERCENT" "$MIN_BASE_QUAL" "$MIN_MAP_QUAL" "$MINOR_COUNT_THRESHOLD" "$MINOR_FRACTION_THRESHOLD" "$MIXED_SITE_DEPTH_THRESHOLD" "$MIXED_POSITIONS_THRESHOLD"
  fi
' >> "$QC"

find "$PER_STRICT" -name "*.fasta" -size +0 -print0 | sort -z | xargs -0 cat > "$STRICT_FASTA" || true
find "$PER_IUPAC" -name "*.fasta" -size +0 -print0 | sort -z | xargs -0 cat > "$IUPAC_FASTA" || true

echo
echo "Done."
echo "QC table: $QC"
echo "Strict single-dominant FASTA: $STRICT_FASTA"
echo "Mixed-aware IUPAC FASTA: $IUPAC_FASTA"
echo
echo "Status counts:"
tail -n +2 "$QC" | cut -f18 | sort | uniq -c

echo
echo "Sequence counts:"
echo -n "strict_single_dominant\t"; grep -c '^>' "$STRICT_FASTA" || true
echo -n "iupac_all_pass\t"; grep -c '^>' "$IUPAC_FASTA" || true
