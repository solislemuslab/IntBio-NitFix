#!/usr/bin/env bash
set -euo pipefail

echo "Mapping all V2 trimmed symbiosis reads to consensus_sequences_50percent reference"

OLD="${OLD:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2}"
OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"

REF="$OUT/reference/consensus_sequences_50percent.fasta"
TRIM="$OLD/symbiosis_trimmed_fastp"
MAP="$OUT/consensus_mapping_full"
LOG="$OUT/consensus_mapping_logs"
META="$OUT/metadata"

JOBS="${JOBS:-4}"
BWA_THREADS="${BWA_THREADS:-4}"
SORT_THREADS="${SORT_THREADS:-2}"
MAX_SAMPLES="${MAX_SAMPLES:-0}"

mkdir -p "$MAP" "$LOG" "$META"

SAMPLES="$META/consensus50_mapping_sample_list.tsv"
MANIFEST="$META/consensus50_mapping_manifest.tsv"
JOBLOG="$LOG/consensus50_parallel_mapping.joblog"

echo "Old V2 folder:      $OLD"
echo "New output folder:  $OUT"
echo "Reference:          $REF"
echo "Trimmed reads:      $TRIM"
echo "Mapping output:     $MAP"
echo "Logs:               $LOG"
echo "GNU parallel jobs:  $JOBS"
echo "BWA threads/job:    $BWA_THREADS"
echo "Sort threads/job:   $SORT_THREADS"
echo "Max samples:        $MAX_SAMPLES"
echo

if [[ ! -s "$REF" ]]; then
  echo "ERROR: reference FASTA not found or empty: $REF" >&2
  exit 1
fi

if [[ ! -d "$TRIM" ]]; then
  echo "ERROR: trimmed-read folder not found: $TRIM" >&2
  echo "This script expects the corrected V2 trimmed folder: symbiosis_trimmed_fastp" >&2
  exit 1
fi

command -v bwa >/dev/null 2>&1 || { echo "ERROR: bwa not found in PATH" >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH" >&2; exit 1; }
command -v parallel >/dev/null 2>&1 || { echo "ERROR: GNU parallel not found in PATH" >&2; exit 1; }

if [[ ! -s "$REF.bwt" ]]; then
  echo "Indexing reference with bwa index..."
  bwa index "$REF"
fi

echo -e "sample\tP1\tP2" > "$SAMPLES"
find "$TRIM" -maxdepth 1 -name "*_P1.fastq.gz" -print \
  | sort \
  | while read -r P1
do
  sample=$(basename "$P1" "_P1.fastq.gz")
  P2="$TRIM/${sample}_P2.fastq.gz"
  if [[ -s "$P2" ]]; then
    echo -e "${sample}\t${P1}\t${P2}"
  else
    echo "WARNING: missing P2 for $sample" >&2
  fi
done >> "$SAMPLES"

sample_count=$(tail -n +2 "$SAMPLES" | wc -l | tr -d ' ')
echo "Paired trimmed samples found: $sample_count"

if [[ "$sample_count" -eq 0 ]]; then
  echo "ERROR: no paired trimmed samples found in $TRIM" >&2
  exit 1
fi

RUN_LIST="$META/consensus50_mapping_sample_list_for_this_run.tsv"
if [[ "$MAX_SAMPLES" != "0" ]]; then
  { head -1 "$SAMPLES"; tail -n +2 "$SAMPLES" | head -n "$MAX_SAMPLES"; } > "$RUN_LIST"
else
  cp "$SAMPLES" "$RUN_LIST"
fi

run_count=$(tail -n +2 "$RUN_LIST" | wc -l | tr -d ' ')
echo "Samples selected for this run: $run_count"
echo

map_one_sample() {
  sample="$1"
  P1="$2"
  P2="$3"

  bam="$MAP/${sample}.bam"
  bai="$MAP/${sample}.bam.bai"
  flagstat="$LOG/${sample}.flagstat.txt"
  bwa_log="$LOG/${sample}.bwa.log"

  if [[ -s "$bam" && -s "$bai" && -s "$flagstat" ]]; then
    echo -e "${sample}\tskipped_existing"
    return 0
  fi

  tmp_bam="$MAP/${sample}.tmp.bam"

  bwa mem -t "$BWA_THREADS" "$REF" "$P1" "$P2" 2> "$bwa_log" \
    | samtools sort -@ "$SORT_THREADS" -o "$tmp_bam"

  mv "$tmp_bam" "$bam"
  samtools index "$bam"
  samtools flagstat "$bam" > "$flagstat"

  echo -e "${sample}\tcompleted"
}

export REF MAP LOG BWA_THREADS SORT_THREADS
export -f map_one_sample

tail -n +2 "$RUN_LIST" \
  | parallel --colsep '\t' --jobs "$JOBS" --joblog "$JOBLOG" map_one_sample {1} {2} {3} \
  > "$MANIFEST"

echo
echo "Mapping complete."
echo "Sample list:       $SAMPLES"
echo "Run list:          $RUN_LIST"
echo "Mapping manifest:  $MANIFEST"
echo "Parallel job log:  $JOBLOG"
echo

echo "Status counts:"
cut -f2 "$MANIFEST" | sort | uniq -c
echo

echo "Completed BAM files:"
find "$MAP" -maxdepth 1 -name "*.bam" ! -name "*.tmp.bam" | wc -l
echo "Completed BAM indexes:"
find "$MAP" -maxdepth 1 -name "*.bam.bai" | wc -l
echo "Completed flagstat reports:"
find "$LOG" -maxdepth 1 -name "*.flagstat.txt" | wc -l
