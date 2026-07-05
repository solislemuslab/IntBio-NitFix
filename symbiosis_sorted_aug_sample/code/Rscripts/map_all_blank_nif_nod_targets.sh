#!/bin/bash

set -uo pipefail

module load SolisLemus-BioPhylo/2026.04.20

BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"
TRIMMED="$BASE/symbiosis_trimmed_ryan_full"
REF="$BASE/nif_nod_target_reference/all_central_nif_nod_targets.fasta"

OUT="$BASE/nif_nod_blank_mapping"
LOG="$BASE/nif_nod_blank_mapping_logs"
SUMMARY="$BASE/nif_nod_blank_mapping_summary.tsv"

mkdir -p "$OUT" "$LOG"

if ! command -v bwa >/dev/null 2>&1 || ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: bwa and/or samtools is unavailable."
    exit 1
fi

if [[ ! -s "$REF" ]]; then
    echo "ERROR: Target reference not found:"
    echo "$REF"
    exit 1
fi

if [[ ! -s "$REF.bwt" ]]; then
    echo "Indexing target reference..."
    bwa index "$REF"
fi

shopt -s nullglob
P1_FILES=("$TRIMMED"/BLAN*_P1.fastq.gz)

if [[ ${#P1_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No BLAN P1 files were found."
    exit 1
fi

echo "Blank-control samples found: ${#P1_FILES[@]}"
echo "Reference: $REF"
echo

for P1 in "${P1_FILES[@]}"; do

    SAMPLE=$(basename "$P1" _P1.fastq.gz)
    P2="$TRIMMED/${SAMPLE}_P2.fastq.gz"

    BAM="$OUT/${SAMPLE}.bam"
    TMP="$OUT/${SAMPLE}.bam.tmp"
    FLAG="$LOG/${SAMPLE}.flagstat.txt"
    IDX="$LOG/${SAMPLE}.idxstats.txt"
    BWA_LOG="$LOG/${SAMPLE}.bwa.log"

    if [[ ! -s "$P2" ]]; then
        echo "ERROR: Missing paired read for $SAMPLE"
        continue
    fi

    if [[ -s "$BAM" && -s "$BAM.bai" && -s "$FLAG" && -s "$IDX" ]]; then
        echo "Already completed, skipping: $SAMPLE"
        continue
    fi

    echo "Mapping blank: $SAMPLE"

    rm -f "$TMP"

    if bwa mem -t 8 "$REF" "$P1" "$P2" 2> "$BWA_LOG" \
        | samtools sort -@ 4 -o "$TMP" -; then

        if samtools quickcheck "$TMP"; then
            mv "$TMP" "$BAM"
            samtools index "$BAM"
            samtools flagstat "$BAM" > "$FLAG"
            samtools idxstats "$BAM" > "$IDX"
            echo "Completed: $SAMPLE"
        else
            echo "FAILED BAM check: $SAMPLE"
            rm -f "$TMP"
        fi

    else
        echo "FAILED mapping: $SAMPLE"
        rm -f "$TMP"
    fi

done

echo -e "sample\ttotal_reads\tmapped_reads\tmapped_percent\tproperly_paired\tproperly_paired_percent\ttop_target\ttop_target_mapped_reads" > "$SUMMARY"

for BAM in "$OUT"/BLAN*.bam; do

    SAMPLE=$(basename "$BAM" .bam)
    FLAG="$LOG/${SAMPLE}.flagstat.txt"
    IDX="$LOG/${SAMPLE}.idxstats.txt"

    TOTAL=$(awk '/in total/ {print $1; exit}' "$FLAG")
    MAPPED=$(awk '/ mapped \(/ {print $1; exit}' "$FLAG")
    MAPPED_PCT=$(awk '/ mapped \(/ {gsub(/[()%]/, "", $5); print $5; exit}' "$FLAG")
    PROPER=$(awk '/ properly paired / {print $1; exit}' "$FLAG")
    PROPER_PCT=$(awk '/ properly paired / {gsub(/[()%]/, "", $6); print $6; exit}' "$FLAG")
    TOP_TARGET=$(sort -k3,3nr "$IDX" | head -1 | cut -f1)
    TOP_READS=$(sort -k3,3nr "$IDX" | head -1 | cut -f3)

    echo -e "${SAMPLE}\t${TOTAL}\t${MAPPED}\t${MAPPED_PCT}\t${PROPER}\t${PROPER_PCT}\t${TOP_TARGET}\t${TOP_READS}" >> "$SUMMARY"

done

echo
echo "Blank-control mapping completed."
echo
column -t -s $'\t' "$SUMMARY"
echo
echo "Output BAM folder: $OUT"
echo "Output log folder: $LOG"
echo "Summary table:     $SUMMARY"
