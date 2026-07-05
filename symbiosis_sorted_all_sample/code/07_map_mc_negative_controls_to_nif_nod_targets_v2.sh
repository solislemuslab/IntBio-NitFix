#!/bin/bash

set -uo pipefail

module load SolisLemus-BioPhylo/2026.04.20
module load samtools-1.9

OUTBASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
TRIMMED_DIR="$OUTBASE/symbiosis_trimmed_fastp"
REF="$OUTBASE/nif_nod_target_reference/all_central_nif_nod_targets.fasta"
OUT_DIR="$OUTBASE/nif_nod_mc_negative_control_mapping"
LOG_DIR="$OUT_DIR/logs"
SUMMARY="$OUT_DIR/nif_nod_mc_negative_control_mapping_summary.tsv"

mkdir -p "$OUT_DIR" "$LOG_DIR"

if ! command -v bwa >/dev/null 2>&1 || ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: bwa and/or samtools is unavailable."
    exit 1
fi

if [[ ! -s "$REF" ]]; then
    echo "ERROR: target reference not found:"
    echo "$REF"
    echo "Run Step 6 target extraction first."
    exit 1
fi

if [[ ! -s "$REF.bwt" || ! -s "$REF.sa" || ! -s "$REF.pac" || ! -s "$REF.ann" || ! -s "$REF.amb" ]]; then
    echo "Indexing nif/nod target reference:"
    echo "$REF"
    bwa index "$REF"
fi

shopt -s nullglob
p1_files=( "$TRIMMED_DIR"/MC*_P1.fastq.gz )

if [[ ${#p1_files[@]} -eq 0 ]]; then
    echo "ERROR: no MC control P1 files found in:"
    echo "$TRIMMED_DIR"
    exit 1
fi

echo "MC negative-control samples found: ${#p1_files[@]}"
echo "Reference: $REF"
echo "Output:    $OUT_DIR"
echo

for p1 in "${p1_files[@]}"; do
    sample=$(basename "$p1" _P1.fastq.gz)
    p2="$TRIMMED_DIR/${sample}_P2.fastq.gz"

    bam="$OUT_DIR/${sample}.bam"
    tmp="$OUT_DIR/${sample}.bam.tmp"
    flag="$LOG_DIR/${sample}.flagstat.txt"
    idx="$LOG_DIR/${sample}.idxstats.txt"
    bwa_log="$LOG_DIR/${sample}.bwa.log"

    if [[ ! -s "$p2" ]]; then
        echo "ERROR: missing paired read for $sample"
        continue
    fi

    if [[ -s "$bam" && -s "$bam.bai" && -s "$flag" && -s "$idx" ]]; then
        echo "Already completed, skipping: $sample"
        continue
    fi

    echo "Mapping MC negative control: $sample"
    rm -f "$tmp"

    if bwa mem -t 8 "$REF" "$p1" "$p2" 2> "$bwa_log" \
        | samtools sort -@ 4 -o "$tmp" -; then
        if samtools quickcheck "$tmp"; then
            mv "$tmp" "$bam"
            samtools index "$bam"
            samtools flagstat "$bam" > "$flag"
            samtools idxstats "$bam" > "$idx"
            echo "Completed: $sample"
        else
            echo "FAILED BAM check: $sample"
            rm -f "$tmp"
        fi
    else
        echo "FAILED mapping: $sample"
        rm -f "$tmp"
    fi
done

echo -e "sample\ttotal_reads\tmapped_reads\tmapped_percent\tproperly_paired\tproperly_paired_percent\ttop_target\ttop_target_mapped_reads" > "$SUMMARY"

for bam in "$OUT_DIR"/MC*.bam; do
    [[ -e "$bam" ]] || continue

    sample=$(basename "$bam" .bam)
    flag="$LOG_DIR/${sample}.flagstat.txt"
    idx="$LOG_DIR/${sample}.idxstats.txt"

    total=$(awk '/in total/ {print $1; exit}' "$flag")
    mapped=$(awk '/ mapped \(/ {print $1; exit}' "$flag")
    mapped_pct=$(awk '/ mapped \(/ {gsub(/[()%]/, "", $5); print $5; exit}' "$flag")
    proper=$(awk '/ properly paired / {print $1; exit}' "$flag")
    proper_pct=$(awk '/ properly paired / {gsub(/[()%]/, "", $6); print $6; exit}' "$flag")
    top_target=$(sort -k3,3nr "$idx" | head -1 | cut -f1)
    top_reads=$(sort -k3,3nr "$idx" | head -1 | cut -f3)

    echo -e "${sample}\t${total}\t${mapped}\t${mapped_pct}\t${proper}\t${proper_pct}\t${top_target}\t${top_reads}" >> "$SUMMARY"
done

echo
echo "MC negative-control nif/nod target mapping completed."
echo
column -t -s $'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
echo
echo "Output BAM folder: $OUT_DIR"
echo "Output log folder: $LOG_DIR"
echo "Summary table:     $SUMMARY"
