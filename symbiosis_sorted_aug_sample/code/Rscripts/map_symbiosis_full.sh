#!/bin/bash

set -uo pipefail

module load SolisLemus-BioPhylo/2026.04.20
module load samtools-1.9

BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"
TRIMMED_DIR="$BASE/symbiosis_trimmed_ryan_full"
REF="$BASE/symbiosis_islands.fasta"
MAPPED_DIR="$BASE/symbiosis_mapped_full"
LOGS_DIR="$BASE/symbiosis_mapping_logs_full"
MANIFEST="$BASE/symbiosis_mapping_full_manifest.tsv"

mkdir -p "$MAPPED_DIR" "$LOGS_DIR"

if ! command -v bwa >/dev/null 2>&1; then
    echo "ERROR: bwa is not available."
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available."
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference file not found: $REF"
    exit 1
fi

# Create the BWA index once, if it is not already present.
if [[ ! -f "$REF.bwt" || ! -f "$REF.sa" || ! -f "$REF.pac" || ! -f "$REF.ann" || ! -f "$REF.amb" ]]; then
    echo "Indexing reference: $REF"
    bwa index "$REF"
fi

if [[ ! -f "$MANIFEST" ]]; then
    printf "sample\tstatus\n" > "$MANIFEST"
fi

shopt -s nullglob
p1_files=( "$TRIMMED_DIR"/*_P1.fastq.gz )

echo "Total trimmed paired samples found: ${#p1_files[@]}"
echo "Reads folder:   $TRIMMED_DIR"
echo "Reference:      $REF"
echo "Mapping output: $MAPPED_DIR"
echo "BLAN control samples are retained for later QC."
echo

for p1 in "${p1_files[@]}"; do
    sample=$(basename "$p1" _P1.fastq.gz)
    p2="$TRIMMED_DIR/${sample}_P2.fastq.gz"

    bam="$MAPPED_DIR/${sample}.bam"
    tmp_bam="$MAPPED_DIR/${sample}.bam.tmp"
    flagstat="$LOGS_DIR/${sample}.flagstat.txt"
    bwa_log="$LOGS_DIR/${sample}.bwa.stderr.log"
    sort_log="$LOGS_DIR/${sample}.samtools_sort.stderr.log"

    if [[ ! -f "$p2" ]]; then
        echo "ERROR: Missing P2 for $sample"
        printf "%s\tmissing_P2\n" "$sample" >> "$MANIFEST"
        continue
    fi

    if [[ -s "$bam" && -s "$bam.bai" && -s "$flagstat" ]]; then
        echo "Already completed, skipping: $sample"
        continue
    fi

    echo "Mapping sample: $sample"

    rm -f "$tmp_bam"

    if bwa mem -t 8 "$REF" "$p1" "$p2" 2> "$bwa_log" \
        | samtools sort -@ 4 -o "$tmp_bam" - 2> "$sort_log"; then

        mv "$tmp_bam" "$bam"
        samtools index "$bam"
        samtools flagstat "$bam" > "$flagstat"

        echo "Completed: $sample"
        printf "%s\tcompleted\n" "$sample" >> "$MANIFEST"
    else
        echo "ERROR: Mapping failed for $sample. Check logs."
        rm -f "$tmp_bam"
        printf "%s\tfailed\n" "$sample" >> "$MANIFEST"
    fi
done

echo
echo "=========================================="
echo "Mapping run finished."
echo "Mapped BAM files: $MAPPED_DIR"
echo "Logs:             $LOGS_DIR"
echo "Manifest:         $MANIFEST"
echo "=========================================="

echo "Completed BAM files:"
find "$MAPPED_DIR" -maxdepth 1 -type f -name '*.bam' | wc -l

echo "Completed BAM indexes:"
find "$MAPPED_DIR" -maxdepth 1 -type f -name '*.bam.bai' | wc -l

echo "Completed flagstat reports:"
find "$LOGS_DIR" -maxdepth 1 -type f -name '*.flagstat.txt' | wc -l
