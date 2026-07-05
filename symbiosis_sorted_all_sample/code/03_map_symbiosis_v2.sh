#!/bin/bash

set -uo pipefail

module load SolisLemus-BioPhylo/2026.04.20
module load samtools-1.9

OUTBASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
TRIMMED_DIR="$OUTBASE/symbiosis_trimmed_fastp"
REF="$OUTBASE/reference/symbiosis_islands.fasta"
MAPPED_DIR="$OUTBASE/symbiosis_mapped_full"
LOGS_DIR="$OUTBASE/symbiosis_mapping_logs"
MAP_MANIFEST="$OUTBASE/metadata/symbiosis_v2_mapping_manifest.tsv"

mkdir -p "$MAPPED_DIR" "$LOGS_DIR" "$OUTBASE/metadata"

if ! command -v bwa >/dev/null 2>&1; then
    echo "ERROR: bwa is not available."
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available."
    exit 1
fi

if [[ ! -s "$REF" ]]; then
    echo "ERROR: V2 reference file not found:"
    echo "$REF"
    echo "Run Rscripts_v2/01_create_symbiosis_v2_manifest.sh first."
    exit 1
fi

if [[ ! -f "$REF.bwt" || ! -f "$REF.sa" || ! -f "$REF.pac" || ! -f "$REF.ann" || ! -f "$REF.amb" ]]; then
    echo "Indexing V2 reference: $REF"
    bwa index "$REF"
fi

printf "sample\tstatus\n" > "$MAP_MANIFEST"

shopt -s nullglob
p1_files=( "$TRIMMED_DIR"/*_P1.fastq.gz )

echo "Mapping V2 trimmed reads"
echo "Trimmed paired samples found: ${#p1_files[@]}"
echo "Reference: $REF"
echo "Output:    $MAPPED_DIR"
echo

for p1 in "${p1_files[@]}"; do
    sample=$(basename "$p1" _P1.fastq.gz)
    p2="$TRIMMED_DIR/${sample}_P2.fastq.gz"

    bam="$MAPPED_DIR/${sample}.bam"
    tmp_bam="$MAPPED_DIR/${sample}.bam.tmp"
    flagstat="$LOGS_DIR/${sample}.flagstat.txt"
    idxstats="$LOGS_DIR/${sample}.idxstats.txt"
    bwa_log="$LOGS_DIR/${sample}.bwa.stderr.log"
    sort_log="$LOGS_DIR/${sample}.samtools_sort.stderr.log"

    if [[ ! -s "$p2" ]]; then
        echo "ERROR: missing P2 for $sample"
        printf "%s\tmissing_P2\n" "$sample" >> "$MAP_MANIFEST"
        continue
    fi

    if [[ -s "$bam" && -s "$bam.bai" && -s "$flagstat" && -s "$idxstats" ]]; then
        echo "Already mapped, skipping: $sample"
        printf "%s\talready_completed\n" "$sample" >> "$MAP_MANIFEST"
        continue
    fi

    echo "Mapping: $sample"
    rm -f "$tmp_bam"

    if bwa mem -t 8 "$REF" "$p1" "$p2" 2> "$bwa_log" \
        | samtools sort -@ 4 -o "$tmp_bam" - 2> "$sort_log"; then

        if samtools quickcheck "$tmp_bam"; then
            mv "$tmp_bam" "$bam"
            samtools index "$bam"
            samtools flagstat "$bam" > "$flagstat"
            samtools idxstats "$bam" > "$idxstats"
            printf "%s\tcompleted\n" "$sample" >> "$MAP_MANIFEST"
        else
            echo "ERROR: quickcheck failed for $sample"
            rm -f "$tmp_bam"
            printf "%s\tquickcheck_failed\n" "$sample" >> "$MAP_MANIFEST"
        fi
    else
        echo "ERROR: mapping failed for $sample"
        rm -f "$tmp_bam"
        printf "%s\tfailed\n" "$sample" >> "$MAP_MANIFEST"
    fi
done

echo
echo "Mapping finished."
echo "Mapping manifest: $MAP_MANIFEST"
echo "Completed BAM files:"
find "$MAPPED_DIR" -maxdepth 1 -type f -name '*.bam' | wc -l
echo "Completed BAM indexes:"
find "$MAPPED_DIR" -maxdepth 1 -type f -name '*.bam.bai' | wc -l
echo "Completed flagstat reports:"
find "$LOGS_DIR" -maxdepth 1 -type f -name '*.flagstat.txt' | wc -l
