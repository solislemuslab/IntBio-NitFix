#!/bin/bash

set -uo pipefail

module load SolisLemus-BioPhylo/2026.04.20
module load samtools-1.9

OUTBASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
TRIMMED_DIR="$OUTBASE/symbiosis_trimmed_fastp"
REF="$OUTBASE/reference/symbiosis_islands.fasta"
MAPPED_DIR="$OUTBASE/symbiosis_mapped_full"
LOGS_DIR="$OUTBASE/symbiosis_mapping_logs"
STATUS_DIR="$LOGS_DIR/per_sample_status"
MAP_MANIFEST="$OUTBASE/metadata/symbiosis_v2_mapping_manifest.tsv"

# Keep this conservative unless Ryan/Solis-Lemus cluster policy says more is OK.
JOBS="${JOBS:-3}"
BWA_THREADS="${BWA_THREADS:-4}"
SORT_THREADS="${SORT_THREADS:-2}"
MAX_SAMPLES="${MAX_SAMPLES:-0}"

mkdir -p "$MAPPED_DIR" "$LOGS_DIR" "$STATUS_DIR" "$OUTBASE/metadata"

if ! command -v bwa >/dev/null 2>&1; then
    echo "ERROR: bwa is not available."
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available."
    exit 1
fi

if ! command -v parallel >/dev/null 2>&1; then
    echo "ERROR: GNU parallel is not available."
    exit 1
fi

if [[ ! -s "$REF" ]]; then
    echo "ERROR: V2 reference file not found:"
    echo "$REF"
    echo "Run Rscripts_v2/01_create_symbiosis_v2_manifest.sh first."
    exit 1
fi

# Build the shared reference index once before any parallel jobs start.
if [[ ! -f "$REF.bwt" || ! -f "$REF.sa" || ! -f "$REF.pac" || ! -f "$REF.ann" || ! -f "$REF.amb" ]]; then
    echo "Indexing V2 reference once before parallel mapping: $REF"
    bwa index "$REF"
fi

shopt -s nullglob
p1_files=( "$TRIMMED_DIR"/*_P1.fastq.gz )
samples=()

for p1 in "${p1_files[@]}"; do
    sample=$(basename "$p1" _P1.fastq.gz)
    samples+=( "$sample" )
done

if [[ "${#samples[@]}" -eq 0 ]]; then
    echo "ERROR: no trimmed P1 files found in $TRIMMED_DIR"
    exit 1
fi

if [[ "$MAX_SAMPLES" -gt 0 ]]; then
    samples=( "${samples[@]:0:$MAX_SAMPLES}" )
fi

echo "Parallel mapping V2 trimmed reads"
echo "Trimmed paired samples found: ${#samples[@]}"
echo "Reference: $REF"
echo "Output:    $MAPPED_DIR"
echo "GNU parallel jobs: $JOBS"
echo "BWA threads per job: $BWA_THREADS"
echo "samtools sort threads per job: $SORT_THREADS"
echo "Max samples for this run: $MAX_SAMPLES"
echo

map_one_sample() {
    sample="$1"

    p1="$TRIMMED_DIR/${sample}_P1.fastq.gz"
    p2="$TRIMMED_DIR/${sample}_P2.fastq.gz"

    bam="$MAPPED_DIR/${sample}.bam"
    tmp_bam="$MAPPED_DIR/${sample}.bam.tmp.$$"
    flagstat="$LOGS_DIR/${sample}.flagstat.txt"
    idxstats="$LOGS_DIR/${sample}.idxstats.txt"
    bwa_log="$LOGS_DIR/${sample}.bwa.stderr.log"
    sort_log="$LOGS_DIR/${sample}.samtools_sort.stderr.log"
    status_file="$STATUS_DIR/${sample}.status.tsv"

    if [[ ! -s "$p1" || ! -s "$p2" ]]; then
        printf "%s\tmissing_trimmed_pair\n" "$sample" > "$status_file"
        echo "ERROR: missing P1/P2 for $sample"
        return 0
    fi

    if [[ -s "$bam" && -s "$bam.bai" && -s "$flagstat" && -s "$idxstats" ]]; then
        printf "%s\talready_completed\n" "$sample" > "$status_file"
        echo "Already mapped, skipping: $sample"
        return 0
    fi

    echo "Mapping: $sample"
    rm -f "$tmp_bam"

    if bwa mem -t "$BWA_THREADS" "$REF" "$p1" "$p2" 2> "$bwa_log" \
        | samtools sort -@ "$SORT_THREADS" -o "$tmp_bam" - 2> "$sort_log"; then

        if samtools quickcheck "$tmp_bam"; then
            mv "$tmp_bam" "$bam"
            samtools index "$bam"
            samtools flagstat "$bam" > "$flagstat"
            samtools idxstats "$bam" > "$idxstats"
            printf "%s\tcompleted\n" "$sample" > "$status_file"
        else
            echo "ERROR: quickcheck failed for $sample"
            rm -f "$tmp_bam"
            printf "%s\tquickcheck_failed\n" "$sample" > "$status_file"
        fi
    else
        echo "ERROR: mapping failed for $sample"
        rm -f "$tmp_bam"
        printf "%s\tfailed\n" "$sample" > "$status_file"
    fi
}

export OUTBASE TRIMMED_DIR REF MAPPED_DIR LOGS_DIR STATUS_DIR
export BWA_THREADS SORT_THREADS
export -f map_one_sample

parallel --will-cite --jobs "$JOBS" --joblog "$LOGS_DIR/parallel_mapping.joblog" map_one_sample ::: "${samples[@]}"

{
    printf "sample\tstatus\n"
    for sample in "${samples[@]}"; do
        status_file="$STATUS_DIR/${sample}.status.tsv"
        if [[ -s "$status_file" ]]; then
            cat "$status_file"
        else
            printf "%s\tmissing_status_file\n" "$sample"
        fi
    done
} > "$MAP_MANIFEST"

echo
echo "Parallel mapping finished."
echo "Mapping manifest: $MAP_MANIFEST"
echo "Parallel job log: $LOGS_DIR/parallel_mapping.joblog"
echo
echo "Status counts:"
tail -n +2 "$MAP_MANIFEST" | cut -f2 | sort | uniq -c
echo
echo "Completed BAM files:"
find "$MAPPED_DIR" -maxdepth 1 -type f -name '*.bam' | wc -l
echo "Completed BAM indexes:"
find "$MAPPED_DIR" -maxdepth 1 -type f -name '*.bam.bai' | wc -l
echo "Completed flagstat reports:"
find "$LOGS_DIR" -maxdepth 1 -type f -name '*.flagstat.txt' | wc -l
