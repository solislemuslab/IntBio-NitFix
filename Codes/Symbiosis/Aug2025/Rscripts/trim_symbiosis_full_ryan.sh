#!/bin/bash

set -uo pipefail

module load SolisLemus-BioPhylo/2026.04.20

BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"
RAW="$BASE/raw_symbiosis_full"
TRIMMED_DIR="$BASE/symbiosis_trimmed_ryan_full"
LOGS_DIR="$BASE/fastp_logs_ryan_full"
MANIFEST="$BASE/symbiosis_trimmed_ryan_full_manifest.tsv"

mkdir -p "$TRIMMED_DIR" "$LOGS_DIR"

if ! command -v fastp >/dev/null 2>&1; then
    echo "ERROR: fastp is not available."
    exit 1
fi

if [[ ! -f "$MANIFEST" ]]; then
    printf "sample\tstatus\n" > "$MANIFEST"
fi

shopt -s nullglob
r1_files=( "$RAW"/*_R1.fq.gz )

echo "Total R1 files found: ${#r1_files[@]}"
echo "Input folder:  $RAW"
echo "Output folder: $TRIMMED_DIR"
echo "BLAN control samples are retained for later QC."
echo

for r1 in "${r1_files[@]}"; do
    sample=$(basename "$r1" _R1.fq.gz)
    r2="$RAW/${sample}_R2.fq.gz"

    p1="$TRIMMED_DIR/${sample}_P1.fastq.gz"
    p2="$TRIMMED_DIR/${sample}_P2.fastq.gz"
    u1="$TRIMMED_DIR/${sample}_U1.fastq.gz"
    u2="$TRIMMED_DIR/${sample}_U2.fastq.gz"
    html="$LOGS_DIR/${sample}.fastp.html"
    json="$LOGS_DIR/${sample}.fastp.json"
    err="$LOGS_DIR/${sample}.fastp.stderr.log"

    if [[ ! -f "$r2" ]]; then
        echo "ERROR: Missing R2 for $sample"
        printf "%s\tmissing_R2\n" "$sample" >> "$MANIFEST"
        continue
    fi

    if [[ -s "$p1" && -s "$p2" && -s "$json" ]]; then
        echo "Already completed, skipping: $sample"
        continue
    fi

    echo "Processing sample: $sample"

    if fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$p1" \
        -O "$p2" \
        --unpaired1 "$u1" \
        --unpaired2 "$u2" \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --trim_poly_g \
        --poly_g_min_len 5 \
        --trim_poly_x \
        --poly_x_min_len 5 \
        --low_complexity_filter \
        --complexity_threshold 20 \
        --cut_right \
        --cut_right_window_size 4 \
        --cut_right_mean_quality 20 \
        --length_required 36 \
        --thread 8 \
        --compression 4 \
        --html "$html" \
        --json "$json" \
        2> "$err"; then

        echo "Completed: $sample"
        printf "%s\tcompleted\n" "$sample" >> "$MANIFEST"
    else
        echo "ERROR: fastp failed for $sample. See $err"
        printf "%s\tfailed\n" "$sample" >> "$MANIFEST"
    fi
done

echo
echo "Trimming run finished."
echo "Trimmed reads: $TRIMMED_DIR"
echo "Logs:          $LOGS_DIR"
echo "Manifest:      $MANIFEST"

echo "Completed output pairs:"
find "$TRIMMED_DIR" -maxdepth 1 -type f -name '*_P1.fastq.gz' | wc -l
