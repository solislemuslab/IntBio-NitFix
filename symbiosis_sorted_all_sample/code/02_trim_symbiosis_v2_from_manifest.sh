#!/bin/bash

set -uo pipefail

module load SolisLemus-BioPhylo/2026.04.20

OUTBASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
MANIFEST="$OUTBASE/metadata/symbiosis_v2_clean_manifest.tsv"
TRIMMED_DIR="$OUTBASE/symbiosis_trimmed_fastp"
LOGS_DIR="$OUTBASE/fastp_logs"
TRIM_MANIFEST="$OUTBASE/metadata/symbiosis_v2_trim_manifest.tsv"

mkdir -p "$TRIMMED_DIR" "$LOGS_DIR" "$OUTBASE/metadata"

if ! command -v fastp >/dev/null 2>&1; then
    echo "ERROR: fastp is not available."
    exit 1
fi

if [[ ! -s "$MANIFEST" ]]; then
    echo "ERROR: clean manifest not found:"
    echo "$MANIFEST"
    echo "Run Rscripts_v2/01_create_symbiosis_v2_manifest.sh first."
    exit 1
fi

printf "sample\tsequencing_folder\tsite\tsample_type\tstatus\n" > "$TRIM_MANIFEST"

echo "Trimming V2 samples from manifest:"
echo "$MANIFEST"
echo "Output folder:"
echo "$TRIMMED_DIR"
echo

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r folder site type sample has_R2 R1 R2 reason; do
    p1="$TRIMMED_DIR/${sample}_P1.fastq.gz"
    p2="$TRIMMED_DIR/${sample}_P2.fastq.gz"
    u1="$TRIMMED_DIR/${sample}_U1.fastq.gz"
    u2="$TRIMMED_DIR/${sample}_U2.fastq.gz"
    html="$LOGS_DIR/${sample}.fastp.html"
    json="$LOGS_DIR/${sample}.fastp.json"
    err="$LOGS_DIR/${sample}.fastp.stderr.log"

    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        echo "ERROR: missing R1/R2 for $sample"
        printf "%s\t%s\t%s\t%s\tmissing_input\n" "$sample" "$folder" "$site" "$type" >> "$TRIM_MANIFEST"
        continue
    fi

    if [[ -s "$p1" && -s "$p2" && -s "$json" ]]; then
        echo "Already trimmed, skipping: $sample"
        printf "%s\t%s\t%s\t%s\talready_completed\n" "$sample" "$folder" "$site" "$type" >> "$TRIM_MANIFEST"
        continue
    fi

    echo "Trimming: $sample"

    if fastp \
        -i "$R1" \
        -I "$R2" \
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
        printf "%s\t%s\t%s\t%s\tcompleted\n" "$sample" "$folder" "$site" "$type" >> "$TRIM_MANIFEST"
    else
        echo "ERROR: fastp failed for $sample. See $err"
        printf "%s\t%s\t%s\t%s\tfailed\n" "$sample" "$folder" "$site" "$type" >> "$TRIM_MANIFEST"
    fi
done

echo
echo "Trimming finished."
echo "Trim manifest: $TRIM_MANIFEST"
echo "Completed P1 files:"
find "$TRIMMED_DIR" -maxdepth 1 -type f -name '*_P1.fastq.gz' | wc -l
