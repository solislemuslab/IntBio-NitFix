#!/bin/bash

set -uo pipefail

module load SolisLemus-BioPhylo/2026.04.20
module load samtools-1.9

OUTBASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
BAM_DIR="$OUTBASE/symbiosis_mapped_full"
REGIONS="$OUTBASE/nif_nod_original_reference_regions/nif_nod_matches_in_original_reference.tsv"
MANIFEST="$OUTBASE/metadata/symbiosis_v2_clean_manifest.tsv"
SCRIPT="$OUTBASE/Rscripts_v2/09_calculate_one_sample_nif_nod_coverage_v2.py"

OUTDIR="$OUTBASE/nif_nod_coverage_existing_mapping_v2"
TMPDIR="$OUTDIR/per_sample_coverage"
RAW_OUT="$OUTDIR/nif_nod_region_coverage_all_samples.tsv"
RUN_LOG="$OUTDIR/nif_nod_region_coverage_run.log"

JOBS="${JOBS:-4}"
MAX_SAMPLES="${MAX_SAMPLES:-0}"

mkdir -p "$OUTDIR" "$TMPDIR"

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available."
    exit 1
fi

if ! command -v parallel >/dev/null 2>&1; then
    echo "ERROR: GNU parallel is not available."
    exit 1
fi

if [[ ! -s "$REGIONS" ]]; then
    echo "ERROR: region table not found:"
    echo "$REGIONS"
    exit 1
fi

if [[ ! -s "$MANIFEST" ]]; then
    echo "ERROR: V2 manifest not found:"
    echo "$MANIFEST"
    exit 1
fi

if [[ ! -s "$SCRIPT" ]]; then
    echo "ERROR: per-sample coverage script not found:"
    echo "$SCRIPT"
    exit 1
fi

shopt -s nullglob
bams=( "$BAM_DIR"/*.bam )

if [[ "${#bams[@]}" -eq 0 ]]; then
    echo "ERROR: no BAM files found in:"
    echo "$BAM_DIR"
    exit 1
fi

if [[ "$MAX_SAMPLES" -gt 0 ]]; then
    bams=( "${bams[@]:0:$MAX_SAMPLES}" )
fi

echo "Calculating V2 nif/nod region coverage"
echo "BAM files selected: ${#bams[@]}"
echo "Region table:       $REGIONS"
echo "Output table:       $RAW_OUT"
echo "GNU parallel jobs:  $JOBS"
echo "Max samples:        $MAX_SAMPLES"
echo

calculate_one() {
    bam="$1"
    sample=$(basename "$bam" .bam)
    out="$TMPDIR/${sample}.coverage.tsv"

    if [[ -s "$out" ]]; then
        echo "Already completed, skipping: $sample"
        return 0
    fi

    echo "Coverage: $sample"
    python3 "$SCRIPT" \
        --bam "$bam" \
        --regions "$REGIONS" \
        --manifest "$MANIFEST" \
        --out "$out"
}

export TMPDIR SCRIPT REGIONS MANIFEST
export -f calculate_one

parallel --will-cite --jobs "$JOBS" --joblog "$OUTDIR/parallel_coverage.joblog" calculate_one ::: "${bams[@]}" \
    | tee "$RUN_LOG"

{
    echo -e "sample\tsequencing_folder\tsite\tsample_type\tis_mc_control\tgene\toriginal_reference\tbed_start\tbed_end\ttarget_id\tstrand\ttarget_length\tcovered_bases\tpercent_covered\tmean_depth\tmax_depth"
    for bam in "${bams[@]}"; do
        sample=$(basename "$bam" .bam)
        file="$TMPDIR/${sample}.coverage.tsv"
        if [[ -s "$file" ]]; then
            cat "$file"
        else
            echo "WARNING: missing coverage output for $sample" >&2
        fi
    done
} > "$RAW_OUT"

echo
echo "Coverage analysis complete."
echo "Output table:"
echo "$RAW_OUT"
echo
echo "Rows including header:"
wc -l "$RAW_OUT"
