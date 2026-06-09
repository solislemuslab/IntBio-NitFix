#!/bin/bash

set -uo pipefail

module load SolisLemus-BioPhylo/2026.04.20

BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"

BAM_DIR="$BASE/symbiosis_mapped_full"
REGIONS="$BASE/nif_nod_original_reference_regions/nif_nod_matches_in_original_reference.tsv"

OUTDIR="$BASE/nif_nod_coverage_existing_mapping"
RAW_OUT="$OUTDIR/nif_nod_region_coverage_all_samples.tsv"
LOG="$OUTDIR/nif_nod_region_coverage_run.log"

mkdir -p "$OUTDIR"

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available."
    exit 1
fi

if [[ ! -s "$REGIONS" ]]; then
    echo "ERROR: Region table not found:"
    echo "$REGIONS"
    exit 1
fi

shopt -s nullglob
BAMS=("$BAM_DIR"/*.bam)

if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found in:"
    echo "$BAM_DIR"
    exit 1
fi

echo "BAM files found: ${#BAMS[@]}"
echo "Region table:    $REGIONS"
echo "Output table:    $RAW_OUT"
echo

echo -e "sample\tsample_type\tblank_control\tgene\toriginal_reference\tbed_start\tbed_end\ttarget_id\tstrand\ttarget_length\tcovered_bases\tpercent_covered\tmean_depth\tmax_depth" > "$RAW_OUT"

for BAM in "${BAMS[@]}"; do

    SAMPLE=$(basename "$BAM" .bam)

    BLANK="no"
    if [[ "$SAMPLE" == BLAN* ]]; then
        BLANK="yes"
    fi

    if [[ "$SAMPLE" == *-No ]]; then
        TYPE="No"
    elif [[ "$SAMPLE" == *-Rh ]]; then
        TYPE="Rh"
    elif [[ "$SAMPLE" == *-Ro ]]; then
        TYPE="Ro"
    else
        TYPE="other"
    fi

    echo "Calculating coverage: $SAMPLE"

    tail -n +2 "$REGIONS" | while IFS=$'\t' read -r REF START END GENE TARGET STRAND LENGTH; do

        REGION="${REF}:$((START + 1))-${END}"

        STATS=$(samtools depth -a -r "$REGION" "$BAM" \
            | awk -v len="$LENGTH" '
                {
                    positions++;
                    sum += $3;
                    if ($3 > 0) covered++;
                    if ($3 > max) max = $3;
                }
                END {
                    if (positions == 0) {
                        printf "0\t0.0000\t0.0000\t0";
                    } else {
                        printf "%d\t%.4f\t%.4f\t%d",
                            covered,
                            (covered / len) * 100,
                            sum / len,
                            max;
                    }
                }'
        )

        echo -e "${SAMPLE}\t${TYPE}\t${BLANK}\t${GENE}\t${REF}\t${START}\t${END}\t${TARGET}\t${STRAND}\t${LENGTH}\t${STATS}" >> "$RAW_OUT"

    done

done

echo
echo "Coverage analysis complete."
echo "Output table:"
echo "$RAW_OUT"
