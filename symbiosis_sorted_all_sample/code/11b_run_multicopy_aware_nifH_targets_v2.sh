#!/bin/bash

set -uo pipefail

OUTBASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
RUNNER="$OUTBASE/Rscripts_v2/11_extract_multicopy_aware_target_sequences_v2.sh"
LOGDIR="$OUTBASE/multicopy_aware_target_sequences_v2/run_logs"

mkdir -p "$LOGDIR"

if [[ ! -s "$RUNNER" ]]; then
    echo "ERROR: target extraction runner not found:"
    echo "$RUNNER"
    exit 1
fi

# Default = MC-clean nifH targets from Step 10 strict 80/10 ranking.
# Set NIFH_TARGET_SET=all to run all 14 nifH references.
NIFH_TARGET_SET="${NIFH_TARGET_SET:-clean}"

clean_targets=(
    'nifH|NC_009937|NC_009937_-_Symbiosis_Island_4|ref56'
    'nifH|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island|ref63'
    'nifH|NC_002678|NC_002678_-_Symbiosis_Island_3|ref52'
    'nifH|NC_015656|NC_015656_-_Symbiosis_Island|ref60'
    'nifH|NC_019848|NC_019848_-_Symbiosis_Island_2|ref62'
)

all_targets=(
    'nifH|CU633751|CU633751_-_Symbiosis_Island|ref51'
    'nifH|NC_002678|NC_002678_-_Symbiosis_Island_3|ref52'
    'nifH|NC_007777|NC_007777_-_Symbiosis_Island|ref53'
    'nifH|NC_008278|NC_008278_-_Symbiosis_Island|ref54'
    'nifH|NC_009937|NC_009937_-_Symbiosis_Island_3|ref55'
    'nifH|NC_009937|NC_009937_-_Symbiosis_Island_4|ref56'
    'nifH|NC_011894|NC_011894_-_Symbiosis_Island|ref57'
    'nifH|NC_012848|NC_012848_-_Symbiosis_Island_2|ref58'
    'nifH|NC_014120|NC_014120_-_Symbiosis_Island_3|ref59'
    'nifH|NC_015656|NC_015656_-_Symbiosis_Island|ref60'
    'nifH|NC_017249|NC_017249_-_Symbiosis_Island|ref61'
    'nifH|NC_019848|NC_019848_-_Symbiosis_Island_2|ref62'
    'nifH|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island|ref63'
    'nifH|NZ_HG938357|NZ_HG938357_-_Symbiosis_Island_2|ref64'
)

if [[ "$NIFH_TARGET_SET" == "all" ]]; then
    targets=( "${all_targets[@]}" )
elif [[ "$NIFH_TARGET_SET" == "clean" ]]; then
    targets=( "${clean_targets[@]}" )
else
    echo "ERROR: NIFH_TARGET_SET must be clean or all."
    exit 1
fi

echo "Running multicopy-aware nifH extraction"
echo "Target set: $NIFH_TARGET_SET"
echo "Targets: ${#targets[@]}"
echo "Sample group: ${SAMPLE_GROUP:-all}"
echo "MIN_BASE_DEPTH: ${MIN_BASE_DEPTH:-1}"
echo "MAX_N_PERCENT: ${MAX_N_PERCENT:-20}"
echo "MAX_SAMPLES: ${MAX_SAMPLES:-0}"
echo

for target in "${targets[@]}"; do
    ref_label="${target##*|}"
    log="$LOGDIR/nifH_${ref_label}_${NIFH_TARGET_SET}_run.log"

    echo "============================================================"
    echo "Target: $target"
    echo "Log:    $log"
    echo "============================================================"

    SAMPLE_GROUP="${SAMPLE_GROUP:-all}" \
    MIN_BASE_DEPTH="${MIN_BASE_DEPTH:-1}" \
    MAX_N_PERCENT="${MAX_N_PERCENT:-20}" \
    MAX_SAMPLES="${MAX_SAMPLES:-0}" \
    TARGET_ID="$target" \
    bash "$RUNNER" | tee "$log"
done

echo
echo "All requested nifH target extractions finished."
echo "Logs:"
echo "$LOGDIR"
