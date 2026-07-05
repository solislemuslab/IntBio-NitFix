#!/bin/bash

set -uo pipefail

module load SolisLemus-BioPhylo/2026.04.20
module load samtools-1.9

OUTBASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"

# Default target is the strongest MC-aware nifH target from Step 10.
TARGET_ID="${TARGET_ID:-nifH|NC_009937|NC_009937_-_Symbiosis_Island_4|ref56}"
SAMPLE_GROUP="${SAMPLE_GROUP:-No}"

# Sample-selection threshold from Step 9 coverage table.
MIN_PERCENT_COVERED="${MIN_PERCENT_COVERED:-80}"
MIN_MEAN_DEPTH="${MIN_MEAN_DEPTH:-10}"

# Per-base consensus and mixed-signal settings.
MIN_BASE_DEPTH="${MIN_BASE_DEPTH:-5}"
MAX_N_PERCENT="${MAX_N_PERCENT:-5}"
MIXED_MIN_DEPTH="${MIXED_MIN_DEPTH:-20}"
MIXED_MINOR_COUNT="${MIXED_MINOR_COUNT:-5}"
MIXED_MINOR_AF="${MIXED_MINOR_AF:-0.20}"
MAX_MIXED_POSITIONS_FOR_SINGLE="${MAX_MIXED_POSITIONS_FOR_SINGLE:-10}"
MAX_SAMPLES="${MAX_SAMPLES:-0}"

export BASE="$OUTBASE"
export TARGET_ID SAMPLE_GROUP
export MIN_PERCENT_COVERED MIN_MEAN_DEPTH
export MIN_BASE_DEPTH MAX_N_PERCENT
export MIXED_MIN_DEPTH MIXED_MINOR_COUNT MIXED_MINOR_AF MAX_MIXED_POSITIONS_FOR_SINGLE
export MAX_SAMPLES

python3 "$OUTBASE/Rscripts_v2/11_extract_multicopy_aware_target_sequences_v2.py"
