#!/bin/bash

set -uo pipefail

# V2 manifest builder.
# Reads original FASTQ files but writes only to processed-data/symbiosis_sorted_v2.

ROOT="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/original-data"
OUTBASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
AUGUST_PROCESSED="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"

MANIFEST_ALL="$OUTBASE/metadata/symbiosis_v2_all_R1_inventory.tsv"
MANIFEST_CLEAN="$OUTBASE/metadata/symbiosis_v2_clean_manifest.tsv"
DUPLICATES="$OUTBASE/metadata/symbiosis_v2_duplicate_samples.tsv"
EXCLUDED="$OUTBASE/metadata/symbiosis_v2_excluded_files.tsv"
COUNTS_FOLDER="$OUTBASE/metadata/symbiosis_v2_counts_by_folder.tsv"
COUNTS_SITE="$OUTBASE/metadata/symbiosis_v2_counts_by_site.tsv"
COUNTS_FOLDER_SITE_TYPE="$OUTBASE/metadata/symbiosis_v2_counts_by_folder_site_type.tsv"

mkdir -p "$OUTBASE/metadata" "$OUTBASE/reference"

if [[ ! -d "$ROOT" ]]; then
    echo "ERROR: original data root not found:"
    echo "$ROOT"
    exit 1
fi

echo "Building V2 symbiosis_sorted manifest"
echo "Original data root: $ROOT"
echo "V2 output folder:   $OUTBASE"
echo

printf "sequencing_folder\tsite\tsample_type\tsample\thas_R2\tR1_path\tR2_path\texclusion_reason\n" > "$MANIFEST_ALL"

find "$ROOT" -path "*/symbiosis_sorted/*_R1.fq.gz" -type f \
| sort \
| while read -r R1; do

    rel="${R1#$ROOT/}"
    folder=$(echo "$rel" | cut -d'/' -f1)
    sample=$(basename "$R1" _R1.fq.gz)
    site=$(echo "$sample" | cut -d'-' -f1)

    if [[ "$sample" == *-No ]]; then
        type="No"
    elif [[ "$sample" == *-Rh ]]; then
        type="Rh"
    elif [[ "$sample" == *-Ro ]]; then
        type="Ro"
    elif [[ "$sample" == MC-* ]]; then
        type="MC"
        site="MC"
    else
        type="other"
    fi

    R2="${R1/_R1.fq.gz/_R2.fq.gz}"
    if [[ -s "$R2" ]]; then
        has_R2="yes"
    else
        has_R2="no"
    fi

    reason="keep_candidate"
    if [[ "$R1" == *"/AAA_delete/"* ]]; then
        reason="exclude_AAA_delete_duplicate_copy"
    elif [[ "$has_R2" == "no" ]]; then
        reason="exclude_missing_R2"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$folder" "$site" "$type" "$sample" "$has_R2" "$R1" "$R2" "$reason"
done >> "$MANIFEST_ALL"

# Keep one row per sample name, excluding known duplicate-copy folders and missing R2.
# If the same sample appears in multiple sequencing folders, prefer December2025 over
# August2025 over May2025 over July2024, and record the excluded duplicate rows.
python3 - "$MANIFEST_ALL" "$MANIFEST_CLEAN" "$DUPLICATES" "$EXCLUDED" <<'PY'
import csv
import sys
from collections import defaultdict

manifest_all, manifest_clean, duplicates, excluded = sys.argv[1:5]

priority = {
    "intbio_reads_NovogeneDecember2025": 4,
    "intbio_reads_NovogeneAugust2025": 3,
    "intbio_reads_NovogeneMay2025": 2,
    "intbio_reads_NovogeneJuly2024": 1,
}

with open(manifest_all, newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

candidate_rows = []
excluded_rows = []
for row in rows:
    if row["exclusion_reason"] == "keep_candidate":
        candidate_rows.append(row)
    else:
        excluded_rows.append(row)

by_sample = defaultdict(list)
for row in candidate_rows:
    by_sample[row["sample"]].append(row)

clean_rows = []
duplicate_rows = []
for sample, sample_rows in sorted(by_sample.items()):
    sample_rows = sorted(
        sample_rows,
        key=lambda r: (priority.get(r["sequencing_folder"], 0), r["R1_path"]),
        reverse=True,
    )
    keep = sample_rows[0].copy()
    keep["exclusion_reason"] = "keep"
    clean_rows.append(keep)

    if len(sample_rows) > 1:
        for row in sample_rows:
            out = row.copy()
            out["exclusion_reason"] = "duplicate_kept" if row is sample_rows[0] else "duplicate_excluded"
            duplicate_rows.append(out)
            if row is not sample_rows[0]:
                excluded_rows.append(out)

fieldnames = [
    "sequencing_folder", "site", "sample_type", "sample",
    "has_R2", "R1_path", "R2_path", "exclusion_reason",
]

for path, out_rows in [
    (manifest_clean, clean_rows),
    (duplicates, duplicate_rows),
    (excluded, excluded_rows),
]:
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)
PY

awk -F'\t' '
    NR > 1 {count[$1]++}
    END {
        print "sequencing_folder\tcount"
        for (k in count) print k "\t" count[k]
    }
' "$MANIFEST_CLEAN" | sort > "$COUNTS_FOLDER"

awk -F'\t' '
    NR > 1 {count[$2]++}
    END {
        print "site\tcount"
        for (k in count) print k "\t" count[k]
    }
' "$MANIFEST_CLEAN" | sort > "$COUNTS_SITE"

awk -F'\t' '
    NR > 1 {count[$1"\t"$2"\t"$3]++}
    END {
        print "sequencing_folder\tsite\tsample_type\tcount"
        for (k in count) print k "\t" count[k]
    }
' "$MANIFEST_CLEAN" | sort > "$COUNTS_FOLDER_SITE_TYPE"

# Copy reference files into the V2 workspace so BWA indexes are written in V2,
# not into the previous August2025 analysis folder.
if [[ -s "$AUGUST_PROCESSED/symbiosis_islands.fasta" ]]; then
    cp -n "$AUGUST_PROCESSED/symbiosis_islands.fasta" "$OUTBASE/reference/symbiosis_islands.fasta"
fi
if [[ -s "$AUGUST_PROCESSED/symbiosis_islands.gb" ]]; then
    cp -n "$AUGUST_PROCESSED/symbiosis_islands.gb" "$OUTBASE/reference/symbiosis_islands.gb"
fi
if [[ -s "$AUGUST_PROCESSED/symbiosis_islands_gene_list.xlsx" ]]; then
    cp -n "$AUGUST_PROCESSED/symbiosis_islands_gene_list.xlsx" "$OUTBASE/reference/symbiosis_islands_gene_list.xlsx"
fi

echo
echo "Manifest complete."
echo "All inventory:       $MANIFEST_ALL"
echo "Clean V2 manifest:   $MANIFEST_CLEAN"
echo "Duplicate samples:   $DUPLICATES"
echo "Excluded files:      $EXCLUDED"
echo
echo "Clean sample count:"
awk 'NR > 1 {count++} END {print count}' "$MANIFEST_CLEAN"
echo
echo "Counts by sequencing folder:"
cat "$COUNTS_FOLDER"
echo
echo "Counts by sample type:"
awk -F'\t' 'NR > 1 {count[$3]++} END {for (k in count) print k "\t" count[k]}' "$MANIFEST_CLEAN" | sort
echo
echo "MC controls:"
awk -F'\t' 'NR == 1 || $3 == "MC"' "$MANIFEST_CLEAN"
