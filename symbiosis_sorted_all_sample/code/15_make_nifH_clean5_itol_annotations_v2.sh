#!/usr/bin/env bash
set -euo pipefail

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2}"
META="$OUT/nifH_clean5_tree_ready_v2/nifH_clean5_tree_ready_metadata.tsv"
ITOL="$OUT/nifH_clean5_itol_annotations_v2"

mkdir -p "$ITOL"

if [[ ! -s "$META" ]]; then
  echo "ERROR: metadata file does not exist or is empty: $META"
  exit 1
fi

echo "Step 15 V2: create iTOL annotation files for nifH clean5 trees"
echo "Metadata: $META"
echo "Output:   $ITOL"
echo

python3 - "$META" "$ITOL" <<'PY'
from pathlib import Path
import csv
import sys

meta = Path(sys.argv[1])
itol = Path(sys.argv[2])
itol.mkdir(parents=True, exist_ok=True)

rows = []
with open(meta, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        rows.append(row)

sample_type_colors = {
    "No": "#2E86AB",
    "Rh": "#F18F01",
    "Ro": "#7A5195",
}

target_ref_colors = {
    "ref52": "#1B9E77",
    "ref56": "#D95F02",
    "ref60": "#7570B3",
    "ref62": "#E7298A",
    "ref63": "#66A61E",
}

status_colors = {
    "pass_single_dominant": "#20854E",
    "pass_mixed_possible_multicopy": "#C44E52",
}

def write_colorstrip(path, label, legend_title, legend_labels, legend_colors, data_rows):
    with open(path, "w") as out:
        out.write("DATASET_COLORSTRIP\n")
        out.write("SEPARATOR TAB\n")
        out.write(f"DATASET_LABEL\t{label}\n")
        out.write("COLOR\t#000000\n")
        out.write(f"LEGEND_TITLE\t{legend_title}\n")
        out.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_labels)) + "\n")
        out.write("LEGEND_COLORS\t" + "\t".join(legend_colors) + "\n")
        out.write("LEGEND_LABELS\t" + "\t".join(legend_labels) + "\n")
        out.write("DATA\n")
        for tree_id, color, value in data_rows:
            out.write(f"{tree_id}\t{color}\t{value}\n")

def rows_for(seq_set, field, colors):
    data = []
    for row in rows:
        if row["sequence_set"] != seq_set:
            continue
        value = row[field]
        color = colors.get(value, "#999999")
        data.append((row["tree_id"], color, value))
    return data

outputs = [
    (
        "strict_itol_sample_type_colorstrip.txt",
        "strict_single_dominant",
        "Sample type",
        "Sample type",
        ["No", "Rh", "Ro"],
        [sample_type_colors["No"], sample_type_colors["Rh"], sample_type_colors["Ro"]],
        rows_for("strict_single_dominant", "sample_type", sample_type_colors),
    ),
    (
        "strict_itol_target_ref_colorstrip.txt",
        "strict_single_dominant",
        "Target ref",
        "Target ref",
        ["ref52", "ref56", "ref60", "ref62", "ref63"],
        [target_ref_colors[x] for x in ["ref52", "ref56", "ref60", "ref62", "ref63"]],
        rows_for("strict_single_dominant", "target_ref", target_ref_colors),
    ),
    (
        "iupac_itol_sample_type_colorstrip.txt",
        "iupac_all_pass",
        "Sample type",
        "Sample type",
        ["No", "Rh", "Ro"],
        [sample_type_colors["No"], sample_type_colors["Rh"], sample_type_colors["Ro"]],
        rows_for("iupac_all_pass", "sample_type", sample_type_colors),
    ),
    (
        "iupac_itol_target_ref_colorstrip.txt",
        "iupac_all_pass",
        "Target ref",
        "Target ref",
        ["ref52", "ref56", "ref60", "ref62", "ref63"],
        [target_ref_colors[x] for x in ["ref52", "ref56", "ref60", "ref62", "ref63"]],
        rows_for("iupac_all_pass", "target_ref", target_ref_colors),
    ),
    (
        "iupac_itol_consensus_status_colorstrip.txt",
        "iupac_all_pass",
        "Consensus status",
        "Consensus status",
        ["single_dominant", "mixed_possible_multicopy"],
        [status_colors["pass_single_dominant"], status_colors["pass_mixed_possible_multicopy"]],
        rows_for("iupac_all_pass", "status", status_colors),
    ),
]

for filename, seq_set, label, legend_title, legend_labels, legend_colors, data_rows in outputs:
    write_colorstrip(itol / filename, label, legend_title, legend_labels, legend_colors, data_rows)

print("iTOL files written:")
for path in sorted(itol.glob("*_itol_*_colorstrip.txt")):
    print(f"{path}\t{path.stat().st_size} bytes")

print()
print("Data row counts:")
for path in sorted(itol.glob("*_itol_*_colorstrip.txt")):
    in_data = False
    count = 0
    with open(path) as handle:
        for line in handle:
            if line.strip() == "DATA":
                in_data = True
                continue
            if in_data and line.strip():
                count += 1
    print(f"{path.name}\t{count}")
PY
