#!/bin/bash

set -uo pipefail

OUTBASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
MAPPED_DIR="$OUTBASE/symbiosis_mapped_full"
LOGS_DIR="$OUTBASE/symbiosis_mapping_logs"
MANIFEST="$OUTBASE/metadata/symbiosis_v2_clean_manifest.tsv"
MAP_MANIFEST="$OUTBASE/metadata/symbiosis_v2_mapping_manifest.tsv"
QC_DIR="$OUTBASE/symbiosis_mapping_QC_checking"

mkdir -p "$QC_DIR"

COMPLETION="$QC_DIR/symbiosis_mapping_completion_summary.txt"
FLAGSTAT_TABLE="$QC_DIR/symbiosis_mapping_flagstat_summary.tsv"
OVERALL="$QC_DIR/symbiosis_mapping_overall_summary.tsv"
BY_TYPE="$QC_DIR/symbiosis_mapping_by_sample_type.tsv"
BY_FOLDER="$QC_DIR/symbiosis_mapping_by_sequencing_folder.tsv"
FIGURE="$QC_DIR/symbiosis_mapping_qc_histograms.svg"

echo "Checking V2 symbiosis mapping quality"
echo "V2 output folder: $OUTBASE"
echo "Mapped BAM folder: $MAPPED_DIR"
echo "Mapping log folder: $LOGS_DIR"
echo "QC output folder: $QC_DIR"
echo

if [[ ! -s "$MANIFEST" ]]; then
    echo "ERROR: clean manifest not found: $MANIFEST"
    exit 1
fi

if [[ ! -d "$MAPPED_DIR" ]]; then
    echo "ERROR: mapped BAM folder not found: $MAPPED_DIR"
    exit 1
fi

if [[ ! -d "$LOGS_DIR" ]]; then
    echo "ERROR: mapping log folder not found: $LOGS_DIR"
    exit 1
fi

{
    echo "Mapping completion summary"
    echo "=========================="
    echo
    echo "Expected samples from manifest:"
    tail -n +2 "$MANIFEST" | wc -l
    echo
    echo "BAM files:"
    find "$MAPPED_DIR" -maxdepth 1 -type f -name '*.bam' | wc -l
    echo
    echo "BAM index files:"
    find "$MAPPED_DIR" -maxdepth 1 -type f -name '*.bam.bai' | wc -l
    echo
    echo "Flagstat reports:"
    find "$LOGS_DIR" -maxdepth 1 -type f -name '*.flagstat.txt' | wc -l
    echo
    echo "Idxstats reports:"
    find "$LOGS_DIR" -maxdepth 1 -type f -name '*.idxstats.txt' | wc -l
    echo
    echo "Zero-size BAM files:"
    find "$MAPPED_DIR" -maxdepth 1 -type f -name '*.bam' -size 0 -print | wc -l
    echo
    echo "Zero-size BAM index files:"
    find "$MAPPED_DIR" -maxdepth 1 -type f -name '*.bam.bai' -size 0 -print | wc -l
    if [[ -s "$MAP_MANIFEST" ]]; then
        echo
        echo "Mapping manifest status counts:"
        tail -n +2 "$MAP_MANIFEST" | cut -f2 | sort | uniq -c
    fi
} | tee "$COMPLETION"

python3 - "$OUTBASE" <<'PY'
import csv
import html
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

outbase = Path(sys.argv[1])
manifest_path = outbase / "metadata" / "symbiosis_v2_clean_manifest.tsv"
logs_dir = outbase / "symbiosis_mapping_logs"
qc_dir = outbase / "symbiosis_mapping_QC_checking"

flagstat_table = qc_dir / "symbiosis_mapping_flagstat_summary.tsv"
overall_path = qc_dir / "symbiosis_mapping_overall_summary.tsv"
by_type_path = qc_dir / "symbiosis_mapping_by_sample_type.tsv"
by_folder_path = qc_dir / "symbiosis_mapping_by_sequencing_folder.tsv"
figure_path = qc_dir / "symbiosis_mapping_qc_histograms.svg"

def mean(values):
    values = [v for v in values if v is not None]
    return statistics.mean(values) if values else None

def median(values):
    values = sorted(v for v in values if v is not None)
    return statistics.median(values) if values else None

def fmt(value, digits=2):
    if value is None:
        return "NA"
    return f"{value:.{digits}f}"

def read_manifest(path):
    rows = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["sample"]] = row
    return rows

def parse_flagstat(path):
    text = path.read_text()
    out = {
        "total_reads": None,
        "mapped_reads": None,
        "mapped_percent": None,
        "properly_paired_reads": None,
        "properly_paired_percent": None,
        "singletons": None,
        "singletons_percent": None,
    }
    for line in text.splitlines():
        if " in total " in line:
            out["total_reads"] = int(line.split()[0])
        elif " mapped (" in line and "mate" not in line:
            out["mapped_reads"] = int(line.split()[0])
            match = re.search(r"\(([^%]+)%", line)
            if match:
                out["mapped_percent"] = float(match.group(1))
        elif " properly paired " in line:
            out["properly_paired_reads"] = int(line.split()[0])
            match = re.search(r"\(([^%]+)%", line)
            if match:
                out["properly_paired_percent"] = float(match.group(1))
        elif " singletons " in line:
            out["singletons"] = int(line.split()[0])
            match = re.search(r"\(([^%]+)%", line)
            if match:
                out["singletons_percent"] = float(match.group(1))
    return out

metadata = read_manifest(manifest_path)
rows = []

for sample, meta in sorted(metadata.items()):
    path = logs_dir / f"{sample}.flagstat.txt"
    metrics = {
        "total_reads": None,
        "mapped_reads": None,
        "mapped_percent": None,
        "properly_paired_reads": None,
        "properly_paired_percent": None,
        "singletons": None,
        "singletons_percent": None,
    }
    status = "missing_flagstat"
    if path.exists():
        try:
            metrics = parse_flagstat(path)
            status = "ok"
        except Exception as exc:
            status = f"parse_error:{type(exc).__name__}"

    rows.append({
        "sample": sample,
        "sequencing_folder": meta.get("sequencing_folder", ""),
        "site": meta.get("site", ""),
        "sample_type": meta.get("sample_type", ""),
        "flagstat_found": "yes" if path.exists() else "no",
        "status": status,
        **metrics,
    })

cols = [
    "sample",
    "sequencing_folder",
    "site",
    "sample_type",
    "flagstat_found",
    "total_reads",
    "mapped_reads",
    "mapped_percent",
    "properly_paired_reads",
    "properly_paired_percent",
    "singletons",
    "singletons_percent",
    "status",
]

with flagstat_table.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(cols)
    for row in rows:
        writer.writerow([
            row["sample"],
            row["sequencing_folder"],
            row["site"],
            row["sample_type"],
            row["flagstat_found"],
            row["total_reads"] if row["total_reads"] is not None else "NA",
            row["mapped_reads"] if row["mapped_reads"] is not None else "NA",
            fmt(row["mapped_percent"]),
            row["properly_paired_reads"] if row["properly_paired_reads"] is not None else "NA",
            fmt(row["properly_paired_percent"]),
            row["singletons"] if row["singletons"] is not None else "NA",
            fmt(row["singletons_percent"]),
            row["status"],
        ])

ok_rows = [r for r in rows if r["status"] == "ok"]

overall = {
    "manifest_samples": len(rows),
    "flagstat_ok": len(ok_rows),
    "flagstat_missing_or_failed": len(rows) - len(ok_rows),
    "mc_control_samples": sum(1 for r in rows if r["sample_type"] == "MC"),
    "mean_mapped_percent": mean([r["mapped_percent"] for r in ok_rows]),
    "median_mapped_percent": median([r["mapped_percent"] for r in ok_rows]),
    "mean_properly_paired_percent": mean([r["properly_paired_percent"] for r in ok_rows]),
    "median_properly_paired_percent": median([r["properly_paired_percent"] for r in ok_rows]),
    "mean_singletons_percent": mean([r["singletons_percent"] for r in ok_rows]),
}

with overall_path.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["metric", "value"])
    for key, value in overall.items():
        writer.writerow([key, fmt(value) if isinstance(value, float) else value])

def grouped_summary(group_key, out_path):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row[group_key]].append(row)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            group_key,
            "samples",
            "flagstat_ok",
            "mean_mapped_percent",
            "median_mapped_percent",
            "mean_properly_paired_percent",
            "median_properly_paired_percent",
        ])
        for group in sorted(grouped):
            group_rows = grouped[group]
            group_ok = [r for r in group_rows if r["status"] == "ok"]
            writer.writerow([
                group,
                len(group_rows),
                len(group_ok),
                fmt(mean([r["mapped_percent"] for r in group_ok])),
                fmt(median([r["mapped_percent"] for r in group_ok])),
                fmt(mean([r["properly_paired_percent"] for r in group_ok])),
                fmt(median([r["properly_paired_percent"] for r in group_ok])),
            ])

grouped_summary("sample_type", by_type_path)
grouped_summary("sequencing_folder", by_folder_path)

def histogram(values, bins):
    counts = [0 for _ in range(len(bins) - 1)]
    for value in values:
        if value is None:
            continue
        for i in range(len(bins) - 1):
            if i == len(bins) - 2:
                if bins[i] <= value <= bins[i + 1]:
                    counts[i] += 1
                    break
            elif bins[i] <= value < bins[i + 1]:
                counts[i] += 1
                break
    return counts

def draw_hist_panel(title, values, mean_value, median_value, color, mean_color, x, y, w, h, y_label_max=None):
    bins = list(range(0, 105, 5))
    counts = histogram(values, bins)
    max_count = y_label_max if y_label_max is not None else max(counts + [1])
    parts = []
    parts.append(f'<text x="{x+w/2}" y="{y-28}" text-anchor="middle" font-size="19" font-family="Arial" font-weight="bold">{html.escape(title)}</text>')
    parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="white" stroke="#172033"/>')
    for frac in [0, 0.5, 1]:
        cy = y + h - frac * h
        label = int(max_count * frac)
        parts.append(f'<line x1="{x}" y1="{cy:.1f}" x2="{x+w}" y2="{cy:.1f}" stroke="#e0e4ea"/>')
        parts.append(f'<text x="{x-12}" y="{cy+4:.1f}" text-anchor="end" font-size="12" font-family="Arial">{label}</text>')
    for tick in [0, 25, 50, 75, 100]:
        tx = x + (tick / 100) * w
        parts.append(f'<line x1="{tx:.1f}" y1="{y+h}" x2="{tx:.1f}" y2="{y+h+5}" stroke="#172033"/>')
        parts.append(f'<text x="{tx:.1f}" y="{y+h+25}" text-anchor="middle" font-size="12" font-family="Arial">{tick}</text>')
    bin_width = w / (len(bins) - 1)
    for i, count in enumerate(counts):
        bar_h = (count / max_count) * h if max_count else 0
        bx = x + i * bin_width + 1
        by = y + h - bar_h
        parts.append(f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bin_width-2:.1f}" height="{bar_h:.1f}" fill="{color}" opacity="0.75"/>')
    if mean_value is not None:
        mx = x + (mean_value / 100) * w
        parts.append(f'<line x1="{mx:.1f}" y1="{y}" x2="{mx:.1f}" y2="{y+h}" stroke="{mean_color}" stroke-width="3"/>')
        parts.append(f'<text x="{mx:.1f}" y="{y-8}" text-anchor="middle" font-size="12" font-family="Arial" fill="{mean_color}">mean {fmt(mean_value, 1)}%</text>')
    parts.append(f'<text x="{x+w/2}" y="{y+h+55}" text-anchor="middle" font-size="14" font-family="Arial">Percent</text>')
    parts.append(f'<text x="{x-58}" y="{y+h/2}" transform="rotate(-90 {x-58},{y+h/2})" text-anchor="middle" font-size="14" font-family="Arial">Number of samples</text>')
    parts.append(f'<text x="{x+w/2}" y="{y+h+82}" text-anchor="middle" font-size="13" font-family="Arial">median {fmt(median_value, 1)}%</text>')
    return parts, max(counts + [1])

mapped_values = [r["mapped_percent"] for r in ok_rows]
paired_values = [r["properly_paired_percent"] for r in ok_rows]
mapped_counts = histogram(mapped_values, list(range(0, 105, 5)))
paired_counts = histogram(paired_values, list(range(0, 105, 5)))
ymax = max(mapped_counts + paired_counts + [1])

svg = []
svg.append('<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760" viewBox="0 0 1200 760">')
svg.append('<rect width="100%" height="100%" fill="white"/>')
svg.append('<text x="600" y="38" text-anchor="middle" font-size="26" font-family="Arial" font-weight="bold">Symbiosis Mapping QC Across Samples V2</text>')
svg.append(f'<text x="600" y="68" text-anchor="middle" font-size="15" font-family="Arial">Samples summarized: {len(ok_rows)}; MC negative controls included: {overall["mc_control_samples"]}</text>')
left, _ = draw_hist_panel("A. Mapped reads", mapped_values, overall["mean_mapped_percent"], overall["median_mapped_percent"], "#64c7bb", "#087f3d", 105, 145, 450, 330, ymax)
right, _ = draw_hist_panel("B. Properly paired reads", paired_values, overall["mean_properly_paired_percent"], overall["median_properly_paired_percent"], "#aaa6d4", "#5b2ca0", 680, 145, 450, 330, ymax)
svg.extend(left)
svg.extend(right)
svg.append('<text x="600" y="690" text-anchor="middle" font-size="14" font-family="Arial">Mapping QC was calculated from samtools flagstat reports after mapping trimmed reads to symbiosis_islands.fasta.</text>')
svg.append('</svg>')
figure_path.write_text("\n".join(svg))

print(f"Wrote flagstat table: {flagstat_table}")
print(f"Wrote overall summary: {overall_path}")
print(f"Wrote sample-type summary: {by_type_path}")
print(f"Wrote sequencing-folder summary: {by_folder_path}")
print(f"Wrote SVG figure: {figure_path}")
PY

echo
echo "Overall mapping QC summary:"
column -t -s $'\t' "$OVERALL" 2>/dev/null || cat "$OVERALL"
echo
echo "Mapping QC outputs:"
echo "$COMPLETION"
echo "$FLAGSTAT_TABLE"
echo "$OVERALL"
echo "$BY_TYPE"
echo "$BY_FOLDER"
echo "$FIGURE"
