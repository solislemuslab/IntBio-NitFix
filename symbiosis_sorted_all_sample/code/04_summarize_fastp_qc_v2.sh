#!/bin/bash

set -uo pipefail

OUTBASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
MANIFEST="$OUTBASE/metadata/symbiosis_v2_clean_manifest.tsv"
LOGS_DIR="$OUTBASE/fastp_logs"
TRIMMED_DIR="$OUTBASE/symbiosis_trimmed_fastp"
QC_DIR="$OUTBASE/trimmed_fastp_QC_checking"

SUMMARY="$QC_DIR/fastp_logs_summary.tsv"
OVERALL="$QC_DIR/fastp_qc_overall_summary.tsv"
BY_FOLDER="$QC_DIR/fastp_qc_by_sequencing_folder.tsv"
BY_TYPE="$QC_DIR/fastp_qc_by_sample_type.tsv"
FIGURE="$QC_DIR/fastp_qc_summary.svg"

mkdir -p "$QC_DIR"

echo "Checking V2 fastp trimming quality"
echo "V2 output folder: $OUTBASE"
echo "Input fastp JSON folder: $LOGS_DIR"
echo "QC output folder: $QC_DIR"
echo

if [[ ! -s "$MANIFEST" ]]; then
    echo "ERROR: manifest not found: $MANIFEST"
    exit 1
fi

if [[ ! -d "$LOGS_DIR" ]]; then
    echo "ERROR: fastp log folder does not exist: $LOGS_DIR"
    exit 1
fi

JSON_COUNT=$(find "$LOGS_DIR" -maxdepth 1 -type f -name '*.fastp.json' | wc -l)
P1_COUNT=$(find "$TRIMMED_DIR" -maxdepth 1 -type f -name '*_P1.fastq.gz' 2>/dev/null | wc -l)
P2_COUNT=$(find "$TRIMMED_DIR" -maxdepth 1 -type f -name '*_P2.fastq.gz' 2>/dev/null | wc -l)

echo "fastp JSON files found: $JSON_COUNT"
echo "Trimmed P1 files found: $P1_COUNT"
echo "Trimmed P2 files found: $P2_COUNT"
echo

if [[ "$JSON_COUNT" -eq 0 ]]; then
    echo "ERROR: no fastp JSON files found. Wait until trimming has produced JSON logs."
    exit 1
fi

python3 - "$OUTBASE" <<'PY'
import csv
import html
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

outbase = Path(sys.argv[1])
manifest = outbase / "metadata" / "symbiosis_v2_clean_manifest.tsv"
logs_dir = outbase / "fastp_logs"
trimmed_dir = outbase / "symbiosis_trimmed_fastp"
qc_dir = outbase / "trimmed_fastp_QC_checking"

summary_path = qc_dir / "fastp_logs_summary.tsv"
overall_path = qc_dir / "fastp_qc_overall_summary.tsv"
by_folder_path = qc_dir / "fastp_qc_by_sequencing_folder.tsv"
by_type_path = qc_dir / "fastp_qc_by_sample_type.tsv"
figure_path = qc_dir / "fastp_qc_summary.svg"

def mean(values):
    values = [v for v in values if v is not None]
    return statistics.mean(values) if values else None

def pct(value):
    return value * 100 if value is not None else None

def fmt(value, digits=2):
    if value is None:
        return "NA"
    return f"{value:.{digits}f}"

def read_manifest(path):
    samples = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            samples[row["sample"]] = row
    return samples

def get_curve(data, section):
    return data.get(section, {}).get("quality_curves", {}).get("mean", [])

def add_curve(collection, curve):
    if not curve:
        return
    while len(collection) < len(curve):
        collection.append([])
    for i, value in enumerate(curve):
        collection[i].append(float(value))

def curve_mean(collection):
    return [mean(values) for values in collection]

samples = read_manifest(manifest)
rows = []

r1_before = []
r1_after = []
r2_before = []
r2_after = []

for sample, meta in sorted(samples.items()):
    json_path = logs_dir / f"{sample}.fastp.json"
    p1_path = trimmed_dir / f"{sample}_P1.fastq.gz"
    p2_path = trimmed_dir / f"{sample}_P2.fastq.gz"

    row = {
        "sample": sample,
        "sequencing_folder": meta.get("sequencing_folder", ""),
        "site": meta.get("site", ""),
        "sample_type": meta.get("sample_type", ""),
        "json_found": "yes" if json_path.exists() else "no",
        "trimmed_P1_found": "yes" if p1_path.exists() else "no",
        "trimmed_P2_found": "yes" if p2_path.exists() else "no",
        "before_total_reads": None,
        "after_total_reads": None,
        "reads_retained_percent": None,
        "before_total_bases": None,
        "after_total_bases": None,
        "bases_retained_percent": None,
        "before_q20_percent": None,
        "after_q20_percent": None,
        "before_q30_percent": None,
        "after_q30_percent": None,
        "before_gc_percent": None,
        "after_gc_percent": None,
        "status": "missing_fastp_json",
    }

    if json_path.exists():
        try:
            data = json.loads(json_path.read_text())
            before = data["summary"]["before_filtering"]
            after = data["summary"]["after_filtering"]

            before_reads = before.get("total_reads", 0)
            after_reads = after.get("total_reads", 0)
            before_bases = before.get("total_bases", 0)
            after_bases = after.get("total_bases", 0)

            row.update({
                "before_total_reads": before_reads,
                "after_total_reads": after_reads,
                "reads_retained_percent": (after_reads / before_reads * 100) if before_reads else None,
                "before_total_bases": before_bases,
                "after_total_bases": after_bases,
                "bases_retained_percent": (after_bases / before_bases * 100) if before_bases else None,
                "before_q20_percent": pct(before.get("q20_rate")),
                "after_q20_percent": pct(after.get("q20_rate")),
                "before_q30_percent": pct(before.get("q30_rate")),
                "after_q30_percent": pct(after.get("q30_rate")),
                "before_gc_percent": pct(before.get("gc_content")),
                "after_gc_percent": pct(after.get("gc_content")),
                "status": "ok",
            })

            add_curve(r1_before, get_curve(data, "read1_before_filtering"))
            add_curve(r1_after, get_curve(data, "read1_after_filtering"))
            add_curve(r2_before, get_curve(data, "read2_before_filtering"))
            add_curve(r2_after, get_curve(data, "read2_after_filtering"))
        except Exception as exc:
            row["status"] = f"json_parse_error:{type(exc).__name__}"

    rows.append(row)

cols = [
    "sample",
    "sequencing_folder",
    "site",
    "sample_type",
    "json_found",
    "trimmed_P1_found",
    "trimmed_P2_found",
    "before_total_reads",
    "after_total_reads",
    "reads_retained_percent",
    "before_total_bases",
    "after_total_bases",
    "bases_retained_percent",
    "before_q20_percent",
    "after_q20_percent",
    "before_q30_percent",
    "after_q30_percent",
    "before_gc_percent",
    "after_gc_percent",
    "status",
]

with summary_path.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(cols)
    for row in rows:
        writer.writerow([
            row["sample"],
            row["sequencing_folder"],
            row["site"],
            row["sample_type"],
            row["json_found"],
            row["trimmed_P1_found"],
            row["trimmed_P2_found"],
            row["before_total_reads"] if row["before_total_reads"] is not None else "NA",
            row["after_total_reads"] if row["after_total_reads"] is not None else "NA",
            fmt(row["reads_retained_percent"]),
            row["before_total_bases"] if row["before_total_bases"] is not None else "NA",
            row["after_total_bases"] if row["after_total_bases"] is not None else "NA",
            fmt(row["bases_retained_percent"]),
            fmt(row["before_q20_percent"]),
            fmt(row["after_q20_percent"]),
            fmt(row["before_q30_percent"]),
            fmt(row["after_q30_percent"]),
            fmt(row["before_gc_percent"]),
            fmt(row["after_gc_percent"]),
            row["status"],
        ])

ok_rows = [r for r in rows if r["status"] == "ok"]
missing_rows = [r for r in rows if r["json_found"] == "no"]
failed_rows = [r for r in rows if r["status"] not in ("ok", "missing_fastp_json")]

overall_metrics = {
    "manifest_samples": len(rows),
    "fastp_json_ok": len(ok_rows),
    "fastp_json_missing": len(missing_rows),
    "fastp_json_parse_failed": len(failed_rows),
    "mean_reads_retained_percent": mean([r["reads_retained_percent"] for r in ok_rows]),
    "mean_bases_retained_percent": mean([r["bases_retained_percent"] for r in ok_rows]),
    "mean_q20_before_percent": mean([r["before_q20_percent"] for r in ok_rows]),
    "mean_q20_after_percent": mean([r["after_q20_percent"] for r in ok_rows]),
    "mean_q30_before_percent": mean([r["before_q30_percent"] for r in ok_rows]),
    "mean_q30_after_percent": mean([r["after_q30_percent"] for r in ok_rows]),
    "mean_gc_before_percent": mean([r["before_gc_percent"] for r in ok_rows]),
    "mean_gc_after_percent": mean([r["after_gc_percent"] for r in ok_rows]),
}

with overall_path.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["metric", "value"])
    for key, value in overall_metrics.items():
        if isinstance(value, float):
            value = fmt(value)
        writer.writerow([key, value])

def grouped_summary(group_key, out_path):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row[group_key]].append(row)

    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            group_key,
            "manifest_samples",
            "fastp_json_ok",
            "fastp_json_missing",
            "mean_reads_retained_percent",
            "mean_q30_before_percent",
            "mean_q30_after_percent",
        ])
        for group in sorted(grouped):
            group_rows = grouped[group]
            group_ok = [r for r in group_rows if r["status"] == "ok"]
            writer.writerow([
                group,
                len(group_rows),
                len(group_ok),
                len([r for r in group_rows if r["json_found"] == "no"]),
                fmt(mean([r["reads_retained_percent"] for r in group_ok])),
                fmt(mean([r["before_q30_percent"] for r in group_ok])),
                fmt(mean([r["after_q30_percent"] for r in group_ok])),
            ])

grouped_summary("sequencing_folder", by_folder_path)
grouped_summary("sample_type", by_type_path)

def points_for_curve(curve, x, y, width, height, ymin=0, ymax=42):
    if not curve:
        return ""
    points = []
    for i, value in enumerate(curve):
        px = x + (i / max(len(curve) - 1, 1)) * width
        py = y + height - ((value - ymin) / (ymax - ymin)) * height
        points.append(f"{px:.1f},{py:.1f}")
    return " ".join(points)

def draw_quality_panel(title, before_curve, after_curve, x, y, width, height):
    parts = []
    parts.append(f'<text x="{x + width / 2}" y="{y - 18}" text-anchor="middle" font-size="17" font-family="Arial" font-weight="bold">{html.escape(title)}</text>')
    parts.append(f'<rect x="{x}" y="{y}" width="{width}" height="{height}" fill="white" stroke="#172033"/>')
    for q in [0, 10, 20, 30, 40]:
        py = y + height - (q / 42) * height
        color = "#d8dee8"
        dash = ""
        if q in (20, 30):
            color = "#c75146"
            dash = ' stroke-dasharray="5,5"'
        parts.append(f'<line x1="{x}" y1="{py:.1f}" x2="{x + width}" y2="{py:.1f}" stroke="{color}"{dash}/>')
        parts.append(f'<text x="{x - 8}" y="{py + 4:.1f}" text-anchor="end" font-size="11" font-family="Arial">{q}</text>')
    for pos in [1, 50, 100, 150]:
        px = x + ((pos - 1) / 149) * width
        parts.append(f'<line x1="{px:.1f}" y1="{y + height}" x2="{px:.1f}" y2="{y + height + 5}" stroke="#172033"/>')
        parts.append(f'<text x="{px:.1f}" y="{y + height + 20}" text-anchor="middle" font-size="11" font-family="Arial">{pos}</text>')
    parts.append(f'<polyline points="{points_for_curve(before_curve, x, y, width, height)}" fill="none" stroke="#8f8f9d" stroke-width="2.4"/>')
    parts.append(f'<polyline points="{points_for_curve(after_curve, x, y, width, height)}" fill="none" stroke="#2474b5" stroke-width="2.4"/>')
    parts.append(f'<text x="{x + width / 2}" y="{y + height + 44}" text-anchor="middle" font-size="13" font-family="Arial">Base position</text>')
    parts.append(f'<text x="{x - 48}" y="{y + height / 2}" transform="rotate(-90 {x - 48},{y + height / 2})" text-anchor="middle" font-size="13" font-family="Arial">Mean quality</text>')
    return parts

def bar(x, y, width, height, value, color, label):
    capped = min(max(value or 0, 0), 100)
    bar_height = height * capped / 100
    parts = [
        f'<rect x="{x}" y="{y}" width="{width}" height="{height}" fill="#f6f7fb" stroke="#d5dce8"/>',
        f'<rect x="{x}" y="{y + height - bar_height:.1f}" width="{width}" height="{bar_height:.1f}" fill="{color}"/>',
        f'<text x="{x + width / 2}" y="{y - 8}" text-anchor="middle" font-size="13" font-family="Arial" font-weight="bold">{fmt(value)}</text>',
        f'<text x="{x + width / 2}" y="{y + height + 18}" text-anchor="middle" font-size="11" font-family="Arial">{html.escape(label)}</text>',
    ]
    return parts

r1_before_mean = curve_mean(r1_before)
r1_after_mean = curve_mean(r1_after)
r2_before_mean = curve_mean(r2_before)
r2_after_mean = curve_mean(r2_after)

svg = []
svg.append('<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="760" viewBox="0 0 1120 760">')
svg.append('<rect width="100%" height="100%" fill="white"/>')
svg.append('<text x="560" y="36" text-anchor="middle" font-size="25" font-family="Arial" font-weight="bold">symbiosis_sorted_v2 fastp QC</text>')
svg.append(f'<text x="560" y="64" text-anchor="middle" font-size="14" font-family="Arial">Manifest samples: {len(rows)}; completed JSON: {len(ok_rows)}; missing JSON: {len(missing_rows)}</text>')
svg.extend(draw_quality_panel("A. Read 1 quality before vs after trimming", r1_before_mean, r1_after_mean, 90, 125, 420, 240))
svg.extend(draw_quality_panel("B. Read 2 quality before vs after trimming", r2_before_mean, r2_after_mean, 625, 125, 420, 240))

svg.append('<line x1="410" y1="430" x2="450" y2="430" stroke="#8f8f9d" stroke-width="3"/>')
svg.append('<text x="460" y="434" font-size="13" font-family="Arial">Before trimming</text>')
svg.append('<line x1="580" y1="430" x2="620" y2="430" stroke="#2474b5" stroke-width="3"/>')
svg.append('<text x="630" y="434" font-size="13" font-family="Arial">After trimming</text>')
svg.append('<line x1="755" y1="430" x2="795" y2="430" stroke="#c75146" stroke-dasharray="5,5"/>')
svg.append('<text x="805" y="434" font-size="13" font-family="Arial">Q20/Q30 guide lines</text>')

svg.append('<text x="560" y="500" text-anchor="middle" font-size="17" font-family="Arial" font-weight="bold">C. Mean sample-level metrics</text>')
metrics = [
    ("reads retained", overall_metrics["mean_reads_retained_percent"], "#76b77a"),
    ("bases retained", overall_metrics["mean_bases_retained_percent"], "#76b77a"),
    ("Q30 before", overall_metrics["mean_q30_before_percent"], "#9e9eaa"),
    ("Q30 after", overall_metrics["mean_q30_after_percent"], "#2474b5"),
    ("Q20 after", overall_metrics["mean_q20_after_percent"], "#2f9aa0"),
]
start_x = 250
for i, (label, value, color) in enumerate(metrics):
    svg.extend(bar(start_x + i * 130, 545, 78, 120, value, color, label))
svg.append('<text x="195" y="548" text-anchor="end" font-size="11" font-family="Arial">100%</text>')
svg.append('<text x="195" y="665" text-anchor="end" font-size="11" font-family="Arial">0%</text>')
svg.append('<line x1="205" y1="545" x2="205" y2="665" stroke="#172033"/>')
svg.append('<line x1="205" y1="665" x2="930" y2="665" stroke="#172033"/>')
svg.append('</svg>')
figure_path.write_text("\n".join(svg))

print(f"Wrote per-sample summary: {summary_path}")
print(f"Wrote overall summary:    {overall_path}")
print(f"Wrote folder summary:     {by_folder_path}")
print(f"Wrote sample-type summary:{by_type_path}")
print(f"Wrote SVG figure:         {figure_path}")
PY

echo
echo "Overall QC summary:"
column -t -s $'\t' "$OVERALL" 2>/dev/null || cat "$OVERALL"
echo
echo "QC outputs:"
echo "$SUMMARY"
echo "$OVERALL"
echo "$BY_FOLDER"
echo "$BY_TYPE"
echo "$FIGURE"
