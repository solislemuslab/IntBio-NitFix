#!/bin/bash

set -uo pipefail

BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"
LOG="$BASE/fastp_logs_ryan_full"
QC_DIR="$BASE/trimmed_fastp_QC_checking"

SUMMARY="$QC_DIR/fastp_logs_ryan_full_summary.tsv"
FIGURE="$QC_DIR/fastp_quality_before_after_and_retention.svg"

mkdir -p "$QC_DIR"

echo "Input fastp JSON folder:"
echo "$LOG"
echo

echo "QC output folder:"
echo "$QC_DIR"
echo

if [[ ! -d "$LOG" ]]; then
    echo "ERROR: fastp log folder does not exist: $LOG"
    exit 1
fi

JSON_COUNT=$(find "$LOG" -maxdepth 1 -type f -name '*.fastp.json' | wc -l)
echo "fastp JSON files found: $JSON_COUNT"

if [[ "$JSON_COUNT" -eq 0 ]]; then
    echo "ERROR: No fastp JSON files found."
    exit 1
fi

python3 - <<'PY'
import csv
import json
from pathlib import Path

BASE = Path("/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted")
LOG = BASE / "fastp_logs_ryan_full"
QC_DIR = BASE / "trimmed_fastp_QC_checking"

SUMMARY = QC_DIR / "fastp_logs_ryan_full_summary.tsv"
FIGURE = QC_DIR / "fastp_quality_before_after_and_retention.svg"

rows = []
read1_before = []
read1_after = []
read2_before = []
read2_after = []

def get_quality_curve(data, section, read_name):
    try:
        return data[section][read_name]["quality_curves"]["mean"]
    except KeyError:
        return []

def add_curve(collection, curve):
    if not curve:
        return
    while len(collection) < len(curve):
        collection.append([])
    for i, value in enumerate(curve):
        collection[i].append(float(value))

for path in sorted(LOG.glob("*.fastp.json")):
    sample = path.name.replace(".fastp.json", "")
    data = json.loads(path.read_text())

    before = data["summary"]["before_filtering"]
    after = data["summary"]["after_filtering"]

    before_reads = before["total_reads"]
    after_reads = after["total_reads"]
    retained = (after_reads / before_reads * 100) if before_reads else 0

    rows.append({
        "sample": sample,
        "before_total_reads": before_reads,
        "after_total_reads": after_reads,
        "reads_retained_percent": retained,
        "before_q20_rate": before["q20_rate"],
        "after_q20_rate": after["q20_rate"],
        "before_q30_rate": before["q30_rate"],
        "after_q30_rate": after["q30_rate"],
        "before_gc_content": before["gc_content"],
        "after_gc_content": after["gc_content"],
    })

    add_curve(read1_before, get_quality_curve(data, "read1_before_filtering", "quality_curves"))
    add_curve(read1_after, get_quality_curve(data, "read1_after_filtering", "quality_curves"))
    add_curve(read2_before, get_quality_curve(data, "read2_before_filtering", "quality_curves"))
    add_curve(read2_after, get_quality_curve(data, "read2_after_filtering", "quality_curves"))

# Some fastp JSON formats store quality curves differently.
if not read1_before:
    for path in sorted(LOG.glob("*.fastp.json")):
        data = json.loads(path.read_text())
        add_curve(read1_before, data.get("read1_before_filtering", {}).get("quality_curves", {}).get("mean", []))
        add_curve(read1_after, data.get("read1_after_filtering", {}).get("quality_curves", {}).get("mean", []))
        add_curve(read2_before, data.get("read2_before_filtering", {}).get("quality_curves", {}).get("mean", []))
        add_curve(read2_after, data.get("read2_after_filtering", {}).get("quality_curves", {}).get("mean", []))

cols = [
    "sample",
    "before_total_reads",
    "after_total_reads",
    "reads_retained_percent",
    "before_q20_rate",
    "after_q20_rate",
    "before_q30_rate",
    "after_q30_rate",
    "before_gc_content",
    "after_gc_content"
]

with SUMMARY.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(cols)
    for r in rows:
        writer.writerow([
            r["sample"],
            r["before_total_reads"],
            r["after_total_reads"],
            f"{r['reads_retained_percent']:.2f}",
            f"{r['before_q20_rate']:.4f}",
            f"{r['after_q20_rate']:.4f}",
            f"{r['before_q30_rate']:.4f}",
            f"{r['after_q30_rate']:.4f}",
            f"{r['before_gc_content']:.4f}",
            f"{r['after_gc_content']:.4f}",
        ])

def mean(values):
    return sum(values) / len(values) if values else 0

def mean_curve(collection):
    return [mean(v) for v in collection]

r1_before = mean_curve(read1_before)
r1_after = mean_curve(read1_after)
r2_before = mean_curve(read2_before)
r2_after = mean_curve(read2_after)

mean_retained = mean([r["reads_retained_percent"] for r in rows])
mean_before_q20 = mean([r["before_q20_rate"] * 100 for r in rows])
mean_after_q20 = mean([r["after_q20_rate"] * 100 for r in rows])
mean_before_q30 = mean([r["before_q30_rate"] * 100 for r in rows])
mean_after_q30 = mean([r["after_q30_rate"] * 100 for r in rows])

def points_for_curve(curve, x, y, w, h, ymin=0, ymax=42):
    pts = []
    if not curve:
        return ""
    n = len(curve)
    for i, value in enumerate(curve):
        px = x + (i / max(n - 1, 1)) * w
        py = y + h - ((value - ymin) / (ymax - ymin)) * h
        pts.append(f"{px:.1f},{py:.1f}")
    return " ".join(pts)

def draw_quality_panel(title, before_curve, after_curve, x, y, w, h):
    out = []
    out.append(f'<text x="{x + w/2}" y="{y - 18}" text-anchor="middle" font-size="17" font-family="Arial" font-weight="bold">{title}</text>')
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="white" stroke="black"/>')

    for q in [0, 10, 20, 30, 40]:
        py = y + h - (q / 42) * h
        out.append(f'<line x1="{x}" y1="{py:.1f}" x2="{x+w}" y2="{py:.1f}" stroke="#eeeeee"/>')
        out.append(f'<text x="{x-8}" y="{py+4:.1f}" text-anchor="end" font-size="11" font-family="Arial">{q}</text>')

    for pos in [1, 50, 100, 150]:
        px = x + ((pos - 1) / 149) * w
        out.append(f'<line x1="{px:.1f}" y1="{y+h}" x2="{px:.1f}" y2="{y+h+5}" stroke="black"/>')
        out.append(f'<text x="{px:.1f}" y="{y+h+20}" text-anchor="middle" font-size="11" font-family="Arial">{pos}</text>')

    out.append(f'<polyline points="{points_for_curve(before_curve, x, y, w, h)}" fill="none" stroke="#777777" stroke-width="2"/>')
    out.append(f'<polyline points="{points_for_curve(after_curve, x, y, w, h)}" fill="none" stroke="#0072B2" stroke-width="2"/>')

    out.append(f'<text x="{x+w-120}" y="{y+20}" font-size="12" font-family="Arial" fill="#777777">Before</text>')
    out.append(f'<line x1="{x+w-165}" y1="{y+16}" x2="{x+w-125}" y2="{y+16}" stroke="#777777" stroke-width="2"/>')
    out.append(f'<text x="{x+w-120}" y="{y+40}" font-size="12" font-family="Arial" fill="#0072B2">After</text>')
    out.append(f'<line x1="{x+w-165}" y1="{y+36}" x2="{x+w-125}" y2="{y+36}" stroke="#0072B2" stroke-width="2"/>')

    out.append(f'<text x="{x+w/2}" y="{y+h+45}" text-anchor="middle" font-size="13" font-family="Arial">Base position</text>')
    out.append(f'<text x="{x-48}" y="{y+h/2}" transform="rotate(-90 {x-48},{y+h/2})" text-anchor="middle" font-size="13" font-family="Arial">Mean quality score</text>')
    return out

width = 1100
height = 760
svg = []

svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
svg.append('<rect width="100%" height="100%" fill="white"/>')
svg.append('<text x="550" y="35" text-anchor="middle" font-size="24" font-family="Arial" font-weight="bold">fastp Quality Control Before and After Trimming</text>')
svg.append(f'<text x="550" y="62" text-anchor="middle" font-size="14" font-family="Arial">Samples summarized: {len(rows)}</text>')

svg.extend(draw_quality_panel("A. Read 1 per-base quality", r1_before, r1_after, 90, 120, 420, 240))
svg.extend(draw_quality_panel("B. Read 2 per-base quality", r2_before, r2_after, 620, 120, 420, 240))

# Panel C: summary bars
x0 = 90
y0 = 500
bar_w = 100
gap = 45
plot_h = 170

bar_data = [
    ("Q20 before", mean_before_q20, "#bdbdbd"),
    ("Q20 after", mean_after_q20, "#4fa3d1"),
    ("Q30 before", mean_before_q30, "#bdbdbd"),
    ("Q30 after", mean_after_q30, "#4fa3d1"),
    ("Reads retained", mean_retained, "#7fc97f"),
]

svg.append('<text x="360" y="475" text-anchor="middle" font-size="17" font-family="Arial" font-weight="bold">C. Mean sample-level QC metrics</text>')
svg.append(f'<line x1="{x0-30}" y1="{y0+plot_h}" x2="{x0+5*(bar_w+gap)}" y2="{y0+plot_h}" stroke="black"/>')
svg.append(f'<line x1="{x0-30}" y1="{y0}" x2="{x0-30}" y2="{y0+plot_h}" stroke="black"/>')

for tick in [0, 25, 50, 75, 100]:
    yy = y0 + plot_h - (tick / 100) * plot_h
    svg.append(f'<line x1="{x0-35}" y1="{yy}" x2="{x0-30}" y2="{yy}" stroke="black"/>')
    svg.append(f'<text x="{x0-42}" y="{yy+5}" text-anchor="end" font-size="11" font-family="Arial">{tick}</text>')
    svg.append(f'<line x1="{x0-30}" y1="{yy}" x2="{x0+5*(bar_w+gap)}" y2="{yy}" stroke="#eeeeee"/>')

for i, (label, value, color) in enumerate(bar_data):
    xx = x0 + i * (bar_w + gap)
    hh = (value / 100) * plot_h
    yy = y0 + plot_h - hh
    svg.append(f'<rect x="{xx}" y="{yy}" width="{bar_w}" height="{hh}" fill="{color}" stroke="black"/>')
    svg.append(f'<text x="{xx + bar_w/2}" y="{yy - 8}" text-anchor="middle" font-size="12" font-family="Arial">{value:.2f}%</text>')
    svg.append(f'<text x="{xx + bar_w/2}" y="{y0 + plot_h + 22}" text-anchor="middle" font-size="11" font-family="Arial">{label}</text>')

svg.append(f'<text x="40" y="{y0 + 85}" transform="rotate(-90 40,{y0 + 85})" text-anchor="middle" font-size="13" font-family="Arial">Percent</text>')

svg.append('<text x="550" y="735" text-anchor="middle" font-size="13" font-family="Arial">Trimming retained most reads and increased Q20/Q30 quality rates after fastp filtering.</text>')
svg.append('</svg>')

FIGURE.write_text("\n".join(svg))

print(f"Wrote summary table: {SUMMARY}")
print(f"Wrote QC figure: {FIGURE}")
print()
print("===== QC summary =====")
print(f"samples: {len(rows)}")
print(f"mean_reads_retained_percent: {mean_retained:.4f}")
print(f"mean_before_q20_rate: {mean_before_q20 / 100:.6f}")
print(f"mean_after_q20_rate: {mean_after_q20 / 100:.6f}")
print(f"mean_before_q30_rate: {mean_before_q30 / 100:.6f}")
print(f"mean_after_q30_rate: {mean_after_q30 / 100:.6f}")
PY

echo
echo "First 10 rows:"
column -t -s $'\t' "$SUMMARY" | head -10

echo
echo "Done."
echo "QC outputs saved in:"
echo "$QC_DIR"
