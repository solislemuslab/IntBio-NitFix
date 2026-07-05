#!/bin/bash

set -uo pipefail

OUTBASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
LOGS_DIR="$OUTBASE/fastp_logs"
QC_DIR="$OUTBASE/trimmed_fastp_QC_checking"

mkdir -p "$QC_DIR"

echo "Making V2 fastp per-base quality profile figures"
echo "Input fastp JSON folder: $LOGS_DIR"
echo "Output folder: $QC_DIR"
echo

if [[ ! -d "$LOGS_DIR" ]]; then
    echo "ERROR: fastp log folder does not exist: $LOGS_DIR"
    exit 1
fi

JSON_COUNT=$(find "$LOGS_DIR" -maxdepth 1 -type f -name '*.fastp.json' | wc -l)
echo "fastp JSON files found: $JSON_COUNT"

if [[ "$JSON_COUNT" -eq 0 ]]; then
    echo "ERROR: no fastp JSON files found."
    exit 1
fi

python3 - "$OUTBASE" <<'PY'
import json
import sys
from pathlib import Path

outbase = Path(sys.argv[1])
logs_dir = outbase / "fastp_logs"
qc_dir = outbase / "trimmed_fastp_QC_checking"

out_before = qc_dir / "fastp_quality_before_trimming_all_samples_lightpurple_mean_blue.svg"
out_after = qc_dir / "fastp_quality_after_trimming_all_samples_lightpurple_mean_blue.svg"

json_files = sorted(logs_dir.glob("*.fastp.json"))

r1_before_all = []
r2_before_all = []
r1_after_all = []
r2_after_all = []

def get_curve(data, section):
    return data.get(section, {}).get("quality_curves", {}).get("mean", [])

def add_curve(curves, curve):
    if curve:
        curves.append([float(x) for x in curve])

for path in json_files:
    data = json.loads(path.read_text())
    add_curve(r1_before_all, get_curve(data, "read1_before_filtering"))
    add_curve(r2_before_all, get_curve(data, "read2_before_filtering"))
    add_curve(r1_after_all, get_curve(data, "read1_after_filtering"))
    add_curve(r2_after_all, get_curve(data, "read2_after_filtering"))

def mean_curve(curves):
    if not curves:
        return []
    max_len = max(len(c) for c in curves)
    out = []
    for i in range(max_len):
        values = [c[i] for c in curves if i < len(c)]
        out.append(sum(values) / len(values))
    return out

def curve_points(curve, x, y, w, h, ymin=0, ymax=42):
    if not curve:
        return ""
    points = []
    for i, q in enumerate(curve):
        px = x + (i / max(len(curve) - 1, 1)) * w
        py = y + h - ((q - ymin) / (ymax - ymin)) * h
        points.append(f"{px:.1f},{py:.1f}")
    return " ".join(points)

def draw_axes(svg, title, x, y, w, h):
    svg.append(f'<text x="{x+w/2}" y="{y-28}" text-anchor="middle" font-size="18" font-family="Arial" font-weight="bold">{title}</text>')
    svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="white" stroke="black"/>')

    for q in [0, 10, 20, 30, 40]:
        py = y + h - (q / 42) * h
        color = "#cfcfcf"
        dash = ""
        if q in [20, 30]:
            color = "#cc0000"
            dash = ' stroke-dasharray="7,6"'
        svg.append(f'<line x1="{x}" y1="{py:.1f}" x2="{x+w}" y2="{py:.1f}" stroke="{color}"{dash}/>')
        svg.append(f'<text x="{x-8}" y="{py+4:.1f}" text-anchor="end" font-size="12" font-family="Arial">{q}</text>')

    for pos in [1, 50, 100, 150]:
        px = x + ((pos - 1) / 149) * w
        svg.append(f'<line x1="{px:.1f}" y1="{y+h}" x2="{px:.1f}" y2="{y+h+5}" stroke="black"/>')
        svg.append(f'<text x="{px:.1f}" y="{y+h+24}" text-anchor="middle" font-size="12" font-family="Arial">{pos}</text>')

    svg.append(f'<text x="{x+w/2}" y="{y+h+58}" text-anchor="middle" font-size="14" font-family="Arial">Base position</text>')
    svg.append(f'<text x="{x-58}" y="{y+h/2}" transform="rotate(-90 {x-58},{y+h/2})" text-anchor="middle" font-size="14" font-family="Arial">Mean quality score</text>')

def add_legend(svg, x, y):
    svg.append(f'<line x1="{x}" y1="{y}" x2="{x+42}" y2="{y}" stroke="#9E9AC8" stroke-width="2" opacity="0.45"/>')
    svg.append(f'<text x="{x+54}" y="{y+5}" font-size="13" font-family="Arial">individual samples</text>')
    y += 24
    svg.append(f'<line x1="{x}" y1="{y}" x2="{x+42}" y2="{y}" stroke="#0072B2" stroke-width="3"/>')
    svg.append(f'<text x="{x+54}" y="{y+5}" font-size="13" font-family="Arial">mean across all samples</text>')
    y += 24
    svg.append(f'<line x1="{x}" y1="{y}" x2="{x+42}" y2="{y}" stroke="#cc0000" stroke-dasharray="7,6"/>')
    svg.append(f'<text x="{x+54}" y="{y+5}" font-size="13" font-family="Arial">Q20 / Q30 thresholds</text>')

def write_figure(out_path, title, subtitle, r1_curves, r2_curves):
    r1_mean = mean_curve(r1_curves)
    r2_mean = mean_curve(r2_curves)

    svg = []
    svg.append('<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="620" viewBox="0 0 1200 620">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append(f'<text x="600" y="36" text-anchor="middle" font-size="26" font-family="Arial" font-weight="bold">{title}</text>')
    svg.append(f'<text x="600" y="64" text-anchor="middle" font-size="15" font-family="Arial">{subtitle}</text>')

    draw_axes(svg, "A. Read 1", 110, 135, 450, 280)
    draw_axes(svg, "B. Read 2", 690, 135, 450, 280)

    # Draw individual samples first, then the mean line on top.
    for curve in r1_curves:
        svg.append(f'<polyline points="{curve_points(curve, 110, 135, 450, 280)}" fill="none" stroke="#9E9AC8" stroke-width="0.45" opacity="0.28"/>')
    for curve in r2_curves:
        svg.append(f'<polyline points="{curve_points(curve, 690, 135, 450, 280)}" fill="none" stroke="#9E9AC8" stroke-width="0.45" opacity="0.28"/>')

    svg.append(f'<polyline points="{curve_points(r1_mean, 110, 135, 450, 280)}" fill="none" stroke="#0072B2" stroke-width="3.2"/>')
    svg.append(f'<polyline points="{curve_points(r2_mean, 690, 135, 450, 280)}" fill="none" stroke="#0072B2" stroke-width="3.2"/>')

    add_legend(svg, 450, 520)
    svg.append('</svg>')
    out_path.write_text("\n".join(svg))

write_figure(
    out_before,
    "Per-base Quality Profiles Before Trimming",
    f"Light purple lines = individual samples; blue line = mean across {len(json_files)} samples",
    r1_before_all,
    r2_before_all,
)

write_figure(
    out_after,
    "Per-base Quality Profiles After Trimming",
    f"Light purple lines = individual samples; blue line = mean across {len(json_files)} samples",
    r1_after_all,
    r2_after_all,
)

print(f"JSON files processed: {len(json_files)}")
print(f"Wrote: {out_before}")
print(f"Wrote: {out_after}")
PY
