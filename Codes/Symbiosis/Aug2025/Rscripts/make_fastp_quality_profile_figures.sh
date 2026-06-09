#!/bin/bash

set -uo pipefail

BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"
LOG="$BASE/fastp_logs_ryan_full"
QC_DIR="$BASE/trimmed_fastp_QC_checking"

mkdir -p "$QC_DIR"

python3 - <<'PY'
import json
from pathlib import Path

BASE = Path("/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted")
LOG = BASE / "fastp_logs_ryan_full"
QC_DIR = BASE / "trimmed_fastp_QC_checking"

MEAN_FIG = QC_DIR / "fastp_quality_mean_only.svg"
ALL_FIG = QC_DIR / "fastp_quality_all_samples_gray_mean_blue.svg"

json_files = sorted(LOG.glob("*.fastp.json"))

r1_before_all = []
r1_after_all = []
r2_before_all = []
r2_after_all = []

def get_curve(data, section):
    return data.get(section, {}).get("quality_curves", {}).get("mean", [])

for path in json_files:
    data = json.loads(path.read_text())

    r1b = get_curve(data, "read1_before_filtering")
    r1a = get_curve(data, "read1_after_filtering")
    r2b = get_curve(data, "read2_before_filtering")
    r2a = get_curve(data, "read2_after_filtering")

    if r1b:
        r1_before_all.append([float(x) for x in r1b])
    if r1a:
        r1_after_all.append([float(x) for x in r1a])
    if r2b:
        r2_before_all.append([float(x) for x in r2b])
    if r2a:
        r2_after_all.append([float(x) for x in r2a])

def mean_curve(curves):
    if not curves:
        return []
    max_len = max(len(c) for c in curves)
    out = []
    for i in range(max_len):
        values = [c[i] for c in curves if i < len(c)]
        out.append(sum(values) / len(values))
    return out

r1_before_mean = mean_curve(r1_before_all)
r1_after_mean = mean_curve(r1_after_all)
r2_before_mean = mean_curve(r2_before_all)
r2_after_mean = mean_curve(r2_after_all)

def curve_points(curve, x, y, w, h, ymin=0, ymax=42):
    if not curve:
        return ""
    pts = []
    n = len(curve)
    for i, q in enumerate(curve):
        px = x + (i / max(n - 1, 1)) * w
        py = y + h - ((q - ymin) / (ymax - ymin)) * h
        pts.append(f"{px:.1f},{py:.1f}")
    return " ".join(pts)

def draw_axes(svg, title, x, y, w, h):
    svg.append(f'<text x="{x+w/2}" y="{y-25}" text-anchor="middle" font-size="18" font-family="Arial" font-weight="bold">{title}</text>')
    svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="white" stroke="black"/>')

    for q in [0, 10, 20, 30, 40]:
        py = y + h - (q / 42) * h
        color = "#bdbdbd"
        dash = ""
        if q in [20, 30]:
            color = "#cc0000"
            dash = ' stroke-dasharray="6,4"'
        svg.append(f'<line x1="{x}" y1="{py:.1f}" x2="{x+w}" y2="{py:.1f}" stroke="{color}"{dash}/>')
        svg.append(f'<text x="{x-8}" y="{py+4:.1f}" text-anchor="end" font-size="12" font-family="Arial">{q}</text>')

    for pos in [1, 50, 100, 150]:
        px = x + ((pos - 1) / 149) * w
        svg.append(f'<line x1="{px:.1f}" y1="{y+h}" x2="{px:.1f}" y2="{y+h+5}" stroke="black"/>')
        svg.append(f'<text x="{px:.1f}" y="{y+h+22}" text-anchor="middle" font-size="12" font-family="Arial">{pos}</text>')

    svg.append(f'<text x="{x+w/2}" y="{y+h+50}" text-anchor="middle" font-size="14" font-family="Arial">Base position</text>')
    svg.append(f'<text x="{x-55}" y="{y+h/2}" transform="rotate(-90 {x-55},{y+h/2})" text-anchor="middle" font-size="14" font-family="Arial">Mean quality score</text>')

def add_legend(svg, x, y, include_gray=True):
    if include_gray:
        svg.append(f'<line x1="{x}" y1="{y}" x2="{x+40}" y2="{y}" stroke="#cccccc" stroke-width="2"/>')
        svg.append(f'<text x="{x+48}" y="{y+4}" font-size="13" font-family="Arial">individual samples</text>')
        y += 20
    svg.append(f'<line x1="{x}" y1="{y}" x2="{x+40}" y2="{y}" stroke="#777777" stroke-width="2"/>')
    svg.append(f'<text x="{x+48}" y="{y+4}" font-size="13" font-family="Arial">before trimming mean</text>')
    y += 20
    svg.append(f'<line x1="{x}" y1="{y}" x2="{x+40}" y2="{y}" stroke="#0072B2" stroke-width="3"/>')
    svg.append(f'<text x="{x+48}" y="{y+4}" font-size="13" font-family="Arial">after trimming mean</text>')
    y += 20
    svg.append(f'<line x1="{x}" y1="{y}" x2="{x+40}" y2="{y}" stroke="#cc0000" stroke-dasharray="6,4"/>')
    svg.append(f'<text x="{x+48}" y="{y+4}" font-size="13" font-family="Arial">Q20 / Q30 thresholds</text>')

def make_mean_figure(outfile):
    svg = []
    svg.append('<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="560" viewBox="0 0 1100 560">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append('<text x="550" y="35" text-anchor="middle" font-size="24" font-family="Arial" font-weight="bold">Mean Per-base Quality Before and After Trimming</text>')
    svg.append(f'<text x="550" y="62" text-anchor="middle" font-size="14" font-family="Arial">Samples summarized: {len(json_files)}</text>')

    draw_axes(svg, "A. Read 1", 90, 120, 420, 280)
    draw_axes(svg, "B. Read 2", 620, 120, 420, 280)

    svg.append(f'<polyline points="{curve_points(r1_before_mean, 90, 120, 420, 280)}" fill="none" stroke="#777777" stroke-width="2"/>')
    svg.append(f'<polyline points="{curve_points(r1_after_mean, 90, 120, 420, 280)}" fill="none" stroke="#0072B2" stroke-width="3"/>')
    svg.append(f'<polyline points="{curve_points(r2_before_mean, 620, 120, 420, 280)}" fill="none" stroke="#777777" stroke-width="2"/>')
    svg.append(f'<polyline points="{curve_points(r2_after_mean, 620, 120, 420, 280)}" fill="none" stroke="#0072B2" stroke-width="3"/>')

    add_legend(svg, 410, 470, include_gray=False)
    svg.append('</svg>')
    outfile.write_text("\n".join(svg))

def make_all_samples_figure(outfile):
    svg = []
    svg.append('<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="560" viewBox="0 0 1100 560">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append('<text x="550" y="35" text-anchor="middle" font-size="24" font-family="Arial" font-weight="bold">Per-base Quality Profiles Across All Samples</text>')
    svg.append(f'<text x="550" y="62" text-anchor="middle" font-size="14" font-family="Arial">Gray lines = individual samples; blue line = mean across {len(json_files)} samples</text>')

    draw_axes(svg, "A. Read 1 after trimming", 90, 120, 420, 280)
    draw_axes(svg, "B. Read 2 after trimming", 620, 120, 420, 280)

    for curve in r1_after_all:
        svg.append(f'<polyline points="{curve_points(curve, 90, 120, 420, 280)}" fill="none" stroke="#cccccc" stroke-width="0.4" opacity="0.25"/>')
    for curve in r2_after_all:
        svg.append(f'<polyline points="{curve_points(curve, 620, 120, 420, 280)}" fill="none" stroke="#cccccc" stroke-width="0.4" opacity="0.25"/>')

    svg.append(f'<polyline points="{curve_points(r1_after_mean, 90, 120, 420, 280)}" fill="none" stroke="#0072B2" stroke-width="3"/>')
    svg.append(f'<polyline points="{curve_points(r2_after_mean, 620, 120, 420, 280)}" fill="none" stroke="#0072B2" stroke-width="3"/>')

    add_legend(svg, 410, 470, include_gray=True)
    svg.append('</svg>')
    outfile.write_text("\n".join(svg))

make_mean_figure(MEAN_FIG)
make_all_samples_figure(ALL_FIG)

print(f"JSON files processed: {len(json_files)}")
print(f"Wrote: {MEAN_FIG}")
print(f"Wrote: {ALL_FIG}")
PY
