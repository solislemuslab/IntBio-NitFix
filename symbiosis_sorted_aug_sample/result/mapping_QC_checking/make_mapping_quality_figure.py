#!/usr/bin/env python3

import csv
from pathlib import Path


HERE = Path(__file__).resolve().parent
INFILE = HERE / "symbiosis_mapping_flagstat_summary.tsv"
OUTFILE = HERE / "symbiosis_mapping_quality_distribution.svg"


rows = []
with INFILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows.append(
            {
                "sample": row["sample"],
                "mapped_percent": float(row["mapped_percent"]),
                "properly_paired_percent": float(row["properly_paired_percent"]),
                "is_blank": row["sample"].startswith("BLAN"),
            }
        )


def mean(values):
    return sum(values) / len(values) if values else 0.0


def median(values):
    values = sorted(values)
    n = len(values)
    if n == 0:
        return 0.0
    mid = n // 2
    if n % 2:
        return values[mid]
    return (values[mid - 1] + values[mid]) / 2


def histogram(values, start=0, end=100, bin_width=5):
    bins = []
    x = start
    while x < end:
        bins.append((x, x + bin_width, 0))
        x += bin_width

    counts = [0 for _ in bins]
    for v in values:
        for i, (lo, hi, _) in enumerate(bins):
            if lo <= v < hi or (hi == end and v <= hi):
                counts[i] += 1
                break

    return [(lo, hi, count) for (lo, hi, _), count in zip(bins, counts)]


def draw_histogram(svg, title, values, x, y, w, h, color, mean_color):
    bins = histogram(values, 0, 100, 5)
    max_count = max(count for _, _, count in bins) or 1
    mean_value = mean(values)
    median_value = median(values)

    svg.append(
        f'<text x="{x + w / 2}" y="{y - 25}" text-anchor="middle" '
        f'font-size="18" font-family="Arial" font-weight="bold">{title}</text>'
    )
    svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="white" stroke="black"/>')

    for pct in [0, 25, 50, 75, 100]:
        px = x + (pct / 100) * w
        svg.append(f'<line x1="{px:.1f}" y1="{y}" x2="{px:.1f}" y2="{y + h}" stroke="#eeeeee"/>')
        svg.append(
            f'<text x="{px:.1f}" y="{y + h + 20}" text-anchor="middle" '
            f'font-size="12" font-family="Arial">{pct}</text>'
        )

    for tick in [0, max_count // 2, max_count]:
        py = y + h - (tick / max_count) * h
        svg.append(f'<line x1="{x - 5}" y1="{py:.1f}" x2="{x}" y2="{py:.1f}" stroke="black"/>')
        svg.append(
            f'<text x="{x - 10}" y="{py + 4:.1f}" text-anchor="end" '
            f'font-size="12" font-family="Arial">{tick}</text>'
        )

    bin_w = w / len(bins)
    for i, (lo, hi, count) in enumerate(bins):
        bh = (count / max_count) * h
        bx = x + i * bin_w
        by = y + h - bh
        svg.append(
            f'<rect x="{bx + 1:.1f}" y="{by:.1f}" width="{bin_w - 2:.1f}" '
            f'height="{bh:.1f}" fill="{color}" stroke="white"/>'
        )

    mx = x + (mean_value / 100) * w
    svg.append(
        f'<line x1="{mx:.1f}" y1="{y}" x2="{mx:.1f}" y2="{y + h}" '
        f'stroke="{mean_color}" stroke-width="3"/>'
    )
    svg.append(
        f'<text x="{mx:.1f}" y="{y - 7}" text-anchor="middle" '
        f'font-size="12" font-family="Arial" fill="{mean_color}">mean {mean_value:.1f}%</text>'
    )

    svg.append(
        f'<text x="{x + w / 2}" y="{y + h + 48}" text-anchor="middle" '
        f'font-size="14" font-family="Arial">Percent</text>'
    )
    svg.append(
        f'<text x="{x - 55}" y="{y + h / 2}" transform="rotate(-90 {x - 55},{y + h / 2})" '
        f'text-anchor="middle" font-size="14" font-family="Arial">Number of samples</text>'
    )
    svg.append(
        f'<text x="{x + w / 2}" y="{y + h + 70}" text-anchor="middle" '
        f'font-size="12" font-family="Arial">median {median_value:.1f}%</text>'
    )


mapped = [r["mapped_percent"] for r in rows]
proper = [r["properly_paired_percent"] for r in rows]
blank_count = sum(1 for r in rows if r["is_blank"])

width = 1100
height = 620
svg = []
svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
svg.append('<rect width="100%" height="100%" fill="white"/>')
svg.append(
    '<text x="550" y="35" text-anchor="middle" font-size="24" '
    'font-family="Arial" font-weight="bold">Symbiosis Mapping QC Across Samples</text>'
)
svg.append(
    f'<text x="550" y="63" text-anchor="middle" font-size="14" font-family="Arial">'
    f'Samples summarized: {len(rows)}; BLAN controls included: {blank_count}</text>'
)

draw_histogram(svg, "A. Mapped reads", mapped, 95, 130, 410, 310, "#8DD3C7", "#006D2C")
draw_histogram(svg, "B. Properly paired reads", proper, 635, 130, 410, 310, "#BEBADA", "#54278F")

svg.append(
    '<text x="550" y="585" text-anchor="middle" font-size="13" font-family="Arial">'
    'Mapping QC was calculated from samtools flagstat reports after mapping trimmed reads to symbiosis_islands.fasta.</text>'
)
svg.append("</svg>")

OUTFILE.write_text("\n".join(svg))

print(f"Samples summarized: {len(rows)}")
print(f"BLAN controls included: {blank_count}")
print(f"Mean mapped percent: {mean(mapped):.4f}")
print(f"Mean properly paired percent: {mean(proper):.4f}")
print(f"Wrote figure: {OUTFILE}")

