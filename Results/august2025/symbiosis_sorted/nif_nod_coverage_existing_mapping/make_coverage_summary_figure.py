#!/usr/bin/env python3

"""
Create a coverage-summary heatmap for central nif/nod genes.

Input:
  nif_nod_gene_coverage_summary_by_sample_type.tsv

Output:
  nif_nod_gene_coverage_by_sample_type_heatmap.svg

The figure shows, for each gene and sample group, the percentage of samples
with at least one target region passing the coverage threshold used in the
pipeline: percent_covered >= 80 and mean_depth >= 10.
"""

from pathlib import Path
import csv
import html


BASE = Path(__file__).resolve().parent
INPUT = BASE / "nif_nod_gene_coverage_summary_by_sample_type.tsv"
OUTPUT = BASE / "nif_nod_gene_coverage_by_sample_type_heatmap.svg"

GROUPS = ["BLAN", "No", "Rh", "Ro"]


def read_table(path):
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            unique_samples = float(row["unique_samples"])
            good_samples = float(row["samples_with_at_least_one_good_target"])
            row["good_sample_percent"] = 0.0 if unique_samples == 0 else 100.0 * good_samples / unique_samples
            rows.append(row)
    return rows


def color_for_percent(value):
    """White-to-blue color scale."""
    value = max(0.0, min(100.0, value))
    t = value / 100.0
    r0, g0, b0 = 245, 247, 251
    r1, g1, b1 = 31, 92, 160
    r = round(r0 + (r1 - r0) * t)
    g = round(g0 + (g1 - g0) * t)
    b = round(b0 + (b1 - b0) * t)
    return f"rgb({r},{g},{b})"


def text_color_for_percent(value):
    return "white" if value >= 55 else "#111111"


def main():
    if not INPUT.exists():
        raise SystemExit(f"Missing input file: {INPUT}")

    rows = read_table(INPUT)
    genes = sorted({row["gene"] for row in rows}, key=lambda x: (not x.startswith("nif"), x))
    data = {(row["gene"], row["sample_group"]): row for row in rows}

    cell_w = 92
    cell_h = 26
    left = 95
    top = 125
    right = 175
    bottom = 105
    width = left + len(GROUPS) * cell_w + right
    height = top + len(genes) * cell_h + bottom

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append('<style>')
    svg.append('text { font-family: Arial, Helvetica, sans-serif; fill: #111; }')
    svg.append('.title { font-size: 24px; font-weight: 700; }')
    svg.append('.subtitle { font-size: 13px; fill: #333; }')
    svg.append('.axis { font-size: 13px; font-weight: 700; }')
    svg.append('.gene { font-size: 12px; font-weight: 700; }')
    svg.append('.value { font-size: 11px; font-weight: 700; }')
    svg.append('.note { font-size: 11px; fill: #444; }')
    svg.append('</style>')

    svg.append(f'<text x="{width/2}" y="38" text-anchor="middle" class="title">Central nif/nod Coverage by Sample Type</text>')
    svg.append(f'<text x="{width/2}" y="61" text-anchor="middle" class="subtitle">Color = percent of samples with at least one good target region per gene</text>')
    svg.append(f'<text x="{width/2}" y="81" text-anchor="middle" class="subtitle">Good coverage threshold: percent covered >= 80 and mean depth >= 10</text>')

    for i, group in enumerate(GROUPS):
        x = left + i * cell_w + cell_w / 2
        svg.append(f'<text x="{x}" y="{top - 18}" text-anchor="middle" class="axis">{html.escape(group)}</text>')

    for j, gene in enumerate(genes):
        y = top + j * cell_h
        svg.append(f'<text x="{left - 12}" y="{y + 18}" text-anchor="end" class="gene">{html.escape(gene)}</text>')
        for i, group in enumerate(GROUPS):
            x = left + i * cell_w
            row = data.get((gene, group))
            value = row["good_sample_percent"] if row else 0.0
            fill = color_for_percent(value)
            txt = text_color_for_percent(value)
            svg.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{fill}" stroke="white" stroke-width="1"/>')
            svg.append(f'<text x="{x + cell_w/2}" y="{y + 17}" text-anchor="middle" class="value" fill="{txt}">{value:.1f}%</text>')

    # Legend
    leg_x = left
    leg_y = top + len(genes) * cell_h + 35
    svg.append(f'<text x="{leg_x}" y="{leg_y - 10}" class="axis">Percent of samples with good coverage</text>')
    legend_w = 250
    legend_h = 14
    for k in range(101):
        x = leg_x + k * legend_w / 100
        svg.append(f'<rect x="{x:.2f}" y="{leg_y}" width="{legend_w/100 + 1:.2f}" height="{legend_h}" fill="{color_for_percent(k)}" stroke="none"/>')
    for tick in [0, 25, 50, 75, 100]:
        x = leg_x + tick * legend_w / 100
        svg.append(f'<line x1="{x}" y1="{leg_y + legend_h}" x2="{x}" y2="{leg_y + legend_h + 5}" stroke="#333"/>')
        svg.append(f'<text x="{x}" y="{leg_y + legend_h + 20}" text-anchor="middle" class="note">{tick}%</text>')

    note = (
        "BLAN controls are shown for QC comparison; biological interpretation should focus on No, Rh, and Ro samples."
    )
    svg.append(f'<text x="{width/2}" y="{height - 25}" text-anchor="middle" class="note">{html.escape(note)}</text>')
    svg.append("</svg>")

    OUTPUT.write_text("\n".join(svg) + "\n")
    print(f"Wrote figure: {OUTPUT}")
    print(f"Genes plotted: {len(genes)}")
    print(f"Sample groups plotted: {', '.join(GROUPS)}")


if __name__ == "__main__":
    main()
