#!/usr/bin/env python3

import csv
import html
from collections import defaultdict
from pathlib import Path


BASE = Path("/Users/rosa/Documents/Ryan Data/symbiosis_sorted_v2")
COV = BASE / "nif_nod_coverage_existing_mapping_v2"
OUT = COV / "step10_figures"
OUT.mkdir(parents=True, exist_ok=True)

GENE_SUMMARY = COV / "nif_nod_gene_coverage_summary_by_sample_type_thresholds.tsv"
OVERALL = COV / "nif_nod_coverage_overall_threshold_summary.tsv"
TARGET_RANKING = COV / "nif_nod_target_mc_aware_ranking_thresholds.tsv"

HEATMAP = OUT / "v2_gene_coverage_heatmap_pct80_depth10.svg"
THRESHOLD_FIG = OUT / "v2_threshold_sensitivity_summary.svg"
TOP_TARGETS = OUT / "v2_top_mc_aware_targets_pct80_depth10.tsv"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def color_scale(value):
    # White to blue-green.
    value = max(0.0, min(1.0, value))
    r0, g0, b0 = 245, 248, 250
    r1, g1, b1 = 42, 145, 156
    r = round(r0 + (r1 - r0) * value)
    g = round(g0 + (g1 - g0) * value)
    b = round(b0 + (b1 - b0) * value)
    return f"rgb({r},{g},{b})"


def save_gene_heatmap(rows):
    threshold = "pct80_depth10"
    groups = ["MC", "No", "Rh", "Ro"]
    genes = sorted({r["gene"] for r in rows if r["threshold"] == threshold})
    values = {}

    for row in rows:
        if row["threshold"] != threshold:
            continue
        values[(row["gene"], row["sample_group"])] = float(row["good_sample_fraction"]) * 100

    cell_w = 98
    cell_h = 30
    left = 95
    top = 95
    width = left + cell_w * len(groups) + 230
    height = top + cell_h * len(genes) + 95

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append('<text x="420" y="32" text-anchor="middle" font-size="24" font-family="Arial" font-weight="bold">V2 nif/nod Gene Coverage by Sample Group</text>')
    svg.append('<text x="420" y="58" text-anchor="middle" font-size="14" font-family="Arial">Cell value = percent of samples with at least one good target; good = percent covered >= 80 and mean depth >= 10</text>')

    for j, group in enumerate(groups):
        x = left + j * cell_w + cell_w / 2
        svg.append(f'<text x="{x}" y="{top-18}" text-anchor="middle" font-size="15" font-family="Arial" font-weight="bold">{group}</text>')

    for i, gene in enumerate(genes):
        y = top + i * cell_h
        svg.append(f'<text x="{left-12}" y="{y+20}" text-anchor="end" font-size="13" font-family="Arial" font-weight="bold">{gene}</text>')
        for j, group in enumerate(groups):
            x = left + j * cell_w
            value = values.get((gene, group), 0.0)
            svg.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{color_scale(value/100)}" stroke="white"/>')
            text_color = "white" if value >= 55 else "#16202a"
            svg.append(f'<text x="{x+cell_w/2}" y="{y+20}" text-anchor="middle" font-size="12" font-family="Arial" fill="{text_color}">{value:.1f}%</text>')

    # Legend
    lx = left + cell_w * len(groups) + 48
    ly = top + 20
    svg.append(f'<text x="{lx}" y="{ly-12}" font-size="13" font-family="Arial" font-weight="bold">Good-sample fraction</text>')
    for k in range(0, 101, 10):
        y = ly + k * 2
        svg.append(f'<rect x="{lx}" y="{y}" width="24" height="20" fill="{color_scale(k/100)}" stroke="none"/>')
    svg.append(f'<text x="{lx+34}" y="{ly+8}" font-size="12" font-family="Arial">0%</text>')
    svg.append(f'<text x="{lx+34}" y="{ly+208}" font-size="12" font-family="Arial">100%</text>')

    svg.append('</svg>')
    HEATMAP.write_text("\n".join(svg))


def save_threshold_figure(rows):
    threshold_order = ["pct80_depth10", "pct70_depth10", "pct70_depth5", "pct50_depth5"]
    labels = {
        "pct80_depth10": "80/10",
        "pct70_depth10": "70/10",
        "pct70_depth5": "70/5",
        "pct50_depth5": "50/5",
    }
    groups = ["MC", "No", "Rh", "Ro"]
    data = {}
    for row in rows:
        samples = float(row["samples"])
        data[(row["threshold"], row["sample_group"], "any")] = float(row["samples_with_any_good_target"]) / samples * 100 if samples else 0
        data[(row["threshold"], row["sample_group"], "nif")] = float(row["samples_with_any_good_nif_target"]) / samples * 100 if samples else 0
        data[(row["threshold"], row["sample_group"], "nod")] = float(row["samples_with_any_good_nod_target"]) / samples * 100 if samples else 0
        data[(row["threshold"], row["sample_group"], "multi")] = float(row["samples_with_multi_target_signal_any_gene"]) / samples * 100 if samples else 0

    width = 1120
    height = 760
    left = 90
    top = 100
    panel_w = 450
    panel_h = 220

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append('<text x="560" y="34" text-anchor="middle" font-size="24" font-family="Arial" font-weight="bold">V2 Coverage Threshold Sensitivity</text>')
    svg.append('<text x="560" y="61" text-anchor="middle" font-size="14" font-family="Arial">Shows how many samples pass under different good-coverage thresholds</text>')

    def panel(title, metric, x, y, colors):
        svg.append(f'<text x="{x+panel_w/2}" y="{y-24}" text-anchor="middle" font-size="17" font-family="Arial" font-weight="bold">{html.escape(title)}</text>')
        svg.append(f'<rect x="{x}" y="{y}" width="{panel_w}" height="{panel_h}" fill="white" stroke="#1b2733"/>')
        for pct in [0, 50, 100]:
            yy = y + panel_h - pct / 100 * panel_h
            svg.append(f'<line x1="{x}" y1="{yy:.1f}" x2="{x+panel_w}" y2="{yy:.1f}" stroke="#e0e5ea"/>')
            svg.append(f'<text x="{x-9}" y="{yy+4:.1f}" text-anchor="end" font-size="11" font-family="Arial">{pct}%</text>')
        group_w = panel_w / len(threshold_order)
        bar_w = group_w / (len(groups) + 1)
        for i, threshold in enumerate(threshold_order):
            gx = x + i * group_w
            svg.append(f'<text x="{gx+group_w/2}" y="{y+panel_h+22}" text-anchor="middle" font-size="12" font-family="Arial">{labels[threshold]}</text>')
            for j, group in enumerate(groups):
                value = data.get((threshold, group, metric), 0)
                bh = value / 100 * panel_h
                bx = gx + 10 + j * bar_w
                by = y + panel_h - bh
                svg.append(f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bar_w-3:.1f}" height="{bh:.1f}" fill="{colors[group]}" opacity="0.82"/>')
        svg.append(f'<text x="{x+panel_w/2}" y="{y+panel_h+50}" text-anchor="middle" font-size="12" font-family="Arial">Threshold: percent covered / mean depth</text>')

    colors = {"MC": "#cf4f56", "No": "#2878b5", "Rh": "#5baa68", "Ro": "#8d6bb8"}
    panel("A. Any nif/nod target", "any", left, top, colors)
    panel("B. Any nif target", "nif", left + 565, top, colors)
    panel("C. Any nod target", "nod", left, top + 345, colors)
    panel("D. Multi-target signal", "multi", left + 565, top + 345, colors)

    lx = 445
    ly = 705
    for i, group in enumerate(groups):
        svg.append(f'<rect x="{lx+i*92}" y="{ly}" width="18" height="12" fill="{colors[group]}"/>')
        svg.append(f'<text x="{lx+24+i*92}" y="{ly+11}" font-size="12" font-family="Arial">{group}</text>')
    svg.append('</svg>')
    THRESHOLD_FIG.write_text("\n".join(svg))


def save_top_targets(rows):
    filtered = [
        r for r in rows
        if r["threshold"] == "pct80_depth10" and r["classification"] == "PROMISING_MC_AWARE"
    ]
    filtered.sort(key=lambda r: (
        -float(r["no_good_samples"]),
        float(r["mc_good_samples"]),
        -float(r["depth_ratio_no_vs_mc"]),
    ))
    fields = [
        "gene", "target_id", "no_good_samples", "mc_good_samples",
        "no_good_fraction", "mc_good_fraction",
        "no_median_mean_depth", "mc_median_mean_depth",
        "depth_ratio_no_vs_mc", "classification"
    ]
    with TOP_TARGETS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in filtered:
            writer.writerow({field: row[field] for field in fields})


def main():
    gene_rows = read_tsv(GENE_SUMMARY)
    overall_rows = read_tsv(OVERALL)
    target_rows = read_tsv(TARGET_RANKING)
    save_gene_heatmap(gene_rows)
    save_threshold_figure(overall_rows)
    save_top_targets(target_rows)
    print("Step 10 V2 figures/tables written:")
    print(HEATMAP)
    print(THRESHOLD_FIG)
    print(TOP_TARGETS)


if __name__ == "__main__":
    main()
