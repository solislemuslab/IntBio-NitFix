#!/usr/bin/env python3

import csv
from pathlib import Path
from statistics import median, mean
from collections import defaultdict

BASE = Path("/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted")

INFILE = BASE / "nif_nod_coverage_existing_mapping" / "nif_nod_region_coverage_all_samples.tsv"
OUTDIR = BASE / "nif_nod_coverage_existing_mapping"

GENE_GROUP_SUMMARY = OUTDIR / "nif_nod_gene_coverage_summary_by_sample_type.tsv"
TARGET_SUMMARY = OUTDIR / "nif_nod_target_coverage_summary.tsv"
SAMPLE_SUMMARY = OUTDIR / "nif_nod_sample_coverage_summary.tsv"

MIN_PERCENT_COVERED = 80.0
MIN_MEAN_DEPTH = 10.0

def sample_group(row):
    if row["blank_control"] == "yes":
        return "BLAN"
    return row["sample_type"]

def good(row):
    return (
        float(row["percent_covered"]) >= MIN_PERCENT_COVERED
        and float(row["mean_depth"]) >= MIN_MEAN_DEPTH
    )

gene_group = defaultdict(list)
target_group = defaultdict(list)
sample_rows = defaultdict(list)

with open(INFILE, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        group = sample_group(row)
        gene = row["gene"]
        target = row["target_id"]
        sample = row["sample"]

        row["_group"] = group
        row["_good"] = good(row)

        gene_group[(gene, group)].append(row)
        target_group[(gene, target, group)].append(row)
        sample_rows[(sample, group)].append(row)

with open(GENE_GROUP_SUMMARY, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "gene", "sample_group", "rows", "unique_samples",
        "targets_tested", "good_rows",
        "samples_with_at_least_one_good_target",
        "median_percent_covered", "median_mean_depth",
        "mean_percent_covered", "mean_mean_depth",
        "max_mean_depth"
    ])

    for (gene, group), rows in sorted(gene_group.items()):
        samples = sorted({r["sample"] for r in rows})
        targets = sorted({r["target_id"] for r in rows})
        good_rows = [r for r in rows if r["_good"]]
        good_samples = sorted({r["sample"] for r in good_rows})

        pct = [float(r["percent_covered"]) for r in rows]
        dep = [float(r["mean_depth"]) for r in rows]

        writer.writerow([
            gene, group, len(rows), len(samples), len(targets),
            len(good_rows), len(good_samples),
            f"{median(pct):.4f}", f"{median(dep):.4f}",
            f"{mean(pct):.4f}", f"{mean(dep):.4f}",
            f"{max(dep):.4f}"
        ])

with open(TARGET_SUMMARY, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "gene", "target_id", "sample_group", "rows",
        "unique_samples", "good_samples",
        "median_percent_covered", "median_mean_depth",
        "max_mean_depth"
    ])

    for (gene, target, group), rows in sorted(target_group.items()):
        samples = sorted({r["sample"] for r in rows})
        good_samples = sorted({r["sample"] for r in rows if r["_good"]})

        pct = [float(r["percent_covered"]) for r in rows]
        dep = [float(r["mean_depth"]) for r in rows]

        writer.writerow([
            gene, target, group, len(rows), len(samples), len(good_samples),
            f"{median(pct):.4f}", f"{median(dep):.4f}",
            f"{max(dep):.4f}"
        ])

with open(SAMPLE_SUMMARY, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "sample", "sample_group", "good_target_regions",
        "good_genes", "median_percent_covered",
        "median_mean_depth", "max_mean_depth"
    ])

    for (sample, group), rows in sorted(sample_rows.items()):
        good_rows = [r for r in rows if r["_good"]]
        good_genes = sorted({r["gene"] for r in good_rows})

        pct = [float(r["percent_covered"]) for r in rows]
        dep = [float(r["mean_depth"]) for r in rows]

        writer.writerow([
            sample, group, len(good_rows), len(good_genes),
            f"{median(pct):.4f}", f"{median(dep):.4f}",
            f"{max(dep):.4f}"
        ])

print("Coverage summary complete.")
print(f"Good coverage threshold: percent_covered >= {MIN_PERCENT_COVERED}, mean_depth >= {MIN_MEAN_DEPTH}")
print()
print(GENE_GROUP_SUMMARY)
print(TARGET_SUMMARY)
print(SAMPLE_SUMMARY)
