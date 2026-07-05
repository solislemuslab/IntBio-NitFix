#!/usr/bin/env python3

import csv
import os
from collections import defaultdict
from pathlib import Path
from statistics import mean, median


BASE = Path(os.environ.get(
    "BASE",
    "/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
))

OUTDIR = BASE / "nif_nod_coverage_existing_mapping_v2"
INFILE = OUTDIR / "nif_nod_region_coverage_all_samples.tsv"

THRESHOLD_TEXT = os.environ.get("THRESHOLDS", "80:10,70:10,70:5,50:5")

GENE_GROUP_SUMMARY = OUTDIR / "nif_nod_gene_coverage_summary_by_sample_type_thresholds.tsv"
TARGET_MC_RANKING = OUTDIR / "nif_nod_target_mc_aware_ranking_thresholds.tsv"
SAMPLE_GENE_COUNTS = OUTDIR / "nif_nod_sample_gene_good_target_counts_thresholds.tsv"
OVERALL = OUTDIR / "nif_nod_coverage_overall_threshold_summary.tsv"


def parse_thresholds(text):
    thresholds = []
    for item in text.split(","):
        pct, depth = item.split(":")
        thresholds.append((float(pct), float(depth), f"pct{pct}_depth{depth}"))
    return thresholds


def group_label(row):
    if row["sample_type"] == "MC" or row["is_mc_control"] == "yes":
        return "MC"
    return row["sample_type"]


def is_good(row, min_percent, min_depth):
    return (
        float(row["percent_covered"]) >= min_percent
        and float(row["mean_depth"]) >= min_depth
    )


def safe_mean(values):
    values = list(values)
    return mean(values) if values else 0.0


def safe_median(values):
    values = list(values)
    return median(values) if values else 0.0


def fmt(value):
    return f"{value:.4f}"


def classify_target(no_fraction, mc_fraction, no_samples, mc_samples, no_depth, mc_depth):
    if no_samples == 0:
        return "NO_NODULE_SAMPLES"
    if no_fraction >= 0.10 and mc_fraction <= 0.05 and no_depth > mc_depth:
        return "PROMISING_MC_AWARE"
    if no_fraction >= 0.10 and mc_fraction > 0.05:
        return "BIOLOGICAL_SIGNAL_BUT_MC_BACKGROUND"
    if no_fraction > 0:
        return "LOW_OR_LIMITED_SIGNAL"
    return "NOT_SUPPORTED"


if not INFILE.exists():
    raise SystemExit(f"ERROR: coverage table not found: {INFILE}")

thresholds = parse_thresholds(THRESHOLD_TEXT)

rows = []
with INFILE.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        row["_group"] = group_label(row)
        row["_percent_covered"] = float(row["percent_covered"])
        row["_mean_depth"] = float(row["mean_depth"])
        rows.append(row)

sample_groups = defaultdict(set)
for row in rows:
    sample_groups[row["_group"]].add(row["sample"])

with GENE_GROUP_SUMMARY.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "threshold",
        "min_percent_covered",
        "min_mean_depth",
        "gene",
        "sample_group",
        "rows",
        "unique_samples",
        "targets_tested",
        "good_rows",
        "samples_with_at_least_one_good_target",
        "good_sample_fraction",
        "median_percent_covered",
        "median_mean_depth",
        "mean_percent_covered",
        "mean_mean_depth",
        "max_mean_depth",
    ])

    for min_percent, min_depth, threshold_name in thresholds:
        gene_group = defaultdict(list)
        for row in rows:
            gene_group[(row["gene"], row["_group"])].append(row)

        for (gene, group), group_rows in sorted(gene_group.items()):
            samples = sorted({r["sample"] for r in group_rows})
            targets = sorted({r["target_id"] for r in group_rows})
            good_rows = [r for r in group_rows if is_good(r, min_percent, min_depth)]
            good_samples = sorted({r["sample"] for r in good_rows})
            pct_values = [r["_percent_covered"] for r in group_rows]
            depth_values = [r["_mean_depth"] for r in group_rows]
            good_fraction = len(good_samples) / len(samples) if samples else 0

            writer.writerow([
                threshold_name,
                min_percent,
                min_depth,
                gene,
                group,
                len(group_rows),
                len(samples),
                len(targets),
                len(good_rows),
                len(good_samples),
                fmt(good_fraction),
                fmt(safe_median(pct_values)),
                fmt(safe_median(depth_values)),
                fmt(safe_mean(pct_values)),
                fmt(safe_mean(depth_values)),
                fmt(max(depth_values) if depth_values else 0),
            ])

with TARGET_MC_RANKING.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "threshold",
        "min_percent_covered",
        "min_mean_depth",
        "gene",
        "target_id",
        "no_samples",
        "mc_samples",
        "no_good_samples",
        "mc_good_samples",
        "no_good_fraction",
        "mc_good_fraction",
        "no_median_mean_depth",
        "mc_median_mean_depth",
        "no_median_percent_covered",
        "mc_median_percent_covered",
        "depth_ratio_no_vs_mc",
        "good_fraction_difference",
        "classification",
    ])

    target_rows = defaultdict(list)
    for row in rows:
        target_rows[(row["gene"], row["target_id"])].append(row)

    for min_percent, min_depth, threshold_name in thresholds:
        for (gene, target_id), target_group_rows in sorted(target_rows.items()):
            no_rows = [r for r in target_group_rows if r["_group"] == "No"]
            mc_rows = [r for r in target_group_rows if r["_group"] == "MC"]

            no_samples = sorted({r["sample"] for r in no_rows})
            mc_samples = sorted({r["sample"] for r in mc_rows})
            no_good = sorted({r["sample"] for r in no_rows if is_good(r, min_percent, min_depth)})
            mc_good = sorted({r["sample"] for r in mc_rows if is_good(r, min_percent, min_depth)})

            no_fraction = len(no_good) / len(no_samples) if no_samples else 0
            mc_fraction = len(mc_good) / len(mc_samples) if mc_samples else 0
            no_depth = safe_median([r["_mean_depth"] for r in no_rows])
            mc_depth = safe_median([r["_mean_depth"] for r in mc_rows])
            depth_ratio = no_depth / mc_depth if mc_depth > 0 else (999999.0 if no_depth > 0 else 0.0)
            classification = classify_target(
                no_fraction,
                mc_fraction,
                len(no_samples),
                len(mc_samples),
                no_depth,
                mc_depth,
            )

            writer.writerow([
                threshold_name,
                min_percent,
                min_depth,
                gene,
                target_id,
                len(no_samples),
                len(mc_samples),
                len(no_good),
                len(mc_good),
                fmt(no_fraction),
                fmt(mc_fraction),
                fmt(no_depth),
                fmt(mc_depth),
                fmt(safe_median([r["_percent_covered"] for r in no_rows])),
                fmt(safe_median([r["_percent_covered"] for r in mc_rows])),
                fmt(depth_ratio),
                fmt(no_fraction - mc_fraction),
                classification,
            ])

with SAMPLE_GENE_COUNTS.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "threshold",
        "min_percent_covered",
        "min_mean_depth",
        "sample",
        "sequencing_folder",
        "site",
        "sample_type",
        "sample_group",
        "gene",
        "good_target_regions_for_gene",
        "good_target_ids",
        "possible_multi_target_signal",
    ])

    sample_gene_rows = defaultdict(list)
    for row in rows:
        sample_gene_rows[(row["sample"], row["gene"])].append(row)

    for min_percent, min_depth, threshold_name in thresholds:
        for (sample, gene), group_rows in sorted(sample_gene_rows.items()):
            good_rows = [r for r in group_rows if is_good(r, min_percent, min_depth)]
            first = group_rows[0]
            good_target_ids = sorted({r["target_id"] for r in good_rows})
            writer.writerow([
                threshold_name,
                min_percent,
                min_depth,
                sample,
                first["sequencing_folder"],
                first["site"],
                first["sample_type"],
                first["_group"],
                gene,
                len(good_target_ids),
                ";".join(good_target_ids),
                "yes" if len(good_target_ids) > 1 else "no",
            ])

with OVERALL.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "threshold",
        "min_percent_covered",
        "min_mean_depth",
        "sample_group",
        "samples",
        "samples_with_any_good_target",
        "samples_with_any_good_nif_target",
        "samples_with_any_good_nod_target",
        "samples_with_multi_target_signal_any_gene",
    ])

    for min_percent, min_depth, threshold_name in thresholds:
        group_samples = defaultdict(set)
        any_good = defaultdict(set)
        any_good_nif = defaultdict(set)
        any_good_nod = defaultdict(set)
        sample_gene_good_targets = defaultdict(set)

        for row in rows:
            group = row["_group"]
            sample = row["sample"]
            gene = row["gene"]
            group_samples[group].add(sample)

            if is_good(row, min_percent, min_depth):
                any_good[group].add(sample)
                if gene.startswith("nif"):
                    any_good_nif[group].add(sample)
                if gene.startswith("nod"):
                    any_good_nod[group].add(sample)
                sample_gene_good_targets[(group, sample, gene)].add(row["target_id"])

        multi_by_group = defaultdict(set)
        for (group, sample, gene), targets in sample_gene_good_targets.items():
            if len(targets) > 1:
                multi_by_group[group].add(sample)

        for group in sorted(group_samples):
            writer.writerow([
                threshold_name,
                min_percent,
                min_depth,
                group,
                len(group_samples[group]),
                len(any_good[group]),
                len(any_good_nif[group]),
                len(any_good_nod[group]),
                len(multi_by_group[group]),
            ])

print("V2 nif/nod coverage summary complete.")
print(f"Input rows: {len(rows)}")
print(f"Thresholds: {', '.join(t[2] for t in thresholds)}")
print()
print(GENE_GROUP_SUMMARY)
print(TARGET_MC_RANKING)
print(SAMPLE_GENE_COUNTS)
print(OVERALL)
