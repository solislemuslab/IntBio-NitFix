#!/usr/bin/env python3

import csv
import os
from pathlib import Path
from statistics import median
from collections import defaultdict

BASE = Path(os.environ.get(
    "BASE",
    "/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"
))
INFILE = BASE / "nif_nod_coverage_existing_mapping" / "nif_nod_region_coverage_all_samples.tsv"
OUTDIR = BASE / "nif_nod_coverage_existing_mapping"

OUTFILE = OUTDIR / "nif_nod_target_blank_aware_ranking_v2.tsv"

MIN_PERCENT_COVERED = 80.0
MIN_MEAN_DEPTH = 10.0
MIN_NO_GOOD_SAMPLES = 50
MAX_BLANK_GOOD_FRACTION = 0.25
MIN_DEPTH_RATIO_NO_VS_BLANK = 2.0
MIN_GOOD_FRACTION_DIFFERENCE = 0.10

def sample_group(row):
    if row["blank_control"] == "yes":
        return "BLAN"
    return row["sample_type"]

def is_good(row):
    return (
        float(row["percent_covered"]) >= MIN_PERCENT_COVERED
        and float(row["mean_depth"]) >= MIN_MEAN_DEPTH
    )

target_group_rows = defaultdict(list)

with open(INFILE, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    expected = [
        "sample", "sample_type", "blank_control", "gene",
        "original_reference", "bed_start", "bed_end",
        "target_id", "strand", "target_length",
        "covered_bases", "percent_covered", "mean_depth", "max_depth"
    ]

    if reader.fieldnames != expected:
        print("ERROR: Header does not match expected 14 columns.")
        print("Observed:")
        print(reader.fieldnames)
        print("Expected:")
        print(expected)
        raise SystemExit(1)

    for row in reader:
        key = (row["gene"], row["target_id"])
        target_group_rows[(key, sample_group(row))].append(row)

all_targets = sorted({key for key, grp in target_group_rows})

with open(OUTFILE, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow([
        "gene", "target_id",
        "no_samples", "blank_samples",
        "no_good_samples", "blank_good_samples",
        "no_good_fraction", "blank_good_fraction",
        "no_median_mean_depth", "blank_median_mean_depth",
        "no_median_percent_covered", "blank_median_percent_covered",
        "depth_ratio_no_vs_blank",
        "good_fraction_difference",
        "status"
    ])

    for gene, target_id in all_targets:
        no_rows = target_group_rows.get(((gene, target_id), "No"), [])
        blank_rows = target_group_rows.get(((gene, target_id), "BLAN"), [])

        no_samples = sorted({r["sample"] for r in no_rows})
        blank_samples = sorted({r["sample"] for r in blank_rows})

        no_good = sorted({r["sample"] for r in no_rows if is_good(r)})
        blank_good = sorted({r["sample"] for r in blank_rows if is_good(r)})

        no_depths = [float(r["mean_depth"]) for r in no_rows]
        blank_depths = [float(r["mean_depth"]) for r in blank_rows]

        no_pcts = [float(r["percent_covered"]) for r in no_rows]
        blank_pcts = [float(r["percent_covered"]) for r in blank_rows]

        no_good_frac = len(no_good) / len(no_samples) if no_samples else 0.0
        blank_good_frac = len(blank_good) / len(blank_samples) if blank_samples else 0.0

        no_med_depth = median(no_depths) if no_depths else 0.0
        blank_med_depth = median(blank_depths) if blank_depths else 0.0

        no_med_pct = median(no_pcts) if no_pcts else 0.0
        blank_med_pct = median(blank_pcts) if blank_pcts else 0.0

        depth_ratio = (no_med_depth + 0.01) / (blank_med_depth + 0.01)
        frac_diff = no_good_frac - blank_good_frac

        if (
            len(no_good) >= MIN_NO_GOOD_SAMPLES
            and blank_good_frac <= MAX_BLANK_GOOD_FRACTION
            and depth_ratio >= MIN_DEPTH_RATIO_NO_VS_BLANK
            and frac_diff >= MIN_GOOD_FRACTION_DIFFERENCE
        ):
            status = "PROMISING"
        elif len(no_good) >= MIN_NO_GOOD_SAMPLES and depth_ratio >= 1:
            status = "BIOLOGICAL_SIGNAL_BUT_BLANK_BACKGROUND"
        elif len(no_good) > 0:
            status = "LOW_OR_LIMITED_SIGNAL"
        else:
            status = "NOT_SUPPORTED"

        writer.writerow([
            gene, target_id,
            len(no_samples), len(blank_samples),
            len(no_good), len(blank_good),
            f"{no_good_frac:.4f}", f"{blank_good_frac:.4f}",
            f"{no_med_depth:.4f}", f"{blank_med_depth:.4f}",
            f"{no_med_pct:.4f}", f"{blank_med_pct:.4f}",
            f"{depth_ratio:.4f}",
            f"{frac_diff:.4f}",
            status
        ])

print("Blank-aware ranking v2 complete.")
print(OUTFILE)
