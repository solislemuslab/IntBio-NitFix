#!/usr/bin/env python3

import argparse
import csv
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def read_manifest(path):
    metadata = {}
    if not path or not Path(path).exists():
        return metadata

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            metadata[row["sample"]] = row
    return metadata


def fallback_metadata(sample):
    if sample.startswith("MC-"):
        return {"sequencing_folder": "unknown", "site": "MC", "sample_type": "MC"}

    pieces = sample.split("-")
    sample_type = pieces[-1] if pieces else "unknown"
    site = pieces[0] if pieces else "unknown"
    return {"sequencing_folder": "unknown", "site": site, "sample_type": sample_type}


def read_regions(path):
    regions = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for i, row in enumerate(reader):
            start = int(row["bed_start"])
            end = int(row["bed_end"])
            length = int(row["target_length"])
            regions.append({
                "region_index": i,
                "original_reference": row["original_reference"],
                "bed_start": start,
                "bed_end": end,
                "gene": row["gene"],
                "target_id": row["extracted_target_id"],
                "strand": row["strand"],
                "target_length": length,
                "sum_depth": 0,
                "covered_bases": 0,
                "max_depth": 0,
            })
    return regions


def write_bed(regions, path):
    with open(path, "w") as handle:
        for region in regions:
            name = f"{region['gene']}|{region['target_id']}|region{region['region_index']}"
            handle.write(
                "\t".join([
                    region["original_reference"],
                    str(region["bed_start"]),
                    str(region["bed_end"]),
                    name,
                    "0",
                    region["strand"],
                ]) + "\n"
            )


def calculate_depth(bam, regions):
    by_ref = defaultdict(list)
    for region in regions:
        by_ref[region["original_reference"]].append(region)

    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tmp:
        bed_path = Path(tmp.name)
    try:
        write_bed(regions, bed_path)
        cmd = ["samtools", "depth", "-a", "-b", str(bed_path), str(bam)]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        assert proc.stdout is not None
        for line in proc.stdout:
            ref, pos_text, depth_text = line.rstrip("\n").split("\t")[:3]
            pos_1based = int(pos_text)
            pos_0based = pos_1based - 1
            depth = int(depth_text)

            for region in by_ref.get(ref, []):
                if region["bed_start"] <= pos_0based < region["bed_end"]:
                    region["sum_depth"] += depth
                    if depth > 0:
                        region["covered_bases"] += 1
                    if depth > region["max_depth"]:
                        region["max_depth"] = depth

        return_code = proc.wait()
        if return_code != 0:
            raise SystemExit(f"ERROR: samtools depth failed for {bam}")
    finally:
        bed_path.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--regions", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    bam = Path(args.bam)
    sample = bam.name.removesuffix(".bam")
    metadata = read_manifest(args.manifest).get(sample, fallback_metadata(sample))
    sequencing_folder = metadata.get("sequencing_folder", "unknown")
    site = metadata.get("site", "unknown")
    sample_type = metadata.get("sample_type", "unknown")
    is_mc_control = "yes" if sample_type == "MC" or sample.startswith("MC-") else "no"

    regions = read_regions(args.regions)
    calculate_depth(bam, regions)

    columns = [
        "sample",
        "sequencing_folder",
        "site",
        "sample_type",
        "is_mc_control",
        "gene",
        "original_reference",
        "bed_start",
        "bed_end",
        "target_id",
        "strand",
        "target_length",
        "covered_bases",
        "percent_covered",
        "mean_depth",
        "max_depth",
    ]

    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for region in regions:
            length = region["target_length"]
            covered = region["covered_bases"]
            mean_depth = region["sum_depth"] / length if length else 0
            percent_covered = (covered / length * 100) if length else 0
            writer.writerow([
                sample,
                sequencing_folder,
                site,
                sample_type,
                is_mc_control,
                region["gene"],
                region["original_reference"],
                region["bed_start"],
                region["bed_end"],
                region["target_id"],
                region["strand"],
                length,
                covered,
                f"{percent_covered:.4f}",
                f"{mean_depth:.4f}",
                region["max_depth"],
            ])


if __name__ == "__main__":
    main()
