#!/usr/bin/env python3

import csv
import os
import re
import subprocess
from collections import Counter
from pathlib import Path


BASE = Path(os.environ.get(
    "BASE",
    "/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
))

TARGET_ID = os.environ.get(
    "TARGET_ID",
    "nifH|NC_009937|NC_009937_-_Symbiosis_Island_4|ref56"
)

SAMPLE_GROUP = os.environ.get("SAMPLE_GROUP", "No")
MIN_PERCENT_COVERED = float(os.environ.get("MIN_PERCENT_COVERED", "80"))
MIN_MEAN_DEPTH = float(os.environ.get("MIN_MEAN_DEPTH", "10"))
MIN_BASE_DEPTH = int(os.environ.get("MIN_BASE_DEPTH", "5"))
MAX_N_PERCENT = float(os.environ.get("MAX_N_PERCENT", "5"))
MAX_SAMPLES = int(os.environ.get("MAX_SAMPLES", "0"))

MIXED_MIN_DEPTH = int(os.environ.get("MIXED_MIN_DEPTH", "20"))
MIXED_MINOR_COUNT = int(os.environ.get("MIXED_MINOR_COUNT", "5"))
MIXED_MINOR_AF = float(os.environ.get("MIXED_MINOR_AF", "0.20"))
MAX_MIXED_POSITIONS_FOR_SINGLE = int(os.environ.get("MAX_MIXED_POSITIONS_FOR_SINGLE", "10"))

COVERAGE = BASE / "nif_nod_coverage_existing_mapping_v2" / "nif_nod_region_coverage_all_samples.tsv"
REGIONS = BASE / "nif_nod_original_reference_regions" / "nif_nod_matches_in_original_reference.tsv"
BAM_DIR = BASE / "symbiosis_mapped_full"
REF = BASE / "reference" / "symbiosis_islands.fasta"


def slugify(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")


OUTDIR = BASE / "multicopy_aware_target_sequences_v2" / slugify(TARGET_ID)
SEQDIR = OUTDIR / "sequences"
LOGDIR = OUTDIR / "logs"
OUTDIR.mkdir(parents=True, exist_ok=True)
SEQDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

SAMPLE_LIST = OUTDIR / "selected_samples.tsv"
DOMINANT_FASTA = OUTDIR / "dominant_consensus.fasta"
IUPAC_FASTA = OUTDIR / "iupac_ambiguity_consensus.fasta"
QC_TABLE = OUTDIR / "sample_consensus_qc.tsv"
MIXED_TABLE = OUTDIR / "mixed_sites.tsv"
README = OUTDIR / "README.txt"


IUPAC = {
    frozenset(["A"]): "A",
    frozenset(["C"]): "C",
    frozenset(["G"]): "G",
    frozenset(["T"]): "T",
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
    frozenset(["A", "C", "G"]): "V",
    frozenset(["A", "C", "T"]): "H",
    frozenset(["A", "G", "T"]): "D",
    frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "C", "G", "T"]): "N",
}

RC = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def reverse_complement(seq):
    return seq.translate(RC)[::-1].upper()


def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def parse_regions():
    matches = []
    with REGIONS.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["extracted_target_id"] == TARGET_ID:
                matches.append(row)
    if not matches:
        raise SystemExit(f"ERROR: target not found in region table: {TARGET_ID}")
    if len(matches) > 1:
        print(f"WARNING: target has {len(matches)} exact locations. Using the first one.")
    row = matches[0]
    row["bed_start"] = int(row["bed_start"])
    row["bed_end"] = int(row["bed_end"])
    row["target_length"] = int(row["target_length"])
    return row


def select_samples():
    selected = []
    with COVERAGE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["target_id"] != TARGET_ID:
                continue
            if SAMPLE_GROUP != "all" and row["sample_type"] != SAMPLE_GROUP:
                continue
            if row["sample_type"] == "MC" or row.get("is_mc_control", "no") == "yes":
                continue
            if float(row["percent_covered"]) < MIN_PERCENT_COVERED:
                continue
            if float(row["mean_depth"]) < MIN_MEAN_DEPTH:
                continue
            selected.append(row)
    selected.sort(key=lambda r: r["sample"])
    return selected


def parse_pileup_bases(base_string, ref_base):
    counts = Counter()
    i = 0
    ref_base = ref_base.upper()
    while i < len(base_string):
        char = base_string[i]
        if char == "^":
            i += 2
            continue
        if char == "$":
            i += 1
            continue
        if char in "+-":
            i += 1
            number = []
            while i < len(base_string) and base_string[i].isdigit():
                number.append(base_string[i])
                i += 1
            indel_len = int("".join(number)) if number else 0
            i += indel_len
            continue
        if char in ".,":
            if ref_base in "ACGT":
                counts[ref_base] += 1
            i += 1
            continue
        base = char.upper()
        if base in "ACGT":
            counts[base] += 1
        i += 1
    return counts


def pileup_sample(sample, region):
    bam = BAM_DIR / f"{sample}.bam"
    if not bam.exists():
        raise FileNotFoundError(f"Missing BAM for {sample}: {bam}")

    samtools_region = f"{region['original_reference']}:{region['bed_start'] + 1}-{region['bed_end']}"
    cmd = [
        "samtools", "mpileup",
        "-aa",
        "-A",
        "-Q", "0",
        "-q", "0",
        "-f", str(REF),
        "-r", samtools_region,
        str(bam),
    ]

    proc = subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    (LOGDIR / f"{sample}.mpileup.stderr.log").write_text(proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"samtools mpileup failed for {sample}")

    by_offset = {}
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        ref_name, pos_text, ref_base, depth_text, bases = parts[:5]
        ref_pos = int(pos_text)
        offset = ref_pos - (region["bed_start"] + 1)
        counts = parse_pileup_bases(bases, ref_base)
        by_offset[offset] = {
            "ref_pos": ref_pos,
            "ref_base": ref_base.upper(),
            "counts": counts,
        }
    return by_offset


def call_consensus(sample, region):
    pileup = pileup_sample(sample, region)
    dominant = []
    iupac = []
    mixed_rows = []
    depths = []
    covered_bases = 0

    for offset in range(region["target_length"]):
        site = pileup.get(offset, {"ref_pos": region["bed_start"] + offset + 1, "ref_base": "N", "counts": Counter()})
        counts = site["counts"]
        depth = sum(counts.values())
        depths.append(depth)

        if depth < MIN_BASE_DEPTH:
            dominant.append("N")
            iupac.append("N")
            continue

        covered_bases += 1
        ordered = counts.most_common()
        major_base, major_count = ordered[0]
        minor_base, minor_count = ordered[1] if len(ordered) > 1 else ("", 0)
        minor_af = minor_count / depth if depth else 0
        dominant.append(major_base)

        mixed_bases = {major_base}
        for base, count in ordered[1:]:
            af = count / depth if depth else 0
            if count >= MIXED_MINOR_COUNT and af >= MIXED_MINOR_AF:
                mixed_bases.add(base)

        iupac.append(IUPAC.get(frozenset(mixed_bases), "N"))

        if (
            depth >= MIXED_MIN_DEPTH
            and minor_count >= MIXED_MINOR_COUNT
            and minor_af >= MIXED_MINOR_AF
        ):
            target_pos = offset + 1 if region["strand"] == "+" else region["target_length"] - offset
            mixed_rows.append({
                "sample": sample,
                "target_id": TARGET_ID,
                "reference": region["original_reference"],
                "ref_pos_1based": site["ref_pos"],
                "target_pos_1based": target_pos,
                "depth": depth,
                "major_base": major_base,
                "major_count": major_count,
                "minor_base": minor_base,
                "minor_count": minor_count,
                "minor_allele_fraction": minor_af,
                "A": counts.get("A", 0),
                "C": counts.get("C", 0),
                "G": counts.get("G", 0),
                "T": counts.get("T", 0),
            })

    dominant_seq = "".join(dominant)
    iupac_seq = "".join(iupac)
    if region["strand"] == "-":
        dominant_seq = reverse_complement(dominant_seq)
        iupac_seq = reverse_complement(iupac_seq)

    n_count = dominant_seq.count("N")
    n_percent = n_count / len(dominant_seq) * 100 if dominant_seq else 100
    mixed_count = len(mixed_rows)
    max_minor_af = max((r["minor_allele_fraction"] for r in mixed_rows), default=0.0)
    mean_depth = sum(depths) / len(depths) if depths else 0
    median_depth = sorted(depths)[len(depths) // 2] if depths else 0
    percent_covered = covered_bases / region["target_length"] * 100 if region["target_length"] else 0

    if n_percent > MAX_N_PERCENT:
        status = "fail_missing_data"
    elif mixed_count > MAX_MIXED_POSITIONS_FOR_SINGLE:
        status = "pass_mixed_possible_multicopy"
    else:
        status = "pass_single_dominant"

    qc = {
        "sample": sample,
        "target_id": TARGET_ID,
        "length": len(dominant_seq),
        "covered_bases": covered_bases,
        "percent_covered_from_pileup": percent_covered,
        "mean_depth_from_pileup": mean_depth,
        "median_depth_from_pileup": median_depth,
        "N_count": n_count,
        "N_percent": n_percent,
        "mixed_positions": mixed_count,
        "max_minor_allele_fraction": max_minor_af,
        "status": status,
    }
    return dominant_seq, iupac_seq, qc, mixed_rows


def main():
    if not COVERAGE.exists():
        raise SystemExit(f"ERROR: coverage table not found: {COVERAGE}")
    if not REGIONS.exists():
        raise SystemExit(f"ERROR: region table not found: {REGIONS}")
    if not REF.exists():
        raise SystemExit(f"ERROR: reference FASTA not found: {REF}")

    region = parse_regions()
    samples = select_samples()
    if MAX_SAMPLES > 0:
        samples = samples[:MAX_SAMPLES]

    print("Multi-copy-aware target extraction V2")
    print(f"Target: {TARGET_ID}")
    print(f"Region: {region['original_reference']}:{region['bed_start'] + 1}-{region['bed_end']} ({region['strand']})")
    print(f"Sample group: {SAMPLE_GROUP}")
    print(f"Selection threshold: percent_covered >= {MIN_PERCENT_COVERED}, mean_depth >= {MIN_MEAN_DEPTH}")
    print(f"Max samples: {MAX_SAMPLES}")
    print(f"Samples selected: {len(samples)}")
    print(f"Output: {OUTDIR}")

    with SAMPLE_LIST.open("w", newline="") as handle:
        fields = ["sample", "sequencing_folder", "site", "sample_type", "percent_covered", "mean_depth", "max_depth"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in samples:
            writer.writerow({field: row[field] for field in fields})

    qc_fields = [
        "sample", "target_id", "length", "covered_bases", "percent_covered_from_pileup",
        "mean_depth_from_pileup", "median_depth_from_pileup", "N_count", "N_percent",
        "mixed_positions", "max_minor_allele_fraction", "status"
    ]
    mixed_fields = [
        "sample", "target_id", "reference", "ref_pos_1based", "target_pos_1based",
        "depth", "major_base", "major_count", "minor_base", "minor_count",
        "minor_allele_fraction", "A", "C", "G", "T"
    ]

    status_counts = Counter()
    with DOMINANT_FASTA.open("w") as dominant_out, \
            IUPAC_FASTA.open("w") as iupac_out, \
            QC_TABLE.open("w", newline="") as qc_handle, \
            MIXED_TABLE.open("w", newline="") as mixed_handle:

        qc_writer = csv.DictWriter(qc_handle, delimiter="\t", fieldnames=qc_fields)
        mixed_writer = csv.DictWriter(mixed_handle, delimiter="\t", fieldnames=mixed_fields)
        qc_writer.writeheader()
        mixed_writer.writeheader()

        for index, row in enumerate(samples, start=1):
            sample = row["sample"]
            print(f"[{index}/{len(samples)}] {sample}")
            try:
                dominant_seq, iupac_seq, qc, mixed_rows = call_consensus(sample, region)
            except Exception as exc:
                qc = {
                    "sample": sample,
                    "target_id": TARGET_ID,
                    "length": region["target_length"],
                    "covered_bases": 0,
                    "percent_covered_from_pileup": 0,
                    "mean_depth_from_pileup": 0,
                    "median_depth_from_pileup": 0,
                    "N_count": region["target_length"],
                    "N_percent": 100,
                    "mixed_positions": 0,
                    "max_minor_allele_fraction": 0,
                    "status": f"error:{type(exc).__name__}",
                }
                dominant_seq = "N" * region["target_length"]
                iupac_seq = dominant_seq
                mixed_rows = []
                (LOGDIR / f"{sample}.error.txt").write_text(str(exc) + "\n")

            status_counts[qc["status"]] += 1
            qc_writer.writerow({
                key: f"{qc[key]:.4f}" if isinstance(qc[key], float) else qc[key]
                for key in qc_fields
            })
            for mixed_row in mixed_rows:
                mixed_writer.writerow({
                    key: f"{mixed_row[key]:.4f}" if isinstance(mixed_row[key], float) else mixed_row[key]
                    for key in mixed_fields
                })

            header = f"{sample}|{qc['status']}|mixed_positions={qc['mixed_positions']}|N_percent={qc['N_percent']:.4f}"
            dominant_out.write(f">{header}\n{wrap(dominant_seq)}\n")
            iupac_out.write(f">{header}\n{wrap(iupac_seq)}\n")

    README.write_text(
        "\n".join([
            "Multi-copy-aware target sequence extraction V2",
            f"Target: {TARGET_ID}",
            f"Region: {region['original_reference']}:{region['bed_start'] + 1}-{region['bed_end']} ({region['strand']})",
            "",
            "This step does not assume one haploid gene copy per sample.",
            "For each target position, it records the dominant base and flags mixed positions",
            f"where minor_count >= {MIXED_MINOR_COUNT}, depth >= {MIXED_MIN_DEPTH}, and minor_allele_fraction >= {MIXED_MINOR_AF}.",
            "",
            "dominant_consensus.fasta contains the major-base consensus.",
            "iupac_ambiguity_consensus.fasta contains IUPAC ambiguity codes at mixed positions.",
            "sample_consensus_qc.tsv classifies samples as pass_single_dominant, pass_mixed_possible_multicopy, or fail_missing_data.",
            "mixed_sites.tsv records position-level mixed-base evidence.",
            "",
            "Status counts:",
            *[f"{status}\t{count}" for status, count in sorted(status_counts.items())],
            "",
        ])
    )

    print()
    print("Done.")
    print("Status counts:")
    for status, count in sorted(status_counts.items()):
        print(f"{status}\t{count}")
    print()
    print(DOMINANT_FASTA)
    print(IUPAC_FASTA)
    print(QC_TABLE)
    print(MIXED_TABLE)


if __name__ == "__main__":
    main()
