#!/usr/bin/env python3
import sys
import subprocess
import re
from collections import Counter

if len(sys.argv) != 14:
    sys.stderr.write(
        "Usage: 13_extract_nifH_consensus_consensus50.py "
        "SAMPLE BAM REF OUT_STRICT_FASTA OUT_IUPAC_FASTA "
        "MIN_BASE_DEPTH MAX_N_PERCENT MIN_BASE_QUAL MIN_MAP_QUAL "
        "MINOR_COUNT_THRESHOLD MINOR_FRACTION_THRESHOLD MIXED_SITE_DEPTH_THRESHOLD MIXED_POSITIONS_THRESHOLD\n"
    )
    sys.exit(2)

(
    sample,
    bam,
    ref,
    out_strict_fa,
    out_iupac_fa,
    min_base_depth,
    max_n_percent,
    min_base_qual,
    min_map_qual,
    minor_count_threshold,
    minor_fraction_threshold,
    mixed_site_depth_threshold,
    mixed_positions_threshold,
) = sys.argv[1:]

min_base_depth = int(min_base_depth)
max_n_percent = float(max_n_percent)
min_base_qual = int(min_base_qual)
min_map_qual = int(min_map_qual)
minor_count_threshold = int(minor_count_threshold)
minor_fraction_threshold = float(minor_fraction_threshold)
mixed_site_depth_threshold = int(mixed_site_depth_threshold)
mixed_positions_threshold = int(mixed_positions_threshold)

target = "nifH.fa"

IUPAC = {
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("CG"): "S",
    frozenset("AT"): "W",
    frozenset("GT"): "K",
    frozenset("AC"): "M",
    frozenset("CGT"): "B",
    frozenset("AGT"): "D",
    frozenset("ACT"): "H",
    frozenset("ACG"): "V",
    frozenset("ACGT"): "N",
}

def read_ref_seq(ref_path, target_name):
    seqs = {}
    name = None
    with open(ref_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.upper())
    if target_name not in seqs:
        raise SystemExit(f"ERROR: {target_name} not found in reference {ref_path}")
    return "".join(seqs[target_name])

def parse_bases(read_bases, ref_base):
    """Count explicit read bases from samtools mpileup.

    Important for consensus_sequences_50percent:
    - If the reference base is A/C/G/T, '.' and ',' mean the read matches that base, so count it.
    - If the reference base is N or an IUPAC ambiguity code, do NOT count '.' or ',' as a
      confident sample base. This avoids copying reference uncertainty into the sample consensus.
    - Explicit read bases A/C/G/T are counted regardless of whether the reference base is A/C/G/T,
      IUPAC, or N.
    """
    counts = Counter()
    i = 0
    while i < len(read_bases):
        c = read_bases[i]
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            i += 1
            m = re.match(r"\d+", read_bases[i:])
            if m:
                n = int(m.group(0))
                i += len(m.group(0)) + n
            continue
        if c in ".,":
            if ref_base in "ACGT":
                counts[ref_base] += 1
        elif c.upper() in "ACGT":
            counts[c.upper()] += 1
        i += 1
    return counts

def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

ref_seq = read_ref_seq(ref, target)
expected_len = len(ref_seq)

dom_seq = []
iupac_seq = []
depths = []
covered = 0
mixed_positions = 0
max_minor_fraction = 0.0

cmd = [
    "samtools", "mpileup",
    "-aa",
    "-A",
    "-Q", str(min_base_qual),
    "-q", str(min_map_qual),
    "-f", ref,
    "-r", target,
    bam,
]
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

for line in proc.stdout:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 5:
        continue
    ref_base = fields[2].upper()
    read_bases = fields[4]

    counts = parse_bases(read_bases, ref_base)
    total = sum(counts.values())
    depths.append(total)

    if total >= min_base_depth and counts:
        covered += 1
        sorted_counts = counts.most_common()
        major_base, major_count = sorted_counts[0]
        dom_seq.append(major_base)

        supported = {
            base for base, count in sorted_counts
            if count >= minor_count_threshold and count / total >= minor_fraction_threshold
        }
        supported.add(major_base)
        iupac_base = IUPAC.get(frozenset(supported), major_base) if len(supported) > 1 else major_base
        iupac_seq.append(iupac_base)

        minor_fraction = 0.0
        if len(sorted_counts) > 1:
            minor_fraction = sorted_counts[1][1] / total
        max_minor_fraction = max(max_minor_fraction, minor_fraction)

        if total >= mixed_site_depth_threshold and len(supported) > 1:
            mixed_positions += 1
    else:
        dom_seq.append("N")
        iupac_seq.append("N")

stderr = proc.stderr.read()
ret = proc.wait()
if ret != 0:
    sys.stderr.write(stderr)
    raise SystemExit(ret)

dom_seq = "".join(dom_seq)
iupac_seq = "".join(iupac_seq)

if len(dom_seq) < expected_len:
    missing = expected_len - len(dom_seq)
    dom_seq += "N" * missing
    iupac_seq += "N" * missing
elif len(dom_seq) > expected_len:
    dom_seq = dom_seq[:expected_len]
    iupac_seq = iupac_seq[:expected_len]

length = len(dom_seq)
n_count = dom_seq.count("N")
n_percent = 100 * n_count / length if length else 100.0
covered_percent = 100 * covered / length if length else 0.0
mean_depth = sum(depths) / len(depths) if depths else 0.0
median_depth = sorted(depths)[len(depths) // 2] if depths else 0

if n_percent > max_n_percent:
    status = "fail_high_N"
elif mixed_positions > mixed_positions_threshold:
    status = "pass_mixed_possible_multitemplate"
else:
    status = "pass_single_dominant"

if n_percent <= max_n_percent:
    # Strict tree: only samples without strong mixed signal.
    if status == "pass_single_dominant":
        with open(out_strict_fa, "w") as out:
            out.write(f">{sample}|nifH|consensus50|strict_single_dominant\n{wrap(dom_seq)}\n")

    # Mixed-aware tree: keep all passing samples; mixed sites are represented as IUPAC codes.
    with open(out_iupac_fa, "w") as out:
        out.write(f">{sample}|nifH|consensus50|{status}\n{wrap(iupac_seq)}\n")

print("\t".join([
    sample,
    target,
    str(length),
    str(covered),
    f"{covered_percent:.4f}",
    f"{mean_depth:.4f}",
    str(median_depth),
    str(n_count),
    f"{n_percent:.4f}",
    str(mixed_positions),
    f"{max_minor_fraction:.4f}",
    str(min_base_qual),
    str(min_map_qual),
    str(minor_count_threshold),
    f"{minor_fraction_threshold:.4f}",
    str(mixed_site_depth_threshold),
    str(mixed_positions_threshold),
    status,
]))
