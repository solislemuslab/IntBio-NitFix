#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict

BASE = Path("/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted")

ORIGINAL_REF = BASE / "symbiosis_islands.fasta"
TARGET_REF = BASE / "nif_nod_target_reference" / "all_central_nif_nod_targets.fasta"

OUTDIR = BASE / "nif_nod_original_reference_regions"
OUTDIR.mkdir(exist_ok=True)

MATCH_TSV = OUTDIR / "nif_nod_matches_in_original_reference.tsv"
MATCH_BED = OUTDIR / "nif_nod_matches_in_original_reference.bed"
SUMMARY_TSV = OUTDIR / "nif_nod_match_summary_by_gene.tsv"
UNMATCHED_TSV = OUTDIR / "nif_nod_targets_without_exact_match.tsv"


def read_fasta(path):
    records = {}
    header = None
    sequence = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(sequence).upper()
                header = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

    if header is not None:
        records[header] = "".join(sequence).upper()

    return records


def reverse_complement(sequence):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1]


def find_all(sequence, target):
    positions = []
    start = 0

    while True:
        index = sequence.find(target, start)
        if index == -1:
            break
        positions.append(index)
        start = index + 1

    return positions


original_records = read_fasta(ORIGINAL_REF)
target_records = read_fasta(TARGET_REF)

matches = []
unmatched = []
gene_target_counts = defaultdict(int)
gene_matched_target_counts = defaultdict(set)
gene_location_counts = defaultdict(int)

for target_header, target_sequence in target_records.items():

    gene = target_header.split("|")[0]
    gene_target_counts[gene] += 1

    target_found = False
    rc_sequence = reverse_complement(target_sequence)

    for original_header, original_sequence in original_records.items():

        for start in find_all(original_sequence, target_sequence):
            end = start + len(target_sequence)
            matches.append([
                original_header,
                str(start),
                str(end),
                gene,
                target_header,
                "+",
                str(len(target_sequence))
            ])
            gene_matched_target_counts[gene].add(target_header)
            gene_location_counts[gene] += 1
            target_found = True

        if rc_sequence != target_sequence:
            for start in find_all(original_sequence, rc_sequence):
                end = start + len(target_sequence)
                matches.append([
                    original_header,
                    str(start),
                    str(end),
                    gene,
                    target_header,
                    "-",
                    str(len(target_sequence))
                ])
                gene_matched_target_counts[gene].add(target_header)
                gene_location_counts[gene] += 1
                target_found = True

    if not target_found:
        unmatched.append([gene, target_header, str(len(target_sequence))])


matches.sort(key=lambda row: (row[3], row[0], int(row[1]), row[4]))
unmatched.sort(key=lambda row: (row[0], row[1]))

with open(MATCH_TSV, "w") as out:
    out.write(
        "original_reference\tbed_start\tbed_end\tgene\t"
        "extracted_target_id\tstrand\ttarget_length\n"
    )
    for row in matches:
        out.write("\t".join(row) + "\n")

with open(MATCH_BED, "w") as out:
    for row in matches:
        region_name = f"{row[3]}|{row[4]}"
        out.write(
            "\t".join([
                row[0],
                row[1],
                row[2],
                region_name,
                "0",
                row[5]
            ]) + "\n"
        )

all_genes = sorted(gene_target_counts)

with open(SUMMARY_TSV, "w") as out:
    out.write(
        "gene\textracted_target_sequences\t"
        "targets_matched_to_original_fasta\t"
        "exact_locations_in_original_fasta\tstatus\n"
    )

    for gene in all_genes:
        extracted = gene_target_counts[gene]
        matched_targets = len(gene_matched_target_counts[gene])
        locations = gene_location_counts[gene]

        if matched_targets == 0:
            status = "NO_EXACT_MATCH"
        elif locations > matched_targets:
            status = "MULTIPLE_REFERENCE_LOCATIONS"
        else:
            status = "MATCHED"

        out.write(
            f"{gene}\t{extracted}\t{matched_targets}\t"
            f"{locations}\t{status}\n"
        )

with open(UNMATCHED_TSV, "w") as out:
    out.write("gene\textracted_target_id\ttarget_length\n")
    for row in unmatched:
        out.write("\t".join(row) + "\n")

print("============================================================")
print("Step 5: nif/nod target matching to original reference")
print("============================================================")
print(f"Original FASTA reference sequences:        {len(original_records)}")
print(f"Extracted nif/nod target sequences:        {len(target_records)}")
print(f"Exact target locations identified:         {len(matches)}")
print(f"Target sequences without an exact match:   {len(unmatched)}")
print()
print("Summary by gene:")
print("gene\ttargets\tmatched_targets\texact_locations\tstatus")

for gene in all_genes:
    extracted = gene_target_counts[gene]
    matched_targets = len(gene_matched_target_counts[gene])
    locations = gene_location_counts[gene]

    if matched_targets == 0:
        status = "NO_EXACT_MATCH"
    elif locations > matched_targets:
        status = "MULTIPLE_REFERENCE_LOCATIONS"
    else:
        status = "MATCHED"

    print(f"{gene}\t{extracted}\t{matched_targets}\t{locations}\t{status}")

print()
print("Output files:")
print(MATCH_TSV)
print(MATCH_BED)
print(SUMMARY_TSV)
print(UNMATCHED_TSV)
