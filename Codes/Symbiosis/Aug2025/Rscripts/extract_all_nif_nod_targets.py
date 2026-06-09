#!/usr/bin/env python3

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path("/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted")
GB = BASE / "symbiosis_islands.gb"
OUTDIR = BASE / "nif_nod_target_reference"
OUTDIR.mkdir(parents=True, exist_ok=True)

NIF_GENES = [
    "nifA", "nifB", "nifD", "nifE", "nifH", "nifJ", "nifK", "nifM",
    "nifN", "nifQ", "nifS", "nifT", "nifU", "nifV", "nifW", "nifX",
    "nifY", "nifZ"
]

NOD_GENES = [
    "nodA", "nodB", "nodC", "nodD", "nodE", "nodF", "nodH", "nodI",
    "nodJ", "nodL", "nodP", "nodQ", "nodS", "nodT", "nodU", "nodX",
    "nodZ"
]

TARGETS = NIF_GENES + NOD_GENES
TARGET_LOOKUP = {g.lower(): g for g in TARGETS}

if not GB.exists():
    raise SystemExit(f"ERROR: GenBank file not found: {GB}")


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def clean_sequence(text):
    return re.sub(r"[^ACGTNacgtn]", "", text).upper()


def extract_location_sequence(location, record_sequence):
    is_reverse = "complement" in location.lower()
    ranges = re.findall(r"<?(\d+)\.\.>?(\d+)", location)

    if not ranges:
        return ""

    pieces = []
    for start, end in ranges:
        start_i = int(start)
        end_i = int(end)
        pieces.append(record_sequence[start_i - 1:end_i])

    sequence = "".join(pieces)

    if is_reverse:
        sequence = reverse_complement(sequence)

    return sequence


def canonical_gene_from_text(gene_text, product_text, label_text):
    gene_clean = gene_text.strip().lower()

    if gene_clean in TARGET_LOOKUP:
        return TARGET_LOOKUP[gene_clean]

    combined = f"{gene_text} {product_text} {label_text}".lower()

    # Prefer explicit gene-family patterns in product/label fields.
    for target in TARGETS:
        if re.search(rf"\b{target.lower()}\b", combined):
            return target

    # Product descriptions sometimes omit the gene-style name.
    inferred_patterns = {
        "nifU": r"nifu family protein",
        "nifQ": r"nitrogen fixation protein nifq",
        "nifZ": r"nitrogen fixation protein nifz",
        "nodA": r"noda family",
        "nodZ": r"nodulation protein nodz",
        "nodU": r"nodulation protein u\b",
    }

    for target, pattern in inferred_patterns.items():
        if re.search(pattern, combined):
            return target

    return ""


def parse_records(path):
    text = path.read_text(errors="replace")
    raw_records = re.split(r"\n//\s*\n?", text)
    extracted = []

    for raw_record in raw_records:
        if "LOCUS" not in raw_record or "ORIGIN" not in raw_record:
            continue

        locus_match = re.search(r"^LOCUS\s+(.+?)\s+\d+\s+bp", raw_record, re.MULTILINE)
        locus = locus_match.group(1).strip() if locus_match else "unknown_record"

        accession_match = re.search(r"^ACCESSION\s+(\S+)", raw_record, re.MULTILINE)
        accession = accession_match.group(1) if accession_match else locus.split()[0]

        definition_match = re.search(r"^DEFINITION\s+(.+?)(?=^ACCESSION|^VERSION)", raw_record,
                                     re.MULTILINE | re.DOTALL)
        definition = " ".join(definition_match.group(1).split()) if definition_match else ""

        before_origin, origin = raw_record.split("ORIGIN", 1)
        record_sequence = clean_sequence(origin)

        features_match = re.search(
            r"FEATURES\s+Location/Qualifiers(.*)$",
            before_origin,
            re.DOTALL
        )

        if not features_match:
            continue

        feature_chunks = re.split(r"\n(?=     \S)", features_match.group(1))

        for chunk in feature_chunks:
            first_line = chunk.splitlines()[0] if chunk.splitlines() else ""
            feature_match = re.match(r"\s{5}(\S+)\s+(.+)", first_line)

            if not feature_match:
                continue

            feature_type = feature_match.group(1)
            location = feature_match.group(2).strip()

            if feature_type != "CDS":
                continue

            gene_match = re.search(r'/gene="([^"]+)"', chunk, re.IGNORECASE)
            product_match = re.search(r'/product="([^"]+)"', chunk, re.IGNORECASE | re.DOTALL)
            label_match = re.search(r'/label="([^"]+)"', chunk, re.IGNORECASE | re.DOTALL)
            locus_tag_match = re.search(r'/locus_tag="([^"]+)"', chunk, re.IGNORECASE)

            raw_gene = gene_match.group(1).strip() if gene_match else ""
            product = " ".join(product_match.group(1).split()) if product_match else ""
            label = " ".join(label_match.group(1).split()) if label_match else ""
            locus_tag = locus_tag_match.group(1).strip() if locus_tag_match else ""

            gene = canonical_gene_from_text(raw_gene, product, label)

            if gene not in TARGETS:
                continue

            sequence = extract_location_sequence(location, record_sequence)

            if not sequence:
                continue

            group = "nif" if gene.startswith("nif") else "nod"
            strand = "-" if "complement" in location.lower() else "+"
            coords = [int(x) for x in re.findall(r"\d+", location)]

            extracted.append({
                "group": group,
                "gene": gene,
                "accession": accession,
                "record": locus,
                "definition": definition,
                "location": location,
                "start_1based": min(coords) if coords else "",
                "end_1based": max(coords) if coords else "",
                "strand": strand,
                "length": len(sequence),
                "locus_tag": locus_tag,
                "product": product,
                "sequence": sequence,
            })

    return extracted


all_records = parse_records(GB)

# Remove exact duplicate sequences within the same gene family to avoid redundant
# target sequences in later gene-specific mapping.
unique_records = []
duplicates_removed = Counter()
seen = set()

for record in all_records:
    key = (record["gene"], record["sequence"])
    if key in seen:
        duplicates_removed[record["gene"]] += 1
        continue
    seen.add(key)
    unique_records.append(record)

unique_records.sort(key=lambda x: (x["group"], x["gene"], x["accession"], x["record"]))

metadata_file = OUTDIR / "all_central_nif_nod_target_records.tsv"
summary_file = OUTDIR / "central_nif_nod_gene_summary.tsv"
combined_fasta = OUTDIR / "all_central_nif_nod_targets.fasta"
nif_fasta = OUTDIR / "all_central_nif_targets.fasta"
nod_fasta = OUTDIR / "all_central_nod_targets.fasta"
absent_file = OUTDIR / "central_genes_not_found_in_genbank.txt"

fields = [
    "group", "gene", "accession", "record", "definition", "location",
    "start_1based", "end_1based", "strand", "length", "locus_tag", "product"
]

with metadata_file.open("w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for record in unique_records:
        writer.writerow({field: record[field] for field in fields})

gene_counts = Counter(record["gene"] for record in unique_records)

with summary_file.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["group", "gene", "unique_reference_sequences", "duplicates_removed", "status"])
    for gene in TARGETS:
        group = "nif" if gene.startswith("nif") else "nod"
        count = gene_counts[gene]
        status = "found" if count > 0 else "not_found_in_GenBank_reference"
        writer.writerow([group, gene, count, duplicates_removed[gene], status])

absent_genes = [gene for gene in TARGETS if gene_counts[gene] == 0]

with absent_file.open("w") as out:
    for gene in absent_genes:
        out.write(gene + "\n")


def fasta_header(record, index):
    record_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", record["record"]).strip("_")
    return f"{record['gene']}|{record['accession']}|{record_name}|ref{index}"


with combined_fasta.open("w") as combined, \
     nif_fasta.open("w") as nif_out, \
     nod_fasta.open("w") as nod_out:

    for index, record in enumerate(unique_records, start=1):
        fasta_entry = f">{fasta_header(record, index)}\n{record['sequence']}\n"
        combined.write(fasta_entry)

        if record["group"] == "nif":
            nif_out.write(fasta_entry)
        else:
            nod_out.write(fasta_entry)

print("============================================================")
print("Extraction of central nif/nod genes from GenBank complete")
print("============================================================")
print()
print(f"Annotated CDS records extracted before deduplication: {len(all_records)}")
print(f"Unique gene reference sequences retained:             {len(unique_records)}")
print()
print("Gene summary:")
for gene in TARGETS:
    status = "FOUND" if gene_counts[gene] > 0 else "NOT FOUND"
    print(f"{gene:5s}  {gene_counts[gene]:3d} unique sequence(s)  {status}")

print()
print("Genes listed in the summary sheet but not found in this GenBank reference:")
if absent_genes:
    print(", ".join(absent_genes))
else:
    print("None")

print()
print("Outputs written to:")
print(OUTDIR)
print(metadata_file.name)
print(summary_file.name)
print(combined_fasta.name)
print(nif_fasta.name)
print(nod_fasta.name)
print(absent_file.name)
