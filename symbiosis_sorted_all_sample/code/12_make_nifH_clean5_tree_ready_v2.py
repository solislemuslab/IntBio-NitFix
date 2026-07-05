#!/usr/bin/env python3

from pathlib import Path
import csv
import os


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                raw_name = line[1:].split()[0]
                name = raw_name.split("|")[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs


def main():
    out = Path(os.environ.get(
        "OUT",
        "/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2",
    ))
    base = out / "multicopy_aware_target_sequences_v2"
    tree = out / "nifH_clean5_tree_ready_v2"
    tree.mkdir(parents=True, exist_ok=True)

    clean_refs = set(os.environ.get(
        "CLEAN_REFS",
        "ref52,ref56,ref60,ref62,ref63",
    ).split(","))

    strict_fa = tree / "nifH_clean5_strict_single_dominant.fasta"
    mixed_fa = tree / "nifH_clean5_iupac_all_pass.fasta"
    meta_tsv = tree / "nifH_clean5_tree_ready_metadata.tsv"

    strict_count = 0
    mixed_count = 0

    with open(strict_fa, "w") as strict_out, open(mixed_fa, "w") as mixed_out, open(meta_tsv, "w", newline="") as meta_out:
        writer = csv.writer(meta_out, delimiter="\t")
        writer.writerow([
            "tree_id",
            "sample",
            "site",
            "sample_type",
            "target_id",
            "target_ref",
            "status",
            "length",
            "covered_bases",
            "percent_covered",
            "mean_depth",
            "median_depth",
            "N_count",
            "N_percent",
            "mixed_positions",
            "max_minor_allele_fraction",
            "sequence_set",
        ])

        for target_dir in sorted(base.glob("nifH_*")):
            qc = target_dir / "sample_consensus_qc.tsv"
            dom = target_dir / "dominant_consensus.fasta"
            iupac = target_dir / "iupac_ambiguity_consensus.fasta"

            if not qc.exists() or not dom.exists() or not iupac.exists():
                continue

            dominant_seqs = read_fasta(dom)
            iupac_seqs = read_fasta(iupac)

            with open(qc) as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    sample = row["sample"]
                    target_id = row["target_id"]
                    status = row["status"]
                    target_ref = target_id.split("|")[-1]

                    if target_ref not in clean_refs:
                        continue
                    if status not in {"pass_single_dominant", "pass_mixed_possible_multicopy"}:
                        continue

                    parts = sample.split("-")
                    site = parts[0]
                    sample_type = parts[-1]
                    tree_id = f"{sample}|{target_ref}|{status}"

                    common_meta = [
                        tree_id,
                        sample,
                        site,
                        sample_type,
                        target_id,
                        target_ref,
                        status,
                        row["length"],
                        row["covered_bases"],
                        row["percent_covered_from_pileup"],
                        row["mean_depth_from_pileup"],
                        row["median_depth_from_pileup"],
                        row["N_count"],
                        row["N_percent"],
                        row["mixed_positions"],
                        row["max_minor_allele_fraction"],
                    ]

                    if status == "pass_single_dominant":
                        seq = dominant_seqs.get(sample)
                        if seq:
                            strict_out.write(f">{tree_id}\n{seq}\n")
                            writer.writerow(common_meta + ["strict_single_dominant"])
                            strict_count += 1

                    seq = iupac_seqs.get(sample)
                    if seq:
                        mixed_out.write(f">{tree_id}\n{seq}\n")
                        writer.writerow(common_meta + ["iupac_all_pass"])
                        mixed_count += 1

    print("Tree-ready nifH clean5 files written")
    print(f"Clean refs: {', '.join(sorted(clean_refs))}")
    print(f"Strict FASTA: {strict_fa}")
    print(f"Mixed-aware FASTA: {mixed_fa}")
    print(f"Metadata TSV: {meta_tsv}")
    print(f"Strict single-dominant sequences: {strict_count}")
    print(f"Mixed-aware IUPAC sequences: {mixed_count}")


if __name__ == "__main__":
    main()
