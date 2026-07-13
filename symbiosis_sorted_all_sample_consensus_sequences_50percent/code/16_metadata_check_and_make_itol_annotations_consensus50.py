#!/usr/bin/env python3
"""
Check IntBio metadata against the consensus_sequences_50percent analysis samples
and create iTOL annotation files for nifH consensus trees.

Run locally after copying cluster result tables into:
  symbiosis_sorted_all_sample_consensus_sequences_50percent/result/tables
"""

from pathlib import Path
from collections import Counter, defaultdict
import argparse
import csv
import difflib
import hashlib
import re


def norm_id(x):
    return re.sub(r"[^A-Z0-9]", "", (x or "").upper())


def sample_type_from_id(sample):
    if sample in {"MC-1", "MC-2"} or sample.startswith("MC-"):
        return "MC"
    if sample.endswith("-No"):
        return "No"
    if sample.endswith("-Rh"):
        return "Rh"
    if sample.endswith("-Ro"):
        return "Ro"
    return "unknown"


def site_from_id(sample):
    if sample.startswith("MC-"):
        return "MC"
    return sample.split("-", 1)[0]


def read_delimited(path, delimiter=None):
    if not path.exists():
        return [], []
    if delimiter is None:
        delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        rows = list(reader)
        return reader.fieldnames or [], rows


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def detect_column(cols, candidates, contains_all=None, contains_any=None):
    lower_map = {c.lower().strip(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    if contains_all or contains_any:
        for c in cols:
            cl = c.lower()
            ok_all = all(s in cl for s in (contains_all or []))
            ok_any = True if not contains_any else any(s in cl for s in contains_any)
            if ok_all and ok_any:
                return c
    return None


def color_for_label(label):
    fixed = {
        "MC": "#e74c3c",
        "No": "#2b8cbe",
        "Rh": "#f39c12",
        "Ro": "#7b4fa3",
        "unknown": "#999999",
        "single_dominant": "#2c7fb8",
        "mixed_possible_multitemplate": "#d95f02",
        "fail_high_N": "#999999",
    }
    if label in fixed:
        return fixed[label]
    h = hashlib.md5(label.encode()).hexdigest()
    return "#" + h[:6]


def read_fasta_headers(path):
    if not path.exists():
        return []
    headers = []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
    return headers


def write_itol_colorstrip(path, label, tip_to_value):
    path.parent.mkdir(parents=True, exist_ok=True)
    values = sorted(set(v for v in tip_to_value.values() if v not in {"", None}))
    colors = {v: color_for_label(v) for v in values}
    with path.open("w") as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write(f"DATASET_LABEL\t{label}\n")
        f.write("COLOR\t#000000\n")
        if values:
            f.write(f"LEGEND_TITLE\t{label}\n")
            f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(values)) + "\n")
            f.write("LEGEND_COLORS\t" + "\t".join(colors[v] for v in values) + "\n")
            f.write("LEGEND_LABELS\t" + "\t".join(values) + "\n")
        f.write("DATA\n")
        for tip, value in tip_to_value.items():
            if value in colors:
                f.write(f"{tip}\t{colors[value]}\t{value}\n")


def summarize_counts(rows, col, outpath, label="value"):
    counts = Counter(r.get(col, "") or "missing" for r in rows)
    out_rows = [{label: k, "count": v} for k, v in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))]
    write_tsv(outpath, out_rows, [label, "count"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--metadata", default=None, help="Metadata CSV path. Default: result/tables/intbio_metadata_draft3.csv")
    args = parser.parse_args()

    base = Path(args.base).expanduser().resolve()
    tables = base / "result" / "tables"
    ann = base / "result" / "annotations"
    seqs = base / "result" / "sequences"

    metadata_path = Path(args.metadata).expanduser().resolve() if args.metadata else tables / "intbio_metadata_draft3.csv"
    coverage_path = tables / "consensus50_gene_coverage_all_samples.tsv"
    qc_path = tables / "nifH_consensus_qc.tsv"

    cov_cols, cov_rows = read_delimited(coverage_path, "\t")
    meta_cols, meta_rows = read_delimited(metadata_path, ",")
    qc_cols, qc_rows = read_delimited(qc_path, "\t")

    if not cov_rows:
        raise SystemExit(f"Missing or empty coverage table: {coverage_path}")
    if not meta_rows:
        raise SystemExit(f"Missing or empty metadata CSV: {metadata_path}")

    meta_sample_col = detect_column(
        meta_cols,
        ["sample", "sample_id", "sampleID", "SampleID", "sample_name", "sample.name", "Sample", "sample_name_original"],
        contains_all=["sample"],
        contains_any=["id", "name"],
    )
    if meta_sample_col is None:
        meta_sample_col = meta_cols[0]

    metadata_samples = sorted({r.get(meta_sample_col, "").strip() for r in meta_rows if r.get(meta_sample_col, "").strip()})
    analysis_samples = sorted({r["sample"].strip() for r in cov_rows if r.get("sample", "").strip()})

    meta_norm = defaultdict(list)
    for s in metadata_samples:
        meta_norm[norm_id(s)].append(s)

    missing_meta_rows = []
    for s in analysis_samples:
        if s in metadata_samples:
            continue
        normalized_hits = meta_norm.get(norm_id(s), [])
        close_hits = difflib.get_close_matches(s, metadata_samples, n=3, cutoff=0.88)
        if normalized_hits:
            status = "formatting_or_case_match"
            suggestion = ";".join(normalized_hits)
        elif close_hits:
            status = "possible_typo_or_name_mismatch"
            suggestion = ";".join(close_hits)
        else:
            status = "not_found_in_metadata"
            suggestion = ""
        missing_meta_rows.append({
            "analysis_sample": s,
            "analysis_site": site_from_id(s),
            "analysis_sample_type": sample_type_from_id(s),
            "match_status": status,
            "metadata_suggestion": suggestion,
        })

    missing_analysis_rows = []
    analysis_set = set(analysis_samples)
    for s in metadata_samples:
        if s not in analysis_set:
            missing_analysis_rows.append({
                "metadata_sample": s,
                "metadata_site_guess": site_from_id(s),
                "metadata_sample_type_guess": sample_type_from_id(s),
                "status": "not_found_in_analysis_coverage_table",
            })

    write_tsv(tables / "metadata_samples_missing_from_metadata.tsv", missing_meta_rows,
              ["analysis_sample", "analysis_site", "analysis_sample_type", "match_status", "metadata_suggestion"])
    write_tsv(tables / "metadata_samples_missing_from_analysis.tsv", missing_analysis_rows,
              ["metadata_sample", "metadata_site_guess", "metadata_sample_type_guess", "status"])

    analysis_summary_rows = []
    for sample in analysis_samples:
        analysis_summary_rows.append({
            "sample": sample,
            "site": site_from_id(sample),
            "sample_type": sample_type_from_id(sample),
            "in_metadata_exact": "yes" if sample in metadata_samples else "no",
        })
    write_tsv(tables / "metadata_analysis_sample_inventory.tsv", analysis_summary_rows,
              ["sample", "site", "sample_type", "in_metadata_exact"])
    summarize_counts(analysis_summary_rows, "sample_type", tables / "metadata_analysis_samples_by_sample_type.tsv", "sample_type")
    summarize_counts(analysis_summary_rows, "site", tables / "metadata_analysis_samples_by_site.tsv", "site")

    # Metadata summaries for biologically useful columns if present.
    useful_columns = {
        "metadata_samples_by_sample_type.tsv": detect_column(meta_cols, ["sample_type", "sampleType", "type"], contains_any=["type"]),
        "metadata_samples_by_site.tsv": detect_column(meta_cols, ["site", "siteID", "NEON_site", "site_code"], contains_any=["site"]),
        "metadata_samples_by_host_genus.tsv": detect_column(meta_cols, ["host_genus", "plant_genus", "genus"], contains_any=["genus"]),
        "metadata_samples_by_host_tribe.tsv": detect_column(meta_cols, ["host_tribe", "plant_tribe", "tribe"], contains_any=["tribe"]),
    }
    for filename, col in useful_columns.items():
        if col:
            summarize_counts(meta_rows, col, tables / filename, col)

    missing_value_rows = []
    for col in meta_cols:
        missing = sum(1 for r in meta_rows if not (r.get(col, "") or "").strip())
        missing_value_rows.append({"column": col, "missing_count": missing, "total_rows": len(meta_rows), "missing_percent": f"{100*missing/len(meta_rows):.2f}"})
    write_tsv(tables / "metadata_missing_values_by_column.tsv", missing_value_rows,
              ["column", "missing_count", "total_rows", "missing_percent"])

    # Build metadata lookup by sample id.
    meta_by_sample = {r.get(meta_sample_col, "").strip(): r for r in meta_rows if r.get(meta_sample_col, "").strip()}
    meta_lookup_norm = {}
    for s, r in meta_by_sample.items():
        n = norm_id(s)
        if n not in meta_lookup_norm:
            meta_lookup_norm[n] = r

    def get_metadata_value(sample, col, fallback="unknown"):
        if not col:
            return fallback
        row = meta_by_sample.get(sample) or meta_lookup_norm.get(norm_id(sample))
        if not row:
            return fallback
        return (row.get(col, "") or fallback).strip() or fallback

    # Tree tips: prefer copied FASTA headers if available; otherwise derive from QC table.
    strict_headers = []
    iupac_headers = []
    for p in list(seqs.glob("*strict*.fasta")) + list(tables.glob("*strict*.fasta")):
        strict_headers.extend(read_fasta_headers(p))
    for p in list(seqs.glob("*iupac*.fasta")) + list(tables.glob("*iupac*.fasta")):
        iupac_headers.extend(read_fasta_headers(p))

    if qc_rows:
        if not strict_headers:
            strict_headers = [f"{r['sample']}|nifH|consensus50|strict_single_dominant" for r in qc_rows if r.get("status") == "pass_single_dominant"]
        if not iupac_headers:
            iupac_headers = []
            for r in qc_rows:
                st = r.get("status", "")
                if st in {"pass_single_dominant", "pass_mixed_possible_multitemplate"}:
                    iupac_headers.append(f"{r['sample']}|nifH|consensus50|{st.replace('pass_', '')}")

    def sample_from_tip(tip):
        return tip.split("|", 1)[0]

    def status_from_tip(tip):
        if "mixed_possible" in tip:
            return "mixed_possible_multitemplate"
        if "single_dominant" in tip:
            return "single_dominant"
        return "unknown"

    site_col = useful_columns["metadata_samples_by_site.tsv"]
    type_col = useful_columns["metadata_samples_by_sample_type.tsv"]
    genus_col = useful_columns["metadata_samples_by_host_genus.tsv"]
    tribe_col = useful_columns["metadata_samples_by_host_tribe.tsv"]

    for prefix, headers in [("strict", strict_headers), ("iupac", iupac_headers)]:
        if not headers:
            continue
        write_itol_colorstrip(ann / f"{prefix}_itol_sample_type_colorstrip.txt", "Sample type",
                              {h: sample_type_from_id(sample_from_tip(h)) for h in headers})
        write_itol_colorstrip(ann / f"{prefix}_itol_site_colorstrip.txt", "Site",
                              {h: site_from_id(sample_from_tip(h)) for h in headers})
        write_itol_colorstrip(ann / f"{prefix}_itol_consensus_status_colorstrip.txt", "Consensus status",
                              {h: status_from_tip(h) for h in headers})
        if genus_col:
            write_itol_colorstrip(ann / f"{prefix}_itol_host_genus_colorstrip.txt", "Host genus",
                                  {h: get_metadata_value(sample_from_tip(h), genus_col) for h in headers})
        if tribe_col:
            write_itol_colorstrip(ann / f"{prefix}_itol_host_tribe_colorstrip.txt", "Host tribe",
                                  {h: get_metadata_value(sample_from_tip(h), tribe_col) for h in headers})
        if type_col:
            write_itol_colorstrip(ann / f"{prefix}_itol_metadata_sample_type_colorstrip.txt", "Metadata sample type",
                                  {h: get_metadata_value(sample_from_tip(h), type_col, sample_type_from_id(sample_from_tip(h))) for h in headers})
        if site_col:
            write_itol_colorstrip(ann / f"{prefix}_itol_metadata_site_colorstrip.txt", "Metadata site",
                                  {h: get_metadata_value(sample_from_tip(h), site_col, site_from_id(sample_from_tip(h))) for h in headers})

    run_summary = [
        {"metric": "metadata_file", "value": str(metadata_path)},
        {"metric": "metadata_sample_column", "value": meta_sample_col},
        {"metric": "metadata_unique_samples", "value": str(len(metadata_samples))},
        {"metric": "analysis_unique_samples_from_coverage", "value": str(len(analysis_samples))},
        {"metric": "analysis_samples_missing_from_metadata", "value": str(len(missing_meta_rows))},
        {"metric": "metadata_samples_missing_from_analysis", "value": str(len(missing_analysis_rows))},
        {"metric": "strict_tree_annotation_tips", "value": str(len(strict_headers))},
        {"metric": "iupac_tree_annotation_tips", "value": str(len(iupac_headers))},
    ]
    write_tsv(tables / "metadata_annotation_run_summary.tsv", run_summary, ["metric", "value"])

    print("Metadata and iTOL annotation step complete.")
    print(f"Base folder: {base}")
    print(f"Metadata sample column: {meta_sample_col}")
    print(f"Analysis samples: {len(analysis_samples)}")
    print(f"Metadata samples: {len(metadata_samples)}")
    print(f"Analysis samples missing from metadata: {len(missing_meta_rows)}")
    print(f"Metadata samples missing from analysis: {len(missing_analysis_rows)}")
    print(f"Annotation folder: {ann}")
    print(f"Summary: {tables / 'metadata_annotation_run_summary.tsv'}")


if __name__ == "__main__":
    main()
