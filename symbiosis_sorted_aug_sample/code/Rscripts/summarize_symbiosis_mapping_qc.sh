#!/bin/bash

set -uo pipefail

BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"

MAPPED="$BASE/symbiosis_mapped_full"
LOGS="$BASE/symbiosis_mapping_logs_full"
MANIFEST="$BASE/symbiosis_mapping_full_manifest.tsv"
QC="$BASE/symbiosis_mapping_QC_checking"

mkdir -p "$QC"

SUMMARY="$QC/symbiosis_mapping_completion_summary.txt"
FLAGSTAT_TABLE="$QC/symbiosis_mapping_flagstat_summary.tsv"

echo "Input mapped BAM folder:"
echo "$MAPPED"
echo

echo "Input flagstat folder:"
echo "$LOGS"
echo

echo "QC output folder:"
echo "$QC"
echo

echo "===== Mapping completion summary =====" | tee "$SUMMARY"

{
echo
echo "BAM files:"
find "$MAPPED" -maxdepth 1 -type f -name '*.bam' | wc -l

echo
echo "BAM index files:"
find "$MAPPED" -maxdepth 1 -type f -name '*.bam.bai' | wc -l

echo
echo "Flagstat reports:"
find "$LOGS" -maxdepth 1 -type f -name '*.flagstat.txt' | wc -l

echo
echo "Failed samples:"
grep -c $'\tfailed' "$MANIFEST"

echo
echo "Missing P2:"
grep -c $'\tmissing_P2' "$MANIFEST"

echo
echo "Zero-size BAM files:"
find "$MAPPED" -maxdepth 1 -type f -name '*.bam' -size 0 -print | wc -l

echo
echo "Zero-size BAM index files:"
find "$MAPPED" -maxdepth 1 -type f -name '*.bam.bai' -size 0 -print | wc -l

echo
echo "BLAN control BAM files:"
find "$MAPPED" -maxdepth 1 -type f -name 'BLAN*.bam' | wc -l
} | tee -a "$SUMMARY"

cp "$MANIFEST" "$QC/symbiosis_mapping_full_manifest.tsv"

python3 - <<'PY'
from pathlib import Path
import re

base = Path("/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted")
logs = base / "symbiosis_mapping_logs_full"
qc = base / "symbiosis_mapping_QC_checking"
out = qc / "symbiosis_mapping_flagstat_summary.tsv"

with out.open("w") as f:
    f.write("sample\ttotal_reads\tmapped_reads\tmapped_percent\tproperly_paired\tproperly_paired_percent\n")

    for path in sorted(logs.glob("*.flagstat.txt")):
        sample = path.name.replace(".flagstat.txt", "")
        text = path.read_text()

        total_reads = mapped_reads = properly_paired = "NA"
        mapped_percent = properly_paired_percent = "NA"

        for line in text.splitlines():
            if " in total " in line:
                total_reads = line.split()[0]
            elif " mapped (" in line and "mate" not in line:
                mapped_reads = line.split()[0]
                m = re.search(r"\(([^%]+)%", line)
                if m:
                    mapped_percent = m.group(1)
            elif " properly paired " in line:
                properly_paired = line.split()[0]
                m = re.search(r"\(([^%]+)%", line)
                if m:
                    properly_paired_percent = m.group(1)

        f.write(f"{sample}\t{total_reads}\t{mapped_reads}\t{mapped_percent}\t{properly_paired}\t{properly_paired_percent}\n")

print(f"Wrote flagstat summary table: {out}")
PY

echo
echo "===== All-sample mapping summary ====="

awk -F'\t' '
NR>1 {
    n++;
    mapped += $4;
    paired += $6;
}
END {
    print "samples:", n;
    print "mean_mapped_percent:", mapped/n;
    print "mean_properly_paired_percent:", paired/n;
}' "$FLAGSTAT_TABLE"

echo
echo "First 10 rows:"
column -t -s $'\t' "$FLAGSTAT_TABLE" | head -10

echo
echo "Done."
echo "Mapping QC outputs saved in:"
echo "$QC"
