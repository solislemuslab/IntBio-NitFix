#!/usr/bin/env bash
set -euo pipefail

echo "Calculating coverage for consensus_sequences_50percent gene reference"

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"

REF="$OUT/reference/consensus_sequences_50percent.fasta"
BAMDIR="$OUT/consensus_mapping_full"
LOGDIR="$OUT/consensus_gene_coverage_logs"
COVDIR="$OUT/consensus_gene_coverage"
META="$OUT/metadata"

JOBS="${JOBS:-4}"
MAX_SAMPLES="${MAX_SAMPLES:-0}"

mkdir -p "$LOGDIR" "$COVDIR" "$META"

GENE_LENGTHS="$META/consensus50_reference_gene_lengths.tsv"
SAMPLE_LIST="$META/consensus50_coverage_sample_list.tsv"
RUN_LIST="$META/consensus50_coverage_sample_list_for_this_run.tsv"
OUT_TABLE="$COVDIR/consensus50_gene_coverage_all_samples.tsv"
JOBLOG="$LOGDIR/consensus50_gene_coverage_parallel.joblog"

if [[ ! -s "$REF" ]]; then
  echo "ERROR: reference FASTA not found: $REF" >&2
  exit 1
fi

if [[ ! -d "$BAMDIR" ]]; then
  echo "ERROR: BAM folder not found: $BAMDIR" >&2
  exit 1
fi

command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH" >&2; exit 1; }
command -v parallel >/dev/null 2>&1 || { echo "ERROR: GNU parallel not found in PATH" >&2; exit 1; }

python3 - "$REF" "$GENE_LENGTHS" <<'PY'
import sys

ref, out = sys.argv[1], sys.argv[2]
name = None
seq = []

with open(out, "w") as w:
    w.write("gene_ref\tgene\tlength\n")
    with open(ref) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    gene = name.replace(".fa", "")
                    w.write(f"{name}\t{gene}\t{len(''.join(seq))}\n")
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            gene = name.replace(".fa", "")
            w.write(f"{name}\t{gene}\t{len(''.join(seq))}\n")
PY

find "$BAMDIR" -maxdepth 1 -name "*.bam" ! -name "*.tmp.bam" -print \
  | sort \
  | sed 's|.*/||; s|\.bam$||' > "$SAMPLE_LIST"

sample_count=$(wc -l < "$SAMPLE_LIST" | tr -d ' ')
echo "BAM samples found: $sample_count"

if [[ "$sample_count" -eq 0 ]]; then
  echo "ERROR: no BAM files found" >&2
  exit 1
fi

if [[ "$MAX_SAMPLES" != "0" ]]; then
  head -n "$MAX_SAMPLES" "$SAMPLE_LIST" > "$RUN_LIST"
else
  cp "$SAMPLE_LIST" "$RUN_LIST"
fi

run_count=$(wc -l < "$RUN_LIST" | tr -d ' ')
echo "Samples selected for coverage: $run_count"
echo "GNU parallel jobs: $JOBS"

echo -e "sample\tgene_ref\tgene\tgene_length\tcovered_bases\tpercent_covered\tmean_depth\tmax_depth" > "$OUT_TABLE"

coverage_one_sample() {
  sample="$1"
  bam="$BAMDIR/${sample}.bam"

  if [[ ! -s "$bam" ]]; then
    echo "WARNING: missing BAM for $sample" >&2
    return 0
  fi

  samtools depth -aa "$bam" \
    | awk -v sample="$sample" -v lengths="$GENE_LENGTHS" '
      BEGIN {
        FS=OFS="\t"
        while ((getline < lengths) > 0) {
          if (NR_lengths++ == 0) continue
          gene_ref=$1
          gene=$2
          len=$3
          gene_name[gene_ref]=gene
          gene_len[gene_ref]=len
          refs[gene_ref]=1
        }
      }
      {
        ref=$1
        depth=$3
        sum_depth[ref] += depth
        if (depth > 0) covered[ref]++
        if (depth > max_depth[ref]) max_depth[ref]=depth
      }
      END {
        for (ref in refs) {
          len=gene_len[ref]
          cov=covered[ref]+0
          pct=(len>0 ? 100*cov/len : 0)
          mean=(len>0 ? sum_depth[ref]/len : 0)
          print sample, ref, gene_name[ref], len, cov, sprintf("%.4f", pct), sprintf("%.4f", mean), max_depth[ref]+0
        }
      }'
}

export BAMDIR GENE_LENGTHS
export -f coverage_one_sample

parallel --jobs "$JOBS" --joblog "$JOBLOG" coverage_one_sample :::: "$RUN_LIST" >> "$OUT_TABLE"

echo
echo "Coverage complete."
echo "Output table: $OUT_TABLE"
echo "Gene lengths:  $GENE_LENGTHS"
echo "Job log:       $JOBLOG"
echo
echo "Rows including header:"
wc -l "$OUT_TABLE"
