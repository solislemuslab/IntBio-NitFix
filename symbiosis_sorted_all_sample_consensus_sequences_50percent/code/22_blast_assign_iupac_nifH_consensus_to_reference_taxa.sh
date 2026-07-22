#!/usr/bin/env bash
set -euo pipefail

# Step 22: BLAST-assign mixed-IUPAC nifH consensus sequences to closest reference taxa.
#
# Purpose:
#   Assign each mixed-aware sample-derived nifH consensus sequence to the closest
#   known/reference nifH sequence from nifH.fa.aln. This is for annotating the
#   mixed-IUPAC sample tree. Because mixed-IUPAC sequences may contain ambiguity
#   codes, these labels should be interpreted as closest-reference labels, not
#   confirmed organism IDs.

module load blast+/2.17.0

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
REF_ALN="${REF_ALN:-$OUT/reference_nifH_tree/nifH.fa.aln}"
QUERY_FASTA="${QUERY_FASTA:-$OUT/nifH_consensus_tree_v2/nifH_consensus50_iupac_all_pass_pct80_depth10_Nle20.fasta}"
WORK="${WORK:-$OUT/nifH_reference_assignment_iupac_tree_v2}"
THREADS="${THREADS:-8}"

BLAST_DIR="$WORK/blast"
mkdir -p "$BLAST_DIR"

REF_UNGAPPED="$BLAST_DIR/nifH_reference_ungapped.fasta"
DB_PREFIX="$BLAST_DIR/nifH_reference_ungapped_db"
BLAST_OUT="$BLAST_DIR/iupac_nifH_vs_reference_blastn.tsv"
BEST_OUT="$BLAST_DIR/iupac_nifH_best_reference_hit_with_taxon.tsv"
SUMMARY_OUT="$BLAST_DIR/iupac_nifH_blast_assignment_summary.tsv"

if [[ ! -s "$REF_ALN" ]]; then
  echo "ERROR: missing reference alignment: $REF_ALN" >&2
  exit 1
fi

if [[ ! -s "$QUERY_FASTA" ]]; then
  echo "ERROR: missing mixed-IUPAC query FASTA: $QUERY_FASTA" >&2
  exit 1
fi

if ! command -v makeblastdb >/dev/null 2>&1 || ! command -v blastn >/dev/null 2>&1; then
  echo "ERROR: BLAST+ commands not available. Try: module load blast+/2.17.0" >&2
  exit 1
fi

echo "BLAST assignment for mixed-IUPAC nifH consensus sequences"
echo "OUT:         $OUT"
echo "Reference:   $REF_ALN"
echo "Query FASTA: $QUERY_FASTA"
echo "Output:      $BLAST_DIR"
echo "Threads:     $THREADS"
echo

# Remove alignment gaps from the reference alignment before BLAST database creation.
awk '
  /^>/ {
    if (name != "") print seq
    name=$0
    print name
    seq=""
    next
  }
  {
    gsub(/-/, "", $0)
    seq=seq $0
  }
  END {
    if (name != "") print seq
  }
' "$REF_ALN" > "$REF_UNGAPPED"

echo "Reference sequences: $(grep -c '^>' "$REF_UNGAPPED")"
echo "Query sequences:     $(grep -c '^>' "$QUERY_FASTA")"

makeblastdb -in "$REF_UNGAPPED" -dbtype nucl -out "$DB_PREFIX" >/dev/null

blastn \
  -query "$QUERY_FASTA" \
  -db "$DB_PREFIX" \
  -task blastn \
  -dust no \
  -soft_masking false \
  -max_target_seqs 10 \
  -num_threads "$THREADS" \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
  > "$BLAST_OUT"

python3 - <<PY
from pathlib import Path
from collections import OrderedDict, Counter

blast_path = Path("$BLAST_OUT")
best_path = Path("$BEST_OUT")
summary_path = Path("$SUMMARY_OUT")

cols = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
]

best = OrderedDict()
with blast_path.open() as f:
    for line in f:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        rec = dict(zip(cols, parts))
        q = rec["qseqid"]
        pident = float(rec["pident"])
        aln_len = int(rec["length"])
        bitscore = float(rec["bitscore"])
        evalue = float(rec["evalue"].replace("e", "E")) if rec["evalue"] != "0.0" else 0.0
        qlen = int(rec["qlen"])
        slen = int(rec["slen"])
        qcov = 100.0 * aln_len / qlen if qlen else 0.0
        scov = 100.0 * aln_len / slen if slen else 0.0
        score_tuple = (bitscore, pident, qcov, scov, -evalue)
        if q not in best or score_tuple > best[q][0]:
            best[q] = (score_tuple, rec, qcov, scov)

def taxon_from_subject(sseqid):
    parts = sseqid.split("|")
    if len(parts) >= 3:
        return parts[1]
    return sseqid

with best_path.open("w") as out:
    out.write("sample_id\tbest_reference_id\tbest_reference_taxon\tpercent_identity\talignment_length\tquery_coverage_percent\treference_coverage_percent\tevalue\tbitscore\n")
    for q, (_score, rec, qcov, scov) in best.items():
        out.write("\t".join([
            q,
            rec["sseqid"],
            taxon_from_subject(rec["sseqid"]),
            rec["pident"],
            rec["length"],
            f"{qcov:.4f}",
            f"{scov:.4f}",
            rec["evalue"],
            rec["bitscore"],
        ]) + "\n")

pidents = []
qcovs = []
scovs = []
taxa = Counter()
for _q, (_score, rec, qcov, scov) in best.items():
    pidents.append(float(rec["pident"]))
    qcovs.append(qcov)
    scovs.append(scov)
    taxa[taxon_from_subject(rec["sseqid"])] += 1

def safe_mean(x):
    return sum(x) / len(x) if x else 0

with summary_path.open("w") as out:
    out.write("metric\tvalue\n")
    out.write(f"assigned_sequences\t{len(best)}\n")
    out.write(f"unique_assigned_taxa\t{len(taxa)}\n")
    out.write(f"mean_percent_identity\t{safe_mean(pidents):.4f}\n")
    out.write(f"min_percent_identity\t{min(pidents) if pidents else 0:.4f}\n")
    out.write(f"max_percent_identity\t{max(pidents) if pidents else 0:.4f}\n")
    out.write(f"mean_query_coverage_percent\t{safe_mean(qcovs):.4f}\n")
    out.write(f"min_query_coverage_percent\t{min(qcovs) if qcovs else 0:.4f}\n")
    out.write(f"max_query_coverage_percent\t{max(qcovs) if qcovs else 0:.4f}\n")
    out.write(f"mean_reference_coverage_percent\t{safe_mean(scovs):.4f}\n")
    out.write(f"min_reference_coverage_percent\t{min(scovs) if scovs else 0:.4f}\n")
    out.write(f"max_reference_coverage_percent\t{max(scovs) if scovs else 0:.4f}\n")

print("Wrote best-hit table:", best_path)
print("Wrote summary:", summary_path)
print("Assigned sequences:", len(best))
print("Top taxa:")
for taxon, count in taxa.most_common(20):
    print(count, taxon)
PY

echo
echo "Done."
echo "Best-hit table: $BEST_OUT"
echo "Summary:        $SUMMARY_OUT"
