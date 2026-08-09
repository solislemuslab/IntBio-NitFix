#!/usr/bin/env bash
set -euo pipefail

# nodL: BLAST strict and mixed-IUPAC sample consensus sequences against nodL taxon references.
# Output labels are closest known reference hits, not confirmed species identity.

OUT="${OUT:-/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent}"
GENE="nodL"
TAXON_REF_DIR="${TAXON_REF_DIR:-$OUT/reference_taxon_fastas_by_gene}"
RUN_DIR="$OUT/gene_trees/nodL/01_consensus"
BLAST_DIR="$OUT/gene_trees/nodL/05_blast_taxon_annotation"
MIN_PCT="${MIN_PCT:-80}"
MIN_MEAN_DEPTH="${MIN_MEAN_DEPTH:-10}"
MAX_N_PERCENT="${MAX_N_PERCENT:-20}"
EVALUE="${EVALUE:-1e-20}"
THREADS="${THREADS:-8}"
mkdir -p "$BLAST_DIR"

module load blast+/2.17.0 || true
if ! command -v makeblastdb >/dev/null 2>&1 || ! command -v blastn >/dev/null 2>&1; then echo "ERROR: BLAST+ not available" >&2; exit 1; fi

REF=""
for candidate in "$TAXON_REF_DIR/nodL.fa" "$TAXON_REF_DIR/nodL.fasta" "$TAXON_REF_DIR/nodL.fa.aln" "$TAXON_REF_DIR/nodL.aln"; do
  if [[ -s "$candidate" ]]; then REF="$candidate"; break; fi
done
if [[ -z "$REF" ]]; then echo "ERROR: no taxon reference found for nodL in $TAXON_REF_DIR" >&2; exit 1; fi

ungap_fasta () { awk '/^>/ {print; next} {gsub(/[-.]/, ""); print toupper($0)}' "$1"; }

best_hits_python () {
python3 - "$1" "$2" "$3" <<'PYB'
import sys, csv
blast_tsv, query_fasta, out_tsv = sys.argv[1:]
qlen={}
name=None; seq=[]
for line in open(query_fasta):
    line=line.strip()
    if not line: continue
    if line.startswith('>'):
        if name: qlen[name]=len(''.join(seq).replace('-',''))
        name=line[1:].split()[0]; seq=[]
    else:
        seq.append(line)
if name: qlen[name]=len(''.join(seq).replace('-',''))

best={}
with open(blast_tsv) as f:
    for row in csv.reader(f, delimiter='	'):
        if len(row) < 12: continue
        q,s,pident,alen,mm,gap,qs,qe,ss,se,e,bits = row[:12]
        score=(float(bits), float(pident), int(alen))
        if q not in best or score > best[q][0]:
            best[q]=(score,row)

with open(out_tsv,'w',newline='') as out:
    w=csv.writer(out, delimiter='	')
    w.writerow(['sample_id','best_reference_id','best_reference_taxon','best_reference_genus','percent_identity','alignment_length','query_coverage_percent','evalue','bitscore'])
    for q in sorted(best):
        row=best[q][1]
        s=row[1]; pident=float(row[2]); alen=int(row[3]); e=row[10]; bits=float(row[11])
        parts=s.split('|')
        taxon=parts[1] if len(parts) >= 2 else s
        genus=taxon.split('_')[0]
        qc=100*alen/qlen[q] if qlen.get(q) else 0
        w.writerow([q,s,taxon,genus,f'{pident:.3f}',alen,f'{qc:.4f}',e,f'{bits:.1f}'])
PYB
}

REF_UNGAP="$BLAST_DIR/nodL_taxon_reference_ungapped.fasta"
ungap_fasta "$REF" > "$REF_UNGAP"
makeblastdb -in "$REF_UNGAP" -dbtype nucl -out "$BLAST_DIR/nodL_taxon_ref_db" >/dev/null

for SET in strict iupac; do
  if [[ "$SET" == "strict" ]]; then
    QUERY="$RUN_DIR/nodL_consensus50_strict_single_dominant_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.fasta"
  else
    QUERY="$RUN_DIR/nodL_consensus50_iupac_all_pass_pct${MIN_PCT}_depth${MIN_MEAN_DEPTH}_Nle${MAX_N_PERCENT}.fasta"
  fi
  if [[ ! -s "$QUERY" ]]; then echo "Skipping $SET; query missing: $QUERY"; continue; fi

  QUERY_UNGAP="$BLAST_DIR/nodL_${SET}_query_ungapped.fasta"
  BLAST_TSV="$BLAST_DIR/nodL_${SET}_blast.tsv"
  BEST_TSV="$BLAST_DIR/nodL_${SET}_best_reference_hit_with_taxon.tsv"
  ungap_fasta "$QUERY" > "$QUERY_UNGAP"
  blastn -query "$QUERY_UNGAP" -db "$BLAST_DIR/nodL_taxon_ref_db"     -out "$BLAST_TSV"     -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'     -evalue "$EVALUE" -num_threads "$THREADS" -max_target_seqs 10
  best_hits_python "$BLAST_TSV" "$QUERY_UNGAP" "$BEST_TSV"
  echo "Wrote: $BEST_TSV"
done
