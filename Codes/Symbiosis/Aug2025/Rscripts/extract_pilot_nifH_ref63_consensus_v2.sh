#!/bin/bash

set -uo pipefail

module load SolisLemus-BioPhylo/2026.04.20

BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"

REF="$BASE/symbiosis_islands.fasta"
BAM_DIR="$BASE/symbiosis_mapped_full"
COVERAGE="$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv"

OUTDIR="$BASE/nifH_ref63_pilot_consensus_v2"
SEQDIR="$OUTDIR/sequences"
LOGDIR="$OUTDIR/logs"
SAMPLE_LIST="$OUTDIR/nifH_ref63_good_nodule_samples.txt"
COMBINED="$OUTDIR/nifH_ref63_consensus_all_samples.fasta"
QC="$OUTDIR/nifH_ref63_consensus_qc.tsv"

TARGET_ID='nifH|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island|ref63'
TARGET_REF='NZ_CP141049_Symbiosis_Island_Microvirga_lotononidis_strain_HAMBI_3237_plasmid_unnamed1'
REGION="${TARGET_REF}:1-894"

mkdir -p "$SEQDIR" "$LOGDIR"

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available."
    exit 1
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "ERROR: bcftools is not available."
    exit 1
fi

echo "Selecting good nodule samples for:"
echo "$TARGET_ID"

awk -F'\t' -v target="$TARGET_ID" '
    NR > 1 &&
    $2 == "No" &&
    $3 == "no" &&
    $8 == target &&
    $12 >= 80 &&
    $13 >= 10 {
        print $1
    }
' "$COVERAGE" | sort -u > "$SAMPLE_LIST"

echo "Selected samples:"
wc -l "$SAMPLE_LIST"

rm -f "$COMBINED"
echo -e "sample\tlength\tN_count\tN_percent\tACGT_count" > "$QC"

while read -r SAMPLE; do

    BAM="$BAM_DIR/${SAMPLE}.bam"

    SAMPLE_DIR="$LOGDIR/$SAMPLE"
    mkdir -p "$SAMPLE_DIR"

    VCF_GZ="$SAMPLE_DIR/${SAMPLE}.nifH_ref63.vcf.gz"
    MASK_BED="$SAMPLE_DIR/${SAMPLE}.nifH_ref63.zero_depth.bed"
    CONSENSUS_FULL="$SAMPLE_DIR/${SAMPLE}.consensus.full_reference.fasta"
    OUT_FASTA="$SEQDIR/${SAMPLE}.nifH_ref63.fasta"

    if [[ ! -s "$BAM" ]]; then
        echo "Missing BAM, skipping: $SAMPLE"
        continue
    fi

    echo "Extracting consensus: $SAMPLE"

    bcftools mpileup -Ou -f "$REF" -r "$REGION" "$BAM" 2> "$SAMPLE_DIR/mpileup.log" \
        | bcftools call -mv --ploidy 1 -Oz -o "$VCF_GZ" 2> "$SAMPLE_DIR/call.log"

    bcftools index -f "$VCF_GZ"

    samtools depth -a -r "$REGION" "$BAM" \
        | awk '$3 == 0 {print $1 "\t" $2-1 "\t" $2}' > "$MASK_BED"

    bcftools consensus -f "$REF" -m "$MASK_BED" "$VCF_GZ" > "$CONSENSUS_FULL" 2> "$SAMPLE_DIR/consensus.log"

    awk -v ref="$TARGET_REF" -v sample="$SAMPLE" '
        BEGIN {printing=0; seq=""}
        /^>/ {
            if (printing == 1) {
                printing=0
            }
            header=$0
            sub(/^>/, "", header)
            split(header, a, " ")
            if (a[1] == ref) {
                printing=1
                next
            }
        }
        printing == 1 {
            seq = seq toupper($0)
        }
        END {
            target = substr(seq, 1, 894)
            print ">" sample
            print target
        }
    ' "$CONSENSUS_FULL" > "$OUT_FASTA"

    cat "$OUT_FASTA" >> "$COMBINED"

    SEQ=$(awk '!/^>/ {printf "%s", toupper($0)}' "$OUT_FASTA")
    LEN=${#SEQ}
    NCOUNT=$(echo "$SEQ" | tr -cd 'Nn' | wc -c)
    ACGT=$(echo "$SEQ" | tr -cd 'ACGTacgt' | wc -c)

    if [[ "$LEN" -gt 0 ]]; then
        NPCT=$(awk -v n="$NCOUNT" -v l="$LEN" 'BEGIN {printf "%.4f", (n/l)*100}')
    else
        NPCT="NA"
    fi

    echo -e "${SAMPLE}\t${LEN}\t${NCOUNT}\t${NPCT}\t${ACGT}" >> "$QC"

done < "$SAMPLE_LIST"

echo
echo "Revised pilot consensus extraction complete."
echo "Sample list:      $SAMPLE_LIST"
echo "Sequence folder:  $SEQDIR"
echo "Combined FASTA:   $COMBINED"
echo "QC table:         $QC"
