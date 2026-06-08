# nifH ref63 Pilot Consensus Sequences

This folder contains the pilot consensus-sequence results for the selected `nifH` target:

```text
nifH|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island|ref63
```

This target was selected from the blank-aware ranking because it had strong nodule (`No`) support and low BLAN-control background. The purpose of this pilot step was to test whether a strong `nif` target could produce usable sample-level consensus sequences for alignment and phylogenetic tree construction.

## Script Used

This step was run on the cluster with:

```bash
BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"

bash "$BASE/Rscripts/extract_pilot_nifH_ref63_consensus_v2.sh"
```

## Input

The script used:

```text
symbiosis_islands.fasta
symbiosis_mapped_full/*.bam
nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv
nif_nod_original_reference_regions/nif_nod_matches_in_original_reference.tsv
```

## Method Summary

The script selected nodule samples with good coverage for the `nifH ref63` target, called variants with `bcftools`, applied variants with `bcftools consensus`, masked zero-depth positions as `N`, and wrote consensus FASTA and QC outputs.

The selected target length was 894 bp.

## Main Result

The pilot extraction recovered 159 full-length consensus sequences of 894 bp. After filtering sequences with `N_percent <= 5`, 104 high-quality consensus sequences remained for alignment and pilot tree construction.

## Files in This Folder

| File | Description | Why it is important |
|---|---|---|
| `nifH_ref63_good_nodule_samples.txt` | List of the 159 nodule samples selected for consensus extraction based on coverage. | Documents which samples were included in the pilot consensus step. |
| `nifH_ref63_consensus_all_samples.fasta` | FASTA file containing all 159 consensus sequences before N filtering. | Complete pilot consensus output before quality filtering. |
| `nifH_ref63_consensus_qc.tsv` | QC table for the 159 consensus sequences, including length, number of `N` bases, percent `N`, and A/C/G/T count. | Shows which sequences passed or failed the N-content filter. |
| `nifH_ref63_consensus_Nle5.fasta` | FASTA file containing the 104 consensus sequences with `N_percent <= 5`. | Final filtered sequence set used for MAFFT alignment and pilot tree construction. |

## QC Summary

| Metric | Value |
|---|---:|
| Selected nodule samples | 159 |
| Consensus sequences recovered | 159 |
| Consensus length | 894 bp |
| Full-length consensus sequences | 159 |
| Sequences with `N_percent <= 5` | 104 |
| Mean `N_percent` across all consensus sequences | 4.79% |

## Notes

The file `nifH_ref63_consensus_Nle5.fasta` is the final filtered FASTA file for this pilot step. The unfiltered FASTA is retained so that filtering decisions remain transparent and reproducible.

## Short Report Text

For the pilot consensus step, the selected `nifH ref63` target produced 159 full-length 894 bp consensus sequences from nodule samples. After filtering for sequences with no more than 5% ambiguous bases (`N_percent <= 5`), 104 high-quality consensus sequences remained for alignment and pilot phylogenetic tree construction.




# nifH ref63 Pilot Consensus, Alignment, and Tree

This folder contains the pilot consensus-sequence, alignment, and phylogenetic-tree results for the selected `nifH` target:

```text
nifH|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island|ref63
```

This target was selected from the blank-aware ranking because it had strong nodule (`No`) support and low BLAN-control background. The purpose of this pilot step was to test whether a strong `nif` target could produce usable sample-level consensus sequences, a multiple sequence alignment, and a pilot maximum-likelihood phylogenetic tree.

## Scripts and Commands Used

Consensus extraction was run on the cluster with:

```bash
BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"

bash "$BASE/Rscripts/extract_pilot_nifH_ref63_consensus_v2.sh"
```

The filtered consensus sequences were then aligned with MAFFT:

```bash
FILTERED="$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_Nle5.fasta"
OUTDIR="$BASE/nifH_ref63_pilot_consensus_v2"

mafft --auto "$FILTERED" > "$OUTDIR/nifH_ref63_consensus_Nle5.mafft.fasta"
```

The pilot tree was inferred with IQ-TREE:

```bash
ALIGN="$OUTDIR/nifH_ref63_consensus_Nle5.mafft.fasta"

iqtree -s "$ALIGN" -m MFP -B 1000 -T AUTO \
  --prefix "$OUTDIR/nifH_ref63_Nle5_iqtree"
```

## Input

The script used:

```text
symbiosis_islands.fasta
symbiosis_mapped_full/*.bam
nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv
nif_nod_original_reference_regions/nif_nod_matches_in_original_reference.tsv
```

## Method Summary

The script selected nodule samples with good coverage for the `nifH ref63` target, called variants with `bcftools`, applied variants with `bcftools consensus`, masked zero-depth positions as `N`, and wrote consensus FASTA and QC outputs.

The selected target length was 894 bp. After filtering consensus sequences by `N_percent <= 5`, MAFFT was used to align the filtered sequences. IQ-TREE was then used for model selection and maximum-likelihood tree inference with ultrafast bootstrap support.

## Main Result

The pilot extraction recovered 159 full-length consensus sequences of 894 bp. After filtering sequences with `N_percent <= 5`, 104 high-quality consensus sequences remained. The MAFFT alignment contained 104 sequences and was 894 bp long. IQ-TREE selected `GTR+F+I+G4` as the best-fit model according to BIC and produced the pilot `nifH` tree.

## Files in This Folder

| File | Description | Why it is important |
|---|---|---|
| `nifH_ref63_good_nodule_samples.txt` | List of the 159 nodule samples selected for consensus extraction based on coverage. | Documents which samples were included in the pilot consensus step. |
| `nifH_ref63_consensus_all_samples.fasta` | FASTA file containing all 159 consensus sequences before N filtering. | Complete pilot consensus output before quality filtering. |
| `nifH_ref63_consensus_qc.tsv` | QC table for the 159 consensus sequences, including length, number of `N` bases, percent `N`, and A/C/G/T count. | Shows which sequences passed or failed the N-content filter. |
| `nifH_ref63_consensus_Nle5.fasta` | FASTA file containing the 104 consensus sequences with `N_percent <= 5`. | Final filtered sequence set used for MAFFT alignment and pilot tree construction. |
| `nifH_ref63_consensus_Nle5.mafft.fasta` | MAFFT alignment of the 104 filtered consensus sequences. | Alignment used as input for IQ-TREE. |
| `nifH_ref63_Nle5_iqtree.treefile` | Final maximum-likelihood tree from IQ-TREE in Newick format. | Main pilot phylogenetic tree file; can be uploaded to iTOL. |
| `nifH_ref63_Nle5_iqtree.iqtree` | IQ-TREE report file, including model-selection and tree-inference details. | Documents the best-fit model and run settings. |
| `nifH_ref63_Nle5_iqtree.log` | IQ-TREE run log. | Useful for reproducibility and troubleshooting. |

## QC Summary

| Metric | Value |
|---|---:|
| Selected nodule samples | 159 |
| Consensus sequences recovered | 159 |
| Consensus length | 894 bp |
| Full-length consensus sequences | 159 |
| Sequences with `N_percent <= 5` | 104 |
| Mean `N_percent` across all consensus sequences | 4.79% |
| MAFFT alignment sequences | 104 |
| MAFFT alignment length | 894 bp |
| IQ-TREE best-fit model | `GTR+F+I+G4` |

## Notes

The file `nifH_ref63_consensus_Nle5.fasta` is the final filtered FASTA file for this pilot step. The file `nifH_ref63_consensus_Nle5.mafft.fasta` is the alignment used for tree inference, and `nifH_ref63_Nle5_iqtree.treefile` is the final pilot tree. The unfiltered FASTA is retained so that filtering decisions remain transparent and reproducible.

## Short Report Text

For the pilot `nifH ref63` analysis, 159 full-length 894 bp consensus sequences were recovered from nodule samples. After filtering for sequences with no more than 5% ambiguous bases (`N_percent <= 5`), 104 high-quality sequences remained. These 104 sequences were aligned with MAFFT, and IQ-TREE produced a pilot maximum-likelihood `nifH` tree using `GTR+F+I+G4` as the best-fit model.


