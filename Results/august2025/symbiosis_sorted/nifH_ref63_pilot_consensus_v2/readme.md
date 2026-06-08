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

