# nif/nod Original Reference Regions

This folder contains the files that connect the extracted central `nif`/`nod` target sequences back to Ryan's original `symbiosis_islands.fasta` reference.

This step is important because the full sample BAM files were mapped to `symbiosis_islands.fasta`, not directly to the extracted target FASTA. Therefore, before calculating gene-level coverage, each extracted `nif`/`nod` target sequence had to be matched back to its exact location in the original reference.

## Step Summary

The extracted target file:

```text
nif_nod_target_reference/all_central_nif_nod_targets.fasta
```

was searched against Ryan's original reference:

```text
symbiosis_islands.fasta
```

This produced coordinate files showing where each central `nif`/`nod` target occurs in the original reference.

## Script Used

This step was run on the cluster with:

```bash
BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"

python3 "$BASE/Rscripts/match_nif_nod_targets_to_original_reference.py"
```

## Main Result

All 231 extracted central `nif`/`nod` target sequences matched Ryan's original `symbiosis_islands.fasta` reference.

The matching step identified **233 exact target locations** in the original reference. This number is slightly larger than 231 because some extracted target sequences occurred in more than one exact location in the original reference, especially `nifQ` and `nifW`.

The file `nif_nod_targets_without_exact_match.tsv` contains only the header line, confirming that no extracted targets were missing from the original FASTA reference.

## Files in This Folder

| File | Description | Why it is important |
|---|---|---|
| `nif_nod_matches_in_original_reference.tsv` | Main matched-coordinate table. Each row shows where an extracted `nif`/`nod` target sequence matched in Ryan's original `symbiosis_islands.fasta`. | Used by the coverage step to measure read support for each `nif`/`nod` gene region from the existing BAM files. |
| `nif_nod_matches_in_original_reference.bed` | BED-format coordinate file for the same matched regions. | Useful for genome-coordinate tools and interval-based analyses. |
| `nif_nod_match_summary_by_gene.tsv` | Gene-level summary showing how many extracted targets matched the original FASTA and how many exact locations were found. | Report-ready summary of whether each gene target was successfully matched. |
| `nif_nod_targets_without_exact_match.tsv` | List of extracted target sequences that did not match the original FASTA. This file has only a header because all targets matched. | Documents that no extracted `nif`/`nod` targets were lost during matching. |

## Column Descriptions

### `nif_nod_matches_in_original_reference.tsv`

| Column | Meaning |
|---|---|
| `original_reference` | Name of the original reference sequence in `symbiosis_islands.fasta` where the target was found. |
| `bed_start` | Start coordinate of the match using BED-style 0-based coordinates. |
| `bed_end` | End coordinate of the match. |
| `gene` | Gene name, such as `nifH`, `nifD`, or `nodC`. |
| `extracted_target_id` | Identifier of the extracted target sequence from `all_central_nif_nod_targets.fasta`. |
| `strand` | Strand where the target matched: `+` or `-`. |
| `target_length` | Length of the target sequence in base pairs. |

### `nif_nod_match_summary_by_gene.tsv`

| Column | Meaning |
|---|---|
| `gene` | Gene name. |
| `extracted_target_sequences` | Number of extracted target sequences for that gene. |
| `targets_matched_to_original_fasta` | Number of extracted target sequences that matched Ryan's original FASTA. |
| `exact_locations_in_original_fasta` | Number of exact locations found in the original FASTA. This can be larger than the number of target sequences if a target appears in multiple locations. |
| `status` | Match status for the gene. |

## How These Files Are Used Next

These matched-coordinate files are used in the next coverage step. The mapped BAM files use coordinates from `symbiosis_islands.fasta`; this folder tells the pipeline which coordinates correspond to central `nif`/`nod` target regions.

In short:

```text
extracted nif/nod target sequences
        ->
matched back to original symbiosis_islands.fasta coordinates
        ->
used with existing BAM files to calculate gene-region coverage
```

## Short Report Text

All 231 extracted central `nif`/`nod` target sequences matched Ryan's original `symbiosis_islands.fasta` reference. The matching step identified 233 exact target locations because some targets occurred in more than one original-reference location. These coordinates were then used to calculate gene-level coverage from the existing BAM files.

