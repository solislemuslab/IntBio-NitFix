# Central nif/nod Target Reference

This folder contains the gene-level target reference files used for the central `nif` and `nod` analysis in the August 2025 `symbiosis_sorted` pipeline.

The files were generated from Ryan's provided symbiosis-island annotation files:

- `symbiosis_islands.gb`: GenBank annotation file with gene names, gene coordinates, strand information, products, and DNA sequences.
- `symbiosis_islands_gene_list.xlsx`: Ryan's summary file listing expected symbiosis genes, including central `nif` and central `nod` genes.

This step does not create a new biological reference. Instead, it extracts the annotated `nif` and `nod` gene regions from Ryan's symbiosis-island reference so that downstream analyses can focus on gene-level coverage, consensus sequences, and phylogenetic trees.

## Why This Step Was Done

The previous mapping step aligned all trimmed `symbiosis_sorted` reads to the full `symbiosis_islands.fasta` reference. That full reference contains whole symbiosis-island regions and individual gene records. For phylogenetic analysis, we need specific functional-gene regions, not the full mixed reference. Therefore, this step uses the GenBank annotation to identify exactly which central `nif` and `nod` genes are present and extracts those gene sequences into smaller target FASTA files.

## Script Used

The target reference files were generated on the cluster with:

```bash
BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"

python3 "$BASE/Rscripts/extract_all_nif_nod_targets.py"
```

The script reads `symbiosis_islands.gb`, searches for the central `nif` and `nod` genes listed in `symbiosis_islands_gene_list.xlsx`, removes duplicate target sequences, and writes FASTA and summary tables.

## Main Result

Ryan's summary sheet listed 18 central `nif` genes and 17 central `nod` genes. From the GenBank annotation, 231 unique central `nif`/`nod` target sequences were extracted:

- 169 `nif` target sequences
- 62 `nod` target sequences
- 231 total unique target sequences

Six genes listed in the summary sheet were not found in the GenBank annotation:

```text
nifM
nifY
nodE
nodF
nodP
nodT
```

## Files in This Folder

| File | Description | Why it is important |
|---|---|---|
| `symbiosis_islands.gb` | Ryan's GenBank annotation file for the symbiosis-island reference. | Provides gene names, coordinates, strand, products, and DNA sequences used to extract target genes. |
| `symbiosis_islands_gene_list.xlsx` | Ryan's gene-list summary file. | Defines which genes are expected central `nif`, central `nod`, and accessory genes. |
| `central_nif_nod_gene_summary.tsv` | Summary table showing each expected central `nif`/`nod` gene, how many unique reference sequences were extracted, whether duplicates were removed, and whether the gene was found. | Main report-ready summary for this extraction step. |
| `central_genes_not_found_in_genbank.txt` | List of expected central genes that were not found in `symbiosis_islands.gb`. | Documents missing targets so they are not treated as failed downstream results. |
| `all_central_nif_nod_target_records.tsv` | Metadata table for all extracted target records. | Links each extracted sequence to its gene, accession, original record, location, strand, length, locus tag, and product. |
| `all_central_nif_nod_targets.fasta` | Combined FASTA file containing all extracted central `nif` and `nod` target sequences. | Main gene-only target reference for combined `nif`/`nod` analyses. |
| `all_central_nif_targets.fasta` | FASTA file containing only central `nif` target sequences. | Used when analyzing nitrogen-fixation genes separately. |
| `all_central_nod_targets.fasta` | FASTA file containing only central `nod` target sequences. | Used when analyzing nodulation genes separately. |

## Table Column Notes

### `central_nif_nod_gene_summary.tsv`

| Column | Meaning |
|---|---|
| `group` | Gene group: `nif` or `nod`. |
| `gene` | Gene name from Ryan's central gene list. |
| `unique_reference_sequences` | Number of unique annotated reference sequences extracted for that gene. For example, `nifA = 11` means 11 unique `nifA` reference sequences were found in the GenBank annotation, not 11 samples. |
| `duplicates_removed` | Number of duplicate sequences removed during target-reference construction. |
| `status` | Whether the gene was found in the GenBank annotation. |

### `all_central_nif_nod_target_records.tsv`

| Column | Meaning |
|---|---|
| `group` | Gene group: `nif` or `nod`. |
| `gene` | Gene name. |
| `accession` | Accession identifier for the source reference record. |
| `record` | Original GenBank record name. |
| `definition` | GenBank record description. |
| `location` | Annotated feature location from the GenBank file. |
| `start_1based` | Gene start coordinate using 1-based coordinates. |
| `end_1based` | Gene end coordinate using 1-based coordinates. |
| `strand` | Gene strand: `+` or `-`. |
| `length` | Extracted gene sequence length in base pairs. |
| `locus_tag` | Locus tag from the annotation, if available. |
| `product` | Annotated gene product. |

## How to Use These Files

These files are used after full-read mapping. The mapped BAM files show where reads align on Ryan's full symbiosis-island reference, while this folder defines the exact central `nif` and `nod` gene targets. Downstream steps use these target definitions to calculate gene-level coverage, identify well-supported targets, generate consensus sequences, and build pilot phylogenetic trees.

## Short Report Text

Ryan's gene-list file identified 18 central `nif` genes and 17 central `nod` genes. Using Ryan's GenBank annotation, we extracted 231 unique central `nif`/`nod` target sequences: 169 `nif` sequences and 62 `nod` sequences. Six expected genes (`nifM`, `nifY`, `nodE`, `nodF`, `nodP`, and `nodT`) were not present in the GenBank annotation, so they were recorded as missing and excluded from downstream target analysis.
