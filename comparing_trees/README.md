# Comparing Trees

This folder contains the GitHub-ready comparison package for the consensus50 10-gene functional-gene tree analysis.

Created: 2026-08-09 11:55

## Folder Layout

- `code/`: reproducible cluster shell/Python pipeline scripts and local R comparison scripts.
- `report/`: final organized report with relative links to figures, tables, and tree files.
- `result/`: report-ready outputs copied from `symbiosis_sorted_all_sample_consensus_sequences_50percent/result`.

## Main Report

Open:

- [`report/analysis_report_10_genes_consensus50.md`](report/analysis_report_10_genes_consensus50.md)

## Main Result Subfolders

- `result/comparative_tree_analysis/figures`: comparison figures across the 10 selected genes.
- `result/comparative_tree_analysis/tables`: coverage, consensus, alignment, tree, BLAST, and overlap summary tables.
- `result/metadata_tree_comparison`: metadata-aware figures and tables.
- `result/gene_trees_full`: per-gene strict and mixed-IUPAC tree outputs plus BLAST annotation tables.
- `result/trees`: selected nifH/nifD tree files and iTOL-exported tree figures.
- `result/reference`: reference checks and source reference files used in the report.

## Cluster Root Used For Reproducibility

`/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent`
