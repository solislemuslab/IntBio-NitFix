# symbiosis_sorted_all_sample Flowchart

This file documents the V2 global analysis in reproducible step order. The analysis uses all available `symbiosis_sorted` sequencing folders, treats sequencing dates as batches rather than biological time points, and keeps MC-1/MC-2 as negative-control/background checks.

![Pipeline overview](./symbiosis_sorted_all_sample_pipeline.svg)

## Base Paths

Cluster working folder:

```bash
OUT="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2"
```


## Step 1. Create Global Manifest

Purpose: Inventory all original `symbiosis_sorted` FASTQ files across sequencing folders and create one clean sample manifest without modifying original data.

Input:

- Original data folders on the cluster, including July2024, May2025, August2025, and December2025 sequencing batches.
- Raw `*_R1.fq.gz` and `*_R2.fq.gz` files inside each `symbiosis_sorted` folder.

Code:

- [01_create_symbiosis_v2_manifest.sh](../code/01_create_symbiosis_v2_manifest.sh)

Command:

```bash
bash "$OUT/Rscripts_v2/01_create_symbiosis_v2_manifest.sh"
```

Result:

- 2,907 unique samples.
- Sample types: 952 No, 989 Rh, 964 Ro, 2 MC.
- Duplicate TALL-82-1 samples were detected across August2025 and December2025; December2025 copies were kept.

## Step 2. Trim Reads With fastp

Purpose: Clean adapters, low-quality bases, poly-G/poly-X artifacts, low-complexity reads, and short reads before mapping.

Input:

- Manifest from Step 1.
- Raw paired FASTQ files.

Code:

- [02_trim_symbiosis_v2_from_manifest.sh](../code/02_trim_symbiosis_v2_from_manifest.sh)

Command:

```bash
bash "$OUT/Rscripts_v2/02_trim_symbiosis_v2_from_manifest.sh"
```

Result:

- 2,907 completed paired trimmed samples.
- Trimmed P1/P2 files and fastp JSON/HTML logs were produced on the cluster.

## Step 3. Check Trimming Quality

Purpose: Confirm that trimmed reads are clean enough for mapping.

Input:

- fastp JSON files from Step 2.

Code:

- [04_summarize_fastp_qc_v2.sh](../code/04_summarize_fastp_qc_v2.sh)
- [04b_make_fastp_quality_profile_figures_v2.sh](../code/04b_make_fastp_quality_profile_figures_v2.sh)

Command:

```bash
bash "$OUT/Rscripts_v2/04_summarize_fastp_qc_v2.sh"
bash "$OUT/Rscripts_v2/04b_make_fastp_quality_profile_figures_v2.sh"
```

Result:

- [fastp_qc_overall_summary.tsv](../result/tables/fastp_qc_overall_summary.tsv)
- [fastp_qc_summary.svg](../result/figures/fastp_qc_summary.svg)
- [fastp_quality_after_trimming_all_samples_lightpurple_mean_blue.svg](../result/figures/fastp_quality_after_trimming_all_samples_lightpurple_mean_blue.svg)

Key metrics:

- 2,907 fastp JSON files parsed.
- Mean reads retained: 96.07%.
- Mean bases retained: 94.98%.
- Mean Q30 improved from 95.40% before trimming to 97.67% after trimming.

## Step 4. Map Reads to the Symbiosis Reference

Purpose: Map each cleaned sample to the functional-gene reference FASTA.

Input:

- Trimmed P1/P2 FASTQ files.
- `symbiosis_islands.fasta` provided by Ryan.

Code:

- [03_map_symbiosis_v2.sh](../code/03_map_symbiosis_v2.sh)
- [03b_map_symbiosis_v2_parallel.sh](../code/03b_map_symbiosis_v2_parallel.sh)

Command:

```bash
bash "$OUT/Rscripts_v2/03b_map_symbiosis_v2_parallel.sh"
```

Result:

- 2,907 BAM files and BAM indexes were generated on the cluster.

## Step 5. Check Mapping Quality

Purpose: Summarize mapping success across all samples.

Input:

- BAM files from Step 4.
- `samtools flagstat` outputs.

Code:

- [05_summarize_mapping_qc_v2.sh](../code/05_summarize_mapping_qc_v2.sh)

Command:

```bash
bash "$OUT/Rscripts_v2/05_summarize_mapping_qc_v2.sh"
```

Result:

- [symbiosis_mapping_overall_summary.tsv](../result/tables/symbiosis_mapping_overall_summary.tsv)
- [symbiosis_mapping_qc_histograms.svg](../result/figures/symbiosis_mapping_qc_histograms.svg)

Key metrics:

- 2,907 samples summarized.
- Mean mapped reads: 94.52%.
- Median mapped reads: 95.11%.
- Mean properly paired reads: 87.62%.

## Step 6. Extract Central nif/nod Target Regions

Purpose: Extract the central nif/nod gene target sequences from the GenBank annotation and gene list.

Input:

- `symbiosis_islands.gb` provided by Ryan.
- `symbiosis_islands_gene_list.xlsx` provided by Ryan.

Code:

- V1 target-extraction products were reused because the reference FASTA/GenBank files are unchanged for V2.

Result:

- [central_nif_nod_gene_summary.tsv](../result/tables/central_nif_nod_gene_summary.tsv)
- [all_central_nif_nod_target_records.tsv](../result/tables/all_central_nif_nod_target_records.tsv)
- 231 unique target sequences: 169 nif targets and 62 nod targets.

## Step 7. Check MC Negative Controls

Purpose: Replace the earlier BLAN-control interpretation with the correct MC mock-community negative controls.

Input:

- MC-1 and MC-2 symbiosis-sorted reads.
- Extracted central nif/nod target FASTA.

Code:

- [07_map_mc_negative_controls_to_nif_nod_targets_v2.sh](../code/07_map_mc_negative_controls_to_nif_nod_targets_v2.sh)

Command:

```bash
bash "$OUT/Rscripts_v2/07_map_mc_negative_controls_to_nif_nod_targets_v2.sh"
```

Result:

- [nif_nod_mc_negative_control_mapping_summary.tsv](../result/tables/nif_nod_mc_negative_control_mapping_summary.tsv)

Key result:

- MC-1 mapped 34.29% to nif/nod target sequences.
- MC-2 mapped 29.62% to nif/nod target sequences.
- Therefore MC is useful as a conservative background screen, not as a zero-signal expectation.

## Step 8. Match nif/nod Targets Back to Original Reference

Purpose: Confirm that extracted targets are exact subsequences of the original reference FASTA.

Input:

- Extracted central nif/nod targets.
- `symbiosis_islands.fasta` provided by Ryan.

Result:

- [nif_nod_matches_in_original_reference.tsv](../result/tables/nif_nod_matches_in_original_reference.tsv)
- [nif_nod_match_summary_by_gene.tsv](../result/tables/nif_nod_match_summary_by_gene.tsv)
- All 231 extracted targets matched the original FASTA; 233 exact reference locations were found.

## Step 9. Measure Coverage for Each Target in Each Sample

Purpose: For every sample and target region, measure covered bases, percent covered, mean depth, and max depth from the existing BAM files.

Input:

- BAM files from Step 4.
- Target coordinate table from Step 8.

Code:

- [09_calculate_one_sample_nif_nod_coverage_v2.py](../code/09_calculate_one_sample_nif_nod_coverage_v2.py)
- [09_calculate_nif_nod_region_coverage_v2.sh](../code/09_calculate_nif_nod_region_coverage_v2.sh)

Command:

```bash
bash "$OUT/Rscripts_v2/09_calculate_nif_nod_region_coverage_v2.sh"
```

Result:

- [nif_nod_region_coverage_all_samples.tsv.gz](../result/tables/nif_nod_region_coverage_all_samples.tsv.gz)
- 677,331 sample-target rows plus header.
- Full uncompressed table is 142 MB locally, so the repository copy is gzip-compressed.

## Step 10. Summarize Coverage by Sample Type and Threshold

Purpose: Summarize coverage patterns across No, Rh, Ro, and MC sample groups using several thresholds.

Input:

- Full coverage table from Step 9.

Code:

- [09b_summarize_nif_nod_coverage_v2.py](../code/09b_summarize_nif_nod_coverage_v2.py)
- [10_make_coverage_summary_figures_v2.py](../code/10_make_coverage_summary_figures_v2.py)

Result:

- [nif_nod_coverage_overall_threshold_summary.tsv](../result/tables/nif_nod_coverage_overall_threshold_summary.tsv)
- [nif_nod_gene_coverage_summary_by_sample_type_thresholds.tsv](../result/tables/nif_nod_gene_coverage_summary_by_sample_type_thresholds.tsv)
- [v2_gene_coverage_heatmap_pct80_depth10.svg](../result/figures/v2_gene_coverage_heatmap_pct80_depth10.svg)
- [v2_threshold_sensitivity_summary.svg](../result/figures/v2_threshold_sensitivity_summary.svg)

At the strict 80% coverage / 10X mean-depth threshold:

- No: 950/952 samples had at least one good target.
- Rh: 979/989 samples had at least one good target.
- Ro: 955/964 samples had at least one good target.

## Step 11. Rank Targets With MC-aware Filtering

Purpose: Select target regions that have strong biological sample signal and low MC background.

Input:

- Coverage summary tables from Step 10.

Result:

- [nif_nod_target_mc_aware_ranking_thresholds.tsv](../result/tables/nif_nod_target_mc_aware_ranking_thresholds.tsv)
- [v2_top_mc_aware_targets_pct80_depth10.tsv](../result/tables/v2_top_mc_aware_targets_pct80_depth10.tsv)
- [v2_nifH_target_good_coverage_by_sample_type_pct80_depth10.tsv](../result/tables/v2_nifH_target_good_coverage_by_sample_type_pct80_depth10.tsv)

Primary clean nifH target set:

- ref52
- ref56
- ref60
- ref62
- ref63

These five targets were selected because they had biological coverage and zero MC-good calls under the 80%/10X screen.

## Step 12. Extract Multi-copy-aware nifH Consensus Sequences

Purpose: Avoid forcing one sample into only one nifH sequence. A sample can contribute one sequence per clean nifH target if that target passes coverage/QC.

Input:

- BAM files from Step 4.
- Five clean nifH targets from Step 11.

Code:

- [11_extract_multicopy_aware_target_sequences_v2.py](../code/11_extract_multicopy_aware_target_sequences_v2.py)
- [11_extract_multicopy_aware_target_sequences_v2.sh](../code/11_extract_multicopy_aware_target_sequences_v2.sh)
- [11b_run_multicopy_aware_nifH_targets_v2.sh](../code/11b_run_multicopy_aware_nifH_targets_v2.sh)
- [12_make_nifH_clean5_tree_ready_v2.py](../code/12_make_nifH_clean5_tree_ready_v2.py)
- [12_make_nifH_clean5_tree_ready_v2.sh](../code/12_make_nifH_clean5_tree_ready_v2.sh)

Result:

- Strict single-dominant tree-ready FASTA: 1,802 sequences.
- Mixed-aware IUPAC FASTA: 4,212 sequences.
- 2,428 unique biological samples had at least one clean nifH target.

## Step 13. Align nifH Sequences With MAFFT

Purpose: Align the target-specific nifH consensus sequences base-by-base.

Input:

- Strict single-dominant FASTA from Step 12.

Code:

- [13_align_nifH_clean5_with_mafft_v2.sh](../code/13_align_nifH_clean5_with_mafft_v2.sh)

Command:

```bash
bash "$OUT/Rscripts_v2/13_align_nifH_clean5_with_mafft_v2.sh"
```

Result:

- Strict clean5 alignment: 1,802 sequences.
- Alignment length: 1,061 columns.

## Step 14. Build nifH Tree With IQ-TREE

Purpose: Infer a maximum-likelihood nifH tree.

Input:

- MAFFT alignment from Step 13.

Code:

- [14_build_nifH_clean5_tree_iqtree_v2.sh](../code/14_build_nifH_clean5_tree_iqtree_v2.sh)
- [14c_build_nifH_clean5_iupac_tree_fast_v2.sh](../code/14c_build_nifH_clean5_iupac_tree_fast_v2.sh)

Result:

- [nifH_clean5_strict_single_dominant.treefile](../result/trees/nifH_clean5_strict_single_dominant.treefile)
- [nifH_clean5_strict_single_dominant.iqtree](../result/trees/nifH_clean5_strict_single_dominant.iqtree)

Strict tree details:

- 1,802 sequences.
- 1,061 nucleotide sites.
- Best-fit model by BIC: GTR+F+R9.

## Step 15. Visualize Tree in iTOL

Purpose: Create report-ready visualizations with metadata rings.

Input:

- IQ-TREE treefile from Step 14.
- iTOL colorstrip annotation files.

Code:

- [15_make_nifH_clean5_itol_annotations_v2.sh](../code/15_make_nifH_clean5_itol_annotations_v2.sh)

Result:

- [single_dominant.svg](../result/figures/single_dominant.svg)
- [strict_itol_sample_type_colorstrip.txt](../result/itol_annotations/strict_itol_sample_type_colorstrip.txt)
- [strict_itol_target_ref_colorstrip.txt](../result/itol_annotations/strict_itol_target_ref_colorstrip.txt)

Interpretation:

- The strict nifH tree contains multiple divergent nifH lineages.
- Target-reference coloring shows that ref60 is especially distinct, while ref52 and ref63 show more overlap.
- The V2 multi-reference strategy captures more nifH diversity than the original one-reference pilot approach.
