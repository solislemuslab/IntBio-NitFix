# symbiosis_sorted_all_sample_consensus_sequences_50percent Report

## Analysis Goal

This analysis reruns the global `symbiosis_sorted` workflow using the new 50% consensus functional-gene reference provided by Pranoti. The goal is to move from many duplicated/diverse reference regions per gene toward a cleaner one-reference-per-gene strategy, then test whether the full set of 2,907 symbiosis-sorted samples can be used to recover high-quality nif/nod functional-gene signal and build a first `nifH` tree.

The main biological goal is to identify which samples contain enough read support for target nitrogen-fixation and nodulation genes, then construct tree-ready consensus sequences for `nifH` while explicitly tracking possible mixed-template signal.

## Folder Structure

| Folder | Purpose |
|---|---|
| [`code`](../code) | Scripts used for mapping, coverage calculation, nifH consensus extraction, tree building, and metadata/iTOL annotation. |
| [`result/reference`](../result/reference) | New consensus reference FASTA used in this analysis. |
| [`result/tables`](../result/tables) | Coverage summaries, metadata checks, nifH consensus QC, and other tabular outputs. |
| [`result/figures`](../result/figures) | Summary plots for gene coverage and per-sample coverage. |
| [`result/annotations`](../result/annotations) | iTOL annotation files for strict and mixed-IUPAC nifH trees. |
| [`result/trees`](../result/trees) | Intended location for copied IQ-TREE outputs. At the time of this report, tree files were still being generated/copied from the cluster. |

## Input Data

| Input | Description |
|---|---|
| Trimmed reads from V2 | The analysis reused the already trimmed FASTQ files from the previous all-sample V2 workflow: `/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2/symbiosis_trimmed_fastp`. |
| Sample set | 2,907 unique samples across all sequencing batches. These include `No`, `Rh`, `Ro`, and two mock-community controls, `MC-1` and `MC-2`. |
| New reference | [`consensus_sequences_50percent.fasta`](../result/reference/consensus_sequences_50percent.fasta), recommended by Ryan. This contains one consensus reference sequence per successfully generated functional gene. |
| Metadata | [`intbio_metadata_draft4.csv`](../result/tables/intbio_metadata_draft4.csv), used to add sample/site/host annotations to the tree. |

## Reference Summary

The new reference is [`consensus_sequences_50percent.fasta`](../result/reference/consensus_sequences_50percent.fasta). It contains 72 functional-gene reference sequences. Each sequence is a gene-level consensus built from multiple source sequences. Unlike the earlier `symbiosis_islands.fasta` reference, this file is not a collection of long symbiosis-island regions; it is a compact gene-reference set.

| Reference property | Value |
|---|---:|
| Number of reference records | 72 |
| Minimum gene length | 168 bp |
| Maximum gene length | 3,528 bp |
| Total reference length | 78,357 bp |
| A/C/G/T bases | 70,208 bp, 89.60% |
| IUPAC ambiguity bases excluding `N` | 6,724 bp, 8.58% |
| `N` bases | 1,425 bp, 1.82% |

Ryan recommended this reference because it chooses a base, including IUPAC ambiguity codes, when at least 50% of the source sequences support that call. This makes it more biologically broad than a single-isolate reference, while avoiding the very high ambiguity burden of the all-ambiguities reference.

## Step 1. Reuse Trimmed Reads

No new trimming was performed for this version. The trimmed reads from the previous all-sample V2 workflow were reused because they were already generated with `fastp` and passed trimming QC.

| Item | Path |
|---|---|
| Trimmed reads on cluster | `/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2/symbiosis_trimmed_fastp` |
| Previous trimming QC | See previous all-sample report in `symbiosis_sorted_all_sample`. |

## Step 2. Map Reads to the Consensus50 Reference

All 2,907 trimmed paired samples were mapped to the new 72-gene consensus reference using BWA-MEM, then BAM files were sorted and indexed with `samtools`.

| Item | Link |
|---|---|
| Mapping script | [`11_map_consensus50_all_samples_parallel.sh`](../code/11_map_consensus50_all_samples_parallel.sh) |
| Cluster output BAM folder | `/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent/consensus_mapping_full` |
| Cluster mapping logs | `/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent/consensus_mapping_logs` |

Important interpretation: the overall mapping percentage is lower than in the earlier symbiosis-island reference workflow because this new reference is much smaller and contains only 72 consensus gene sequences. Therefore, total mapped-read percentage is not the main QC metric here. The more meaningful metric is per-gene coverage and depth.

## Step 3. Calculate Per-Gene Coverage

Coverage was calculated for every sample-gene pair using the BAM files from Step 2. A gene was considered well covered in a sample if it met both thresholds:

- percent covered >= 80%
- mean depth >= 10

| Item | Link |
|---|---|
| Coverage script | [`12_calculate_consensus_gene_coverage_v2.sh`](../code/12_calculate_consensus_gene_coverage_v2.sh) |
| Full coverage table | [`consensus50_gene_coverage_all_samples.tsv`](../result/tables/consensus50_gene_coverage_all_samples.tsv) |
| Coverage-bin table | [`consensus50_gene_coverage_bins_depth10.tsv`](../result/tables/consensus50_gene_coverage_bins_depth10.tsv) |

### Coverage Summary

| Metric | Value |
|---|---:|
| Samples analyzed | 2,907 |
| Genes in reference | 72 |
| Samples with at least one good gene | 2,052 |
| Good-gene rule | percent covered >= 80%; mean depth >= 10 |

The gene-level coverage figure shows that `nifH` has strong recovery and is a good first target for tree construction.

![Gene coverage bins](../result/figures/consensus50_gene_coverage_bins_depth10_R.png)

## Step 4. Visualize Gene Recovery Per Sample

To understand how broad the functional-gene signal is within each sample, I counted how many of the 72 consensus-reference genes passed the good-coverage threshold in each sample.

| Item | Link |
|---|---|
| Per-sample gene-count figure | [`consensus50_genes_per_sample_pct80.png`](../result/figures/consensus50_genes_per_sample_pct80.png) |
| Per-sample gene-count by sample type | [`consensus50_covered_gene_count_by_sample_type_blocks_with_MC.png`](../result/figures/consensus50_covered_gene_count_by_sample_type_blocks_with_MC.png) |

The grouped figure is useful because it separates `No`, `Rh`, `Ro`, and `MC` samples. Gray points mark samples with zero covered genes. The MC controls are shown separately because they are mock-community controls and should not be treated as field samples.

![Genes covered per sample type](../result/figures/consensus50_covered_gene_count_by_sample_type_blocks_with_MC.png)

## Step 5. Check MC Controls

Ryan clarified that BLAN is a real NEON field site, not a blank/negative control. Therefore, this version uses only `MC-1` and `MC-2` as mock-community negative controls.

The MC controls are limited to n=2, so they are useful as a warning/background check, not as a strong statistical negative-control set.

| MC result | Value |
|---|---|
| MC samples included | 2 |
| MC samples with good gene signal | 2 |
| Gene recovered in both MC controls | `nifJ` |
| MC-1 `nifJ` coverage | 90.08%; mean depth 14.38 |
| MC-2 `nifJ` coverage | 90.05%; mean depth 10.15 |

Because the MC controls are not zero-signal, the analysis should avoid using MC as a simple presence/absence filter. Instead, MC results should be reported as a caveat and used conservatively when ranking candidate genes or interpreting background signal.

## Step 6. Build nifH Consensus Sequences

The first tree was built for `nifH` because `nifH` is a standard marker for nitrogen fixation and had strong coverage in the consensus50 reference analysis.

| Item | Link |
|---|---|
| nifH consensus extraction Python script | [`13_extract_nifH_consensus_consensus50.py`](../code/13_extract_nifH_consensus_consensus50.py) |
| nifH consensus wrapper script | [`13_run_nifH_consensus_consensus50.sh`](../code/13_run_nifH_consensus_consensus50.sh) |
| nifH consensus QC table | [`nifH_consensus_qc.tsv`](../result/tables/nifH_consensus_qc.tsv) |

### Sample Selection for nifH Consensus

A sample was selected for nifH consensus construction if `nifH` passed:

- percent covered >= 80%
- mean depth >= 10

Then per-base consensus calling was done from the BAM pileup.

### Pileup and Base-Calling Rules

| Rule | Value |
|---|---:|
| Minimum base quality | 20 |
| Minimum mapping quality | 10 |
| Maximum allowed `N` percent per consensus | 20% |
| Mixed-site minimum depth | 20 |
| Mixed-site minor allele count | >= 5 |
| Mixed-site minor allele fraction | >= 0.20 |
| Mixed-sequence rule | mixed positions > 10 |

The code avoids simply copying ambiguity from the reference. If the reference base is an IUPAC ambiguity code or `N`, the sample consensus is still called from the actual read bases observed in that sample. For example, if the reference has `N` but the reads strongly support `A`, the sample consensus becomes `A`, not `N`.

### nifH Consensus Output

| Status | Count | Meaning |
|---|---:|---|
| `pass_single_dominant` | 387 | Samples with a clear dominant base pattern and <=20% `N`. These are used in the strict tree. |
| `pass_mixed_possible_multitemplate` | 1,151 | Samples passing coverage/N filters but with more than 10 mixed positions. These are included only in the mixed-IUPAC tree. |
| `fail_high_N` | 34 | Samples with too many unknown bases after pileup consensus. These are excluded from both trees. |
| Total nifH samples evaluated | 1,572 | Samples with enough nifH coverage to attempt consensus calling. |

The strict FASTA contains 387 sequences. The mixed-IUPAC FASTA contains 1,538 sequences, which equals `pass_single_dominant` plus `pass_mixed_possible_multitemplate`.

## Step 7. Align nifH Consensus Sequences

The nifH consensus FASTA files were aligned with MAFFT before tree inference.

| Item | Cluster/local status |
|---|---|
| Strict aligned FASTA | Generated on cluster; 387 sequences. |
| Mixed-IUPAC aligned FASTA | Generated on cluster; 1,538 sequences. |
| Local tree folder | Tree outputs should be copied to [`result/trees`](../result/trees) once complete. |

The aligned sequences may contain few or no gaps because all sequences are reconstructed against the same `nifH.fa` consensus reference coordinate system and have similar length. MAFFT is still useful because it ensures all sequences are in a consistent multiple-alignment format for IQ-TREE.

## Step 8. Build nifH Trees

Two nifH trees are planned/generated:

| Tree | Input sequences | Purpose | Script |
|---|---:|---|---|
| Strict single-dominant tree | 387 | Conservative tree using only clear dominant consensus calls. | [`14_build_nifH_consensus50_strict_tree_iqtree_v2.sh`](../code/14_build_nifH_consensus50_strict_tree_iqtree_v2.sh) |
| Mixed-IUPAC tree | 1,538 | More inclusive tree allowing ambiguity codes where samples show possible mixed-template signal. | [`14b_build_nifH_consensus50_iupac_tree_iqtree_v2.sh`](../code/14b_build_nifH_consensus50_iupac_tree_iqtree_v2.sh) |

The strict tree is the most conservative result. The mixed-IUPAC tree is broader and better represents uncertainty, but it should not be interpreted as fully resolving multiple organisms or haplotypes within a sample.

At the time this report was written, local tree files were not yet present in [`result/trees`](../result/trees). Once IQ-TREE finishes on the cluster, copy the `.treefile`, `.iqtree`, `.log`, and support files into that folder and update this section.

## Step 9. Metadata Matching and iTOL Annotation

Metadata were used to annotate tree tips by sample type, site, host genus, and host tribe.

| Item | Link |
|---|---|
| Metadata file | [`intbio_metadata_draft4.csv`](../result/tables/intbio_metadata_draft4.csv) |
| Metadata annotation script | [`16_metadata_check_and_make_itol_annotations_consensus50.py`](../code/16_metadata_check_and_make_itol_annotations_consensus50.py) |
| Metadata run summary | [`metadata_annotation_run_summary.tsv`](../result/tables/metadata_annotation_run_summary.tsv) |
| Strict iTOL annotations | [`result/annotations`](../result/annotations) |
| IUPAC iTOL annotations | [`result/annotations`](../result/annotations) |

### Metadata Matching Summary

After harmonizing the `OAES_19` and `OAES_24` naming format, the corrected metadata file contained 2,898 unique sample IDs, while the consensus50 analysis contained 2,907 unique samples.

| Metadata comparison | Count |
|---|---:|
| Metadata unique samples | 2,898 |
| Analysis unique samples | 2,907 |
| Analysis samples missing from metadata | 24 |
| Metadata samples missing from analysis | 15 |
| Net difference | 9 |

The earlier count of 21 metadata-only samples came from the pre-correction metadata naming, where six OAES samples used hyphens instead of the underscore format present in the raw sequencing files and analysis outputs.

### Annotation Files

| Annotation type | Strict tree | Mixed-IUPAC tree |
|---|---|---|
| Sample type | [`strict_itol_sample_type_colorstrip.txt`](../result/annotations/strict_itol_sample_type_colorstrip.txt) | [`iupac_itol_sample_type_colorstrip.txt`](../result/annotations/iupac_itol_sample_type_colorstrip.txt) |
| Site | [`strict_itol_site_colorstrip.txt`](../result/annotations/strict_itol_site_colorstrip.txt) | [`iupac_itol_site_colorstrip.txt`](../result/annotations/iupac_itol_site_colorstrip.txt) |
| Host genus | [`strict_itol_host_genus_colorstrip.txt`](../result/annotations/strict_itol_host_genus_colorstrip.txt) | [`iupac_itol_host_genus_colorstrip.txt`](../result/annotations/iupac_itol_host_genus_colorstrip.txt) |
| Host tribe | [`strict_itol_host_tribe_colorstrip.txt`](../result/annotations/strict_itol_host_tribe_colorstrip.txt) | [`iupac_itol_host_tribe_colorstrip.txt`](../result/annotations/iupac_itol_host_tribe_colorstrip.txt) |
| Consensus status | [`strict_itol_consensus_status_colorstrip.txt`](../result/annotations/strict_itol_consensus_status_colorstrip.txt) | [`iupac_itol_consensus_status_colorstrip.txt`](../result/annotations/iupac_itol_consensus_status_colorstrip.txt) |

## Key Interpretation

The consensus50 reference makes the analysis simpler and biologically cleaner than the earlier multi-reference `nifH` strategy because there is now one `nifH` reference sequence instead of multiple `nifH` target regions. This means each sample contributes at most one nifH consensus sequence to the tree, which is easier to interpret than the previous sample-target design.

However, this is still a reference-guided method. It can recover genes similar enough to the consensus reference to map and pass coverage thresholds, but it can miss highly divergent or novel nif/nod genes. Such genes would usually appear as low or failed coverage, not as an explicit error.

## Important Caveats

1. MC controls are limited to n=2, so they should be interpreted cautiously.
2. BLAN is a real NEON field site and should not be treated as a negative control.
3. Mixed-IUPAC calls are threshold-based, using the exact rules listed above.
4. Mixed signal suggests possible multiple templates, multicopy signal, or strain mixture, but it is not direct proof of multiple organisms.
5. The strict tree tips represent sample-level nifH consensus sequences. The mixed-IUPAC tree tips also represent sample-level nifH consensus sequences, but include ambiguity codes at mixed positions.
6. Reference-guided mapping can miss novel or highly divergent nif/nod genes that do not align well to the reference.
7. BWA-MEM does not biologically interpret IUPAC ambiguity as a set of possible bases during mapping. Ambiguity handling is applied later during pileup-based consensus calling.

## Current Conclusions

1. The consensus50 reference includes 72 functional genes and is a cleaner starting point than the older multi-region reference for gene-tree construction.
2. Across 2,907 samples, 2,052 samples have at least one gene with >=80% coverage and mean depth >=10.
3. `nifH` is strongly recovered and remains the best first gene for tree construction.
4. The nifH strict tree uses 387 conservative single-dominant sequences.
5. The nifH mixed-IUPAC tree uses 1,538 sequences and preserves possible mixed-template signal as ambiguity codes.
6. Metadata matching is mostly complete, but 24 analysis-only and 15 metadata-only sample IDs remain to confirm.

## Next Steps

1. Copy completed IQ-TREE output files into [`result/trees`](../result/trees).
2. Upload the strict and mixed-IUPAC trees to iTOL and apply the annotation files in [`result/annotations`](../result/annotations).
3. Confirm the remaining metadata mismatches with Ryan/Pranoti.
4. Add a short diagnostic comparing coverage at A/C/G/T, IUPAC, and `N` positions in the consensus50 reference, especially for `nifH`.
5. Repeat the same consensus/tree workflow for other high-coverage genes, such as `nifD`, `nifK`, or `nifA`, if desired.

