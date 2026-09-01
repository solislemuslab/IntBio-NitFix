# Phylocom And Nearest-Neighbor Tree Comparison

This section describes two related analyses that use the functional-gene trees to ask whether samples with the same label are close together in the tree.

The labels tested here include the closest BLAST genus assigned to each sample consensus sequence. The most important question is:

> If several samples are assigned to the same closest BLAST genus, are those samples close together in the gene tree?

The analyses were run separately for each gene and each tree set:

- Strict single-dominant trees
- Mixed-IUPAC trees

Genes included:

`nifH`, `nifD`, `nifK`, `nifJ`, `nodL`, `nolG`, `nolF`, `noeA`, `noeB`, `nodX`

## Part 1. Phylocom-Style Clustering Test

### Purpose

We used the `phylocomr` package to test whether samples with the same label are more clustered in the tree than expected by random labels.

For example, for:

```text
gene = nifH
tree set = Mixed-IUPAC
closest BLAST genus = Mesorhizobium
n = 597 samples
```

the Phylocom-style test asks:

```text
Are the 597 Mesorhizobium-assigned nifH tips closer together in the nifH tree
than expected if genus labels were randomly placed on the same tree?
```

This is a group-level test. It gives one p-value for the whole group, not one p-value for each sample.

### Package, Function, Input, And Output

Package:

```r
library(phylocomr)
```

Function used:

```r
phylocomr::ph_comstruct()
```

One-line comment for how the package runs:

```r
# ph_comstruct() takes a phylogenetic tree and a sample/group table, randomizes tip labels, and tests whether tips in each group are closer together than expected by random labels.
```

Main input to `ph_comstruct()`:

| Input | Meaning in this analysis |
|---|---|
| `phylo` | One functional-gene tree, for example the Mixed-IUPAC `nifH` tree |
| `sample` | A table saying which tree tips belong to each group, for example which tips are `Mesorhizobium` |
| `null_model = 0` | Randomly reshuffle labels/tips to create the null expectation |
| `randomizations = 999` | Run 999 random label tests |
| `abundance = FALSE` | Treat each tip as present once; do not use abundance weights |

Main output from `ph_comstruct()`:

| Output | Meaning |
|---|---|
| `ntaxa` | Number of tips tested in that group |
| `mpd` | Observed mean pairwise phylogenetic distance among tips in the group |
| `mpd_random` | Mean MPD from randomized groups |
| `nri` | Net Relatedness Index; positive values mean same-label tips are closer than random |
| `p_mpd_cluster` | p-value for clustering based on MPD |
| `fdr_mpd` | BH/FDR-corrected p-value |

### NRI Meaning

NRI means Net Relatedness Index.

It is calculated from MPD, the mean pairwise tree distance among samples in a group:

```text
NRI = -1 * (observed MPD - mean random MPD) / sd random MPD
```

Interpretation:

| NRI value | Meaning |
|---|---|
| Positive NRI | Same-label tips are closer together than random |
| NRI near 0 | Same-label tips are about as close as random |
| Negative NRI | Same-label tips are more spread out than random |

In this report, we interpret a group as having supported clustering when:

```text
NRI > 0 and FDR(MPD) < 0.05
```

### Important Sample-Set Check

The Phylocom test must use the same sample set as the BLAST genus heatmap. To make this correct, the updated Step 10 starts from the same deduplicated BLAST assignment table used for the BLAST heatmap:

[11b_blast_assignment_rows_deduplicated_strict_iupac.tsv](../result/comparative_tree_analysis/tables/11b_blast_assignment_rows_deduplicated_strict_iupac.tsv)

The new code checks that:

```text
BLAST heatmap n = Phylocom ntaxa
```

Check table:

[10_blast_count_vs_phylocom_ntaxa_check.tsv](../result/phylocom_clustering/tables/10_blast_count_vs_phylocom_ntaxa_check.tsv)

Result:

```text
All matched BLAST genus counts agree with Phylocom ntaxa.
```

For example:

| Tree set | Gene | Closest BLAST genus | BLAST n | Phylocom ntaxa | NRI | FDR(MPD) |
|---|---|---:|---:|---:|---:|---:|
| Mixed-IUPAC | nifH | Mesorhizobium | 597 | 597 | 5.5921 | 0.0021 |

This means the Phylocom test for Mixed-IUPAC `nifH` + `Mesorhizobium` used the same 597 samples shown in the BLAST heatmap.

### Phylocom Result Files

Code:

- [10_phylocom_clustering_matched_to_blast_heatmap.R](../code/Rcode/10_phylocom_clustering_matched_to_blast_heatmap.R)
- [11_plot_blast_heatmap_with_matched_phylocom_NRI.R](../code/Rcode/11_plot_blast_heatmap_with_matched_phylocom_NRI.R)

Main tables:

- [10_phylocom_report_table_matched.tsv](../result/phylocom_clustering/tables/10_phylocom_report_table_matched.tsv)
- [10_blast_closest_genus_summary_matched_min5.tsv](../result/phylocom_clustering/tables/10_blast_closest_genus_summary_matched_min5.tsv)
- [10_blast_count_vs_phylocom_ntaxa_check.tsv](../result/phylocom_clustering/tables/10_blast_count_vs_phylocom_ntaxa_check.tsv)
- [11_blast_genus_counts_with_matched_phylocom_NRI.tsv](../result/phylocom_clustering/tables/11_blast_genus_counts_with_matched_phylocom_NRI.tsv)
- [11_matched_heatmap_join_check.tsv](../result/phylocom_clustering/tables/11_matched_heatmap_join_check.tsv)

Main figure:

![Gene-by-genus BLAST assignment summary with matched Phylocom NRI](../result/phylocom_clustering/figures/11_blast_genus_by_gene_heatmap_min5_with_matched_phylocom_NRI.png)

Figure file:

[11_blast_genus_by_gene_heatmap_min5_with_matched_phylocom_NRI.png](../result/phylocom_clustering/figures/11_blast_genus_by_gene_heatmap_min5_with_matched_phylocom_NRI.png)

### Example Phylocom Results

Some strong positive NRI results were:

| Tree set | Gene | Closest BLAST genus | n | NRI | FDR(MPD) | Meaning |
|---|---|---:|---:|---:|---:|---|
| Mixed-IUPAC | nifK | Mesorhizobium | 382 | 21.0949 | 0.0021 | Strong clustering |
| Mixed-IUPAC | nifK | Sinorhizobium | 354 | 16.6006 | 0.0021 | Strong clustering |
| Mixed-IUPAC | nifH | Bradyrhizobium | 416 | 16.5468 | 0.0021 | Strong clustering |
| Mixed-IUPAC | nifK | Microvirga | 176 | 14.1392 | 0.0021 | Strong clustering |
| Mixed-IUPAC | nodL | Sinorhizobium | 211 | 10.5318 | 0.0021 | Strong clustering |
| Mixed-IUPAC | nifH | Mesorhizobium | 597 | 5.5921 | 0.0021 | Same-genus tips are closer than random |

## Part 2. Sample-Level Nearest-Neighbor Analysis

### Purpose

The Phylocom result is useful, but it is not a per-sample result. Because of that, we added a second analysis to answer a more direct question:

```text
For each sample, is its nearest neighbor in the tree assigned to the same closest BLAST genus?
```

This produces one row per sample/tip.

For each sample, the code finds:

1. The closest tip to that sample in the tree.
2. The closest BLAST genus of the sample.
3. The closest BLAST genus of the nearest-neighbor tip.
4. Whether the two genera are the same.

### What Is A Nearest Neighbor In The Tree?

Each tree has branch lengths. The distance between two tips is the total branch length connecting them.

For one sample/tip, the nearest neighbor is the other tip with the smallest tree distance.

Example:

| Sample | Distance to sample B | Distance to sample C | Distance to sample D | Nearest neighbor |
|---|---:|---:|---:|---|
| sample A | 0.01 | 0.20 | 0.35 | sample B |

If sample A and sample B are both assigned to `Mesorhizobium`, then sample A is counted as:

```text
nearest_neighbor_same_blast_genus = TRUE
```

If sample A is `Mesorhizobium` but sample B is `Bradyrhizobium`, then sample A is counted as:

```text
nearest_neighbor_same_blast_genus = FALSE
```

### Why This Analysis Is Easier To Explain

This analysis gives a direct fraction.

For example:

```text
454 / 597 nifH Mesorhizobium samples have a nearest neighbor also assigned to Mesorhizobium.
```

That is:

```text
76.0%
```

So the interpretation is:

```text
Most Mixed-IUPAC nifH samples assigned to Mesorhizobium are nearest to another
Mesorhizobium-assigned nifH sample in the tree.
```

### Random Test For The Nearest-Neighbor Result

We also asked whether the observed nearest-neighbor percentage is higher than expected by random labels.

For each gene and tree set, the code:

1. Keeps the same tree.
2. Keeps the same number of tips.
3. Keeps the same number of samples assigned to each genus.
4. Randomly shuffles the genus labels across tips.
5. Recalculates how many samples have a same-genus nearest neighbor.
6. Repeats this 999 times.

This gives a random expectation.

For example, for Mixed-IUPAC `nifH` + `Mesorhizobium`:

| Value | Result |
|---|---:|
| Total samples assigned to Mesorhizobium | 597 |
| Samples whose nearest neighbor is also Mesorhizobium | 454 |
| Observed percent | 76.0% |
| Random expected percent | 38.8% |
| FDR | 0.0016 |

Interpretation:

```text
Mesorhizobium-assigned nifH samples are much more often nearest to other
Mesorhizobium-assigned nifH samples than expected by random labels.
```

### Nearest-Neighbor Result Files

Code:

- [12_nearest_neighbor_same_blast_genus_in_trees.R](../code/Rcode/12_nearest_neighbor_same_blast_genus_in_trees.R)

Main tables:

- [12_nearest_neighbor_by_sample.tsv](../result/tree_nearest_neighbor/tables/12_nearest_neighbor_by_sample.tsv)
- [12_nearest_neighbor_summary_by_gene_genus.tsv](../result/tree_nearest_neighbor/tables/12_nearest_neighbor_summary_by_gene_genus.tsv)
- [12_nearest_neighbor_summary_with_random_test.tsv](../result/tree_nearest_neighbor/tables/12_nearest_neighbor_summary_with_random_test.tsv)
- [12_nearest_neighbor_n_matches_blast_heatmap_check.tsv](../result/tree_nearest_neighbor/tables/12_nearest_neighbor_n_matches_blast_heatmap_check.tsv)

Main figure:

![Nearest-neighbor same BLAST genus heatmap](../result/tree_nearest_neighbor/figures/12_nearest_neighbor_same_blast_genus_heatmap.png)

Figure file:

[12_nearest_neighbor_same_blast_genus_heatmap.png](../result/tree_nearest_neighbor/figures/12_nearest_neighbor_same_blast_genus_heatmap.png)

### Meaning Of Columns In The Sample-Level Table

File:

[12_nearest_neighbor_by_sample.tsv](../result/tree_nearest_neighbor/tables/12_nearest_neighbor_by_sample.tsv)

| Column | Meaning |
|---|---|
| `tree_set` | `strict` or `iupac` |
| `tree_set_label` | Clear tree-set name |
| `gene` | Functional gene |
| `sample_id` | Sample name |
| `tip_label` | Full tip label in the tree |
| `closest_genus_report` | Closest BLAST genus assigned to that sample consensus sequence |
| `nearest_neighbor_tip` | Closest tip to this sample in the tree |
| `nearest_neighbor_sample_id` | Sample name of the nearest-neighbor tip |
| `nearest_neighbor_closest_genus` | Closest BLAST genus assigned to the nearest-neighbor tip |
| `nearest_neighbor_distance` | Tree distance from this sample to its nearest neighbor |
| `nearest_neighbor_same_blast_genus` | `TRUE` if sample and nearest neighbor have the same closest BLAST genus |
| `nearest_neighbor_tie_count` | Number of equally closest tips if there is a tie |
| `sample_type`, `site`, `state`, `host_genus`, `host_tribe`, `native_status` | Metadata for the sample |

### Meaning Of Columns In The Summary Table

File:

[12_nearest_neighbor_summary_with_random_test.tsv](../result/tree_nearest_neighbor/tables/12_nearest_neighbor_summary_with_random_test.tsv)

| Column | Meaning |
|---|---|
| `tree_set` | `strict` or `iupac` |
| `gene` | Functional gene |
| `closest_genus_report` | Closest BLAST genus |
| `n` | Number of samples/tips in that gene-genus group |
| `same_nearest_neighbor_n` | Number of samples whose nearest neighbor has the same closest BLAST genus |
| `same_nearest_neighbor_percent` | Percent of samples whose nearest neighbor has the same closest BLAST genus |
| `random_mean_same_n` | Mean same-genus nearest-neighbor count from randomized labels |
| `random_mean_same_percent` | Mean same-genus nearest-neighbor percent from randomized labels |
| `nearest_neighbor_z` | How far the observed value is from the random mean, in random standard deviations |
| `p_more_same_than_random` | p-value for having more same-genus nearest neighbors than random |
| `fdr_more_same_than_random` | FDR-corrected p-value |
| `interpretation` | Simple interpretation of the result |

### Example Nearest-Neighbor Results

Some strong sample-level nearest-neighbor results were:

| Tree set | Gene | Closest BLAST genus | n | Same-genus nearest neighbors | Observed % | Random expected % | FDR | Meaning |
|---|---|---:|---:|---:|---:|---:|---:|---|
| Mixed-IUPAC | nodL | Sinorhizobium | 211 | 207 | 98.1% | 91.7% | 0.0016 | More same-genus nearest neighbors than random |
| Mixed-IUPAC | nifK | Rhizobium | 140 | 120 | 85.7% | 13.2% | 0.0016 | Strong sample-level clustering |
| Mixed-IUPAC | nifD | Sinorhizobium | 247 | 211 | 85.4% | 16.8% | 0.0016 | Strong sample-level clustering |
| Mixed-IUPAC | nifK | Sinorhizobium | 354 | 302 | 85.3% | 33.7% | 0.0016 | Strong sample-level clustering |
| Mixed-IUPAC | nifH | Bradyrhizobium | 416 | 331 | 79.6% | 26.9% | 0.0016 | Strong sample-level clustering |
| Mixed-IUPAC | nifK | Mesorhizobium | 382 | 302 | 79.1% | 36.2% | 0.0016 | Strong sample-level clustering |
| Mixed-IUPAC | nifH | Mesorhizobium | 597 | 454 | 76.0% | 38.8% | 0.0016 | Strong sample-level clustering |
| Mixed-IUPAC | nifD | Mesorhizobium | 455 | 334 | 73.4% | 31.1% | 0.0016 | Strong sample-level clustering |
| Mixed-IUPAC | nifD | Microvirga | 564 | 390 | 69.1% | 38.6% | 0.0016 | Strong sample-level clustering |

## How To Report These Two Analyses Together

The two analyses answer related but different questions.

| Analysis | Question | Result type |
|---|---|---|
| Phylocom NRI | Are all tips in a label group closer together than random? | One NRI and one p-value per group |
| Nearest-neighbor analysis | For each sample, is the closest tree neighbor the same BLAST genus? | One row per sample, plus a percent summary per group |

Suggested report sentence:

```text
We first used phylocomr::ph_comstruct to test whether tips assigned to the same
closest BLAST genus were phylogenetically clustered within each functional-gene
tree. We then performed a sample-level nearest-neighbor analysis to make the
tree pattern easier to interpret. For each sample, we identified the closest
tip in the tree and asked whether that nearest neighbor had the same closest
BLAST genus assignment. This provided an interpretable percentage for each
gene-genus group.
```

Example result sentence:

```text
For the Mixed-IUPAC nifH tree, 597 samples were assigned to Mesorhizobium by
BLAST. Of these, 454 samples had a nearest tree neighbor also assigned to
Mesorhizobium (76.0%), compared with a random expectation of 38.8%
(FDR = 0.0016). This indicates that Mesorhizobium-assigned nifH samples are
not randomly distributed across the tree, but tend to occur near other
Mesorhizobium-assigned nifH samples.
```

## Reproducibility

Run Step 10 first:

```bash
Rscript comparing_trees/code/Rcode/10_phylocom_clustering_matched_to_blast_heatmap.R
```

Then run Step 11:

```bash
Rscript comparing_trees/code/Rcode/11_plot_blast_heatmap_with_matched_phylocom_NRI.R
```

Then run Step 12:

```bash
Rscript comparing_trees/code/Rcode/12_nearest_neighbor_same_blast_genus_in_trees.R
```

The result files are under:

[comparing_trees/result](../result)

On GitHub, the result folder will be available here after pushing:

https://github.com/solislemuslab/IntBio-NitFix/tree/main/comparing_trees/result

## Reference

Webb, C. O., Ackerly, D. D., & Kembel, S. W. (2008). Phylocom: software for the analysis of phylogenetic community structure and trait evolution. *Bioinformatics*, 24(18), 2098-2100. DOI: [10.1093/bioinformatics/btn358](https://doi.org/10.1093/bioinformatics/btn358)
