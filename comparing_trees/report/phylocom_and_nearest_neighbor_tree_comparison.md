

# Phylocom And Nearest-Neighbor Tree Comparison

This section describes two related analyses that use the functional-gene trees to ask whether samples with the same label are close together in the tree.

The labels tested here include the closest BLAST genus assigned to each sample consensus sequence. The most important question is:

> If several samples are assigned to the same closest BLAST genus, are those samples close together in the gene tree?

The analyses were run separately for each gene and each tree set:

- Strict single-dominant trees
- Mixed-IUPAC trees


## Part 1. Phylocom-Style Clustering Test

### Purpose

I used the `phylocomr` package to test whether samples with the same label are more clustered in the tree than expected by random labels.

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

This is a group-level test. It gives one p-value for the whole group, **not one p-value for each sample**.

### Package, Function, Input, And Output

Package:

```r
library(phylocomr)
phylocomr::ph_comstruct()
```
 
 ph_comstruct() takes a phylogenetic tree and a sample/group table, randomizes tip labels, and tests whether tips in each group are closer together than expected by random labels.



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
| `nri` | Net Relatedness Index; **positive values mean same-label tips are closer than random** |
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

In this report, I interpret a group as having supported clustering when:

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

![Gene-by-genus BLAST assignment summary with matched Phylocom NRI](../result/phylocom_clustering/figures/06_blast_genus_by_gene_heatmap_min5_with_phylocom_NRI.png)

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

The Phylocom result is useful, but it is not a per-sample result. Because of that, I added a second analysis to answer a more direct question:

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

For one sample/tip, the nearest neighbor is the other tip with the smallest tree distance (dist_mat <- cophenetic.phylo(tr)).

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

I also asked whether the observed nearest-neighbor percentage is higher than expected by random labels.

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

I kept the tree fixed and randomly shuffled only the closest BLAST genus labels across the tips. Then we compared the real result to the randomized results to test whether samples assigned to the same genus were closer together than expected by chance.
```text
Mesorhizobium-assigned nifH samples are much more often nearest to other
Mesorhizobium-assigned nifH samples than expected by random labels.
```

### Nearest-Neighbor Result Files

Code:

- [12_nearest_neighbor_same_blast_genus_in_trees.R](../code/Rcode/12_nearest_neighbor_same_blast_genus_in_trees.R)

Result:

- [12_nearest_neighbor_by_sample.tsv](../result/tree_nearest_neighbor/tables/12_nearest_neighbor_by_sample.tsv)
- [12_nearest_neighbor_summary_by_gene_genus.tsv](../result/tree_nearest_neighbor/tables/12_nearest_neighbor_summary_by_gene_genus.tsv)
- [12_nearest_neighbor_summary_with_random_test.tsv](../result/tree_nearest_neighbor/tables/12_nearest_neighbor_summary_with_random_test.tsv)


![Nearest-neighbor same BLAST genus heatmap](../result/tree_nearest_neighbor/figures/12_nearest_neighbor_same_blast_genus_heatmap.png)


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


The two analyses answer related but different questions.

| Analysis | Question | Result type |
|---|---|---|
| Phylocom NRI | Are all tips in a label group closer together than random? | One NRI and one p-value per group |
| Nearest-neighbor analysis | For each sample, is the closest tree neighbor the same BLAST genus? | One row per sample, plus a percent summary per group |



I first used phylocomr::ph_comstruct to test whether tips assigned to the same
closest BLAST genus were phylogenetically clustered within each functional-gene
tree. I then performed a sample-level nearest-neighbor analysis to make the
tree pattern easier to interpret. For each sample, I identified the closest
tip in the tree and asked whether that nearest neighbor had the same closest
BLAST genus assignment. This provided an interpretable percentage for each
gene-genus group.


Example result sentence:


For the Mixed-IUPAC nifH tree, 597 samples were assigned to Mesorhizobium by
BLAST. Of these, 454 samples had a nearest tree neighbor also assigned to
Mesorhizobium (76.0%), compared with a random expectation of 38.8%
(FDR = 0.0016). This indicates that Mesorhizobium-assigned nifH samples are
not randomly distributed across the tree, but tend to occur near other
Mesorhizobium-assigned nifH samples.


## Running code

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

## Reference

Webb, C. O., Ackerly, D. D., & Kembel, S. W. (2008). Phylocom: software for the analysis of phylogenetic community structure and trait evolution. *Bioinformatics*, 24(18), 2098-2100. DOI: [10.1093/bioinformatics/btn358](https://doi.org/10.1093/bioinformatics/btn358)



## Literature Support For Closest BLAST Genus Patterns

The closest BLAST genus assignments should be interpreted as the closest available reference sequence, not as confirmed species identity. However, the main genera recovered in the functional-gene trees are consistent with known nitrogen-fixing and nodulating bacterial groups reported in the literature.

Mesorhizobium is a well-known legume symbiont genus. Colombi et al. (2023) showed that Mesorhizobium can become nitrogen-fixing symbionts through horizontal transfer of symbiosis-gene-carrying integrative and conjugative elements (ICESyms), supporting the interpretation of Mesorhizobium-associated nif/nod signals in our results. Colombi et al. (2021) also showed that Mesorhizobium symbiosis genes are often carried on mobile ICE elements.

Bradyrhizobium is another major nitrogen-fixing legume symbiont genus. Ormeño-Orrillo and Martínez-Romero (2019) describe Bradyrhizobium as a diverse genus that includes nitrogen-fixing nodule-forming bacteria. The Genisteae rhizobial review by Stępkowski et al. (2018) also shows that nifD and nifH are commonly used symbiotic-gene markers for Bradyrhizobium and related rhizobia.

Microvirga is also supported as a nitrogen-fixing root-nodule bacterial genus. Ardley et al. (2012) described Microvirga lupini, Microvirga lotononidis, and Microvirga zambiensis as alphaproteobacterial root-nodule bacteria that nodulate and fix nitrogen, and used concatenated nifD and nifH sequences in their phylogenetic analysis.

Sinorhizobium is strongly supported as a symbiotic nitrogen-fixing genus. Barnett et al. (2001) reported the full sequence of the Sinorhizobium meliloti pSymA megaplasmid, which contains many genes involved in symbiosis, including nitrogen-fixation genes.

For nodulation-associated genes, Nod factor biosynthesis literature supports the biological relevance of genes such as nodL and nodX. NodL is involved in Nod factor acetylation, and nodX is a known host-specificity gene in Rhizobium leguminosarum bv. viciae, especially in the Afghanistan pea nodulation system.

Together, these references support our interpretation that repeated clustering of samples assigned to Mesorhizobium, Bradyrhizobium, Microvirga, Sinorhizobium, and Rhizobium in nif/nodulation-associated gene trees is biologically meaningful, while still requiring cautious wording because BLAST genus labels represent closest reference matches rather than confirmed isolate identities.

## References

Colombi, E., Hill, Y., Lines, R., Sullivan, J. T., Kohlmeier, M. G., Christophersen, C. T., Ronson, C. W., Terpolilli, J. J., & Ramsay, J. P. (2023). Population genomics of Australian indigenous Mesorhizobium reveals diverse nonsymbiotic genospecies capable of nitrogen-fixing symbioses following horizontal gene transfer. Microbial Genomics, 9(1), 000918. https://doi.org/10.1099/mgen.0.000918

Colombi, E., Perry, B. J., Sullivan, J. T., Bekuma, A. A., Terpolilli, J. J., Ronson, C. W., & Ramsay, J. P. (2021). Comparative analysis of integrative and conjugative mobile genetic elements in the genus Mesorhizobium. Microbial Genomics, 7(10), 000657. https://doi.org/10.1099/mgen.0.000657

Ormeño-Orrillo, E., & Martínez-Romero, E. (2019). A genomotaxonomy view of the Bradyrhizobium genus. Frontiers in Microbiology, 10, 1334. https://doi.org/10.3389/fmicb.2019.01334

Stępkowski, T., Banasiewicz, J., Granada, C. E., Andrews, M., & Passaglia, L. M. P. (2018). Phylogeny and phylogeography of rhizobial symbionts nodulating legumes of the tribe Genisteae. Genes, 9(3), 163. https://doi.org/10.3390/genes9030163

Ardley, J. K., Parker, M. A., De Meyer, S. E., Trengove, R. D., O'Hara, G. W., Reeve, W. G., Yates, R. J., Dilworth, M. J., Willems, A., & Howieson, J. G. (2012). Microvirga lupini sp. nov., Microvirga lotononidis sp. nov., and Microvirga zambiensis sp. nov. are alphaproteobacterial root nodule bacteria that specifically nodulate and fix nitrogen with geographically and taxonomically separate legume hosts. International Journal of Systematic and Evolutionary Microbiology, 62, 2579-2588. https://doi.org/10.1099/ijs.0.035097-0

Barnett, M. J., Fisher, R. F., Jones, T., et al. (2001). Nucleotide sequence and predicted functions of the entire Sinorhizobium meliloti pSymA megaplasmid. Proceedings of the National Academy of Sciences, 98(17), 9883-9888. https://doi.org/10.1073/pnas.161294798

López-Lara, I. M., Kafetzopoulos, D., Spaink, H. P., & Thomas-Oates, J. E. (2001). Rhizobial NodL O-acetyl transferase and NodS N-methyl transferase functionally interfere in production of modified Nod factors. Journal of Bacteriology, 183(11), 3408-3416. https://doi.org/10.1128/JB.183.11.3408-3416.2001

Hogg, B., Davies, A. E., Wilson, K. E., Bisseling, T., & Downie, J. A. (2002). Competitive nodulation blocking of cv. Afghanistan pea is related to high levels of nodulation factors made by some strains of Rhizobium leguminosarum bv. viciae. Molecular Plant-Microbe Interactions, 15(1), 60-68. https://doi.org/10.1094/MPMI.2002.15.1.60




# Literature Confirmation of Closest-BLAST-Genus Results

Closest-BLAST-genus assignments represent the nearest available reference sequence, not confirmed isolate identity. The table below cross-checks every genus recovered at meaningful frequency (≥5 samples, either tree set) against the published literature on nitrogen fixation and legume/actinorhizal nodulation, using the actual sample counts recovered by this pipeline (`10_blast_count_vs_phylocom_ntaxa_check.tsv`).

## Confirmation table

| Gene | Group | Closest BLAST genus | n (max, strict/mixed) | Confirmed in literature | Key citation |
|---|---|---|---:|---|---|
| `nifH` | nif | Mesorhizobium | 597 | Yes — legume symbiont; acquires nif/nod genes via ICESym horizontal transfer | Colombi et al. (2023) |
| `nifH` | nif | Bradyrhizobium | 416 | Yes — major nitrogen-fixing, nodule-forming legume symbiont genus | Ormeño-Orrillo & Martínez-Romero (2019) |
| `nifH` | nif | Microvirga | 220 | Yes — alphaproteobacterial root-nodule bacteria; nifD/nifH used in original species description | Ardley et al. (2012) |
| `nifH` | nif | Rhizobium | 128 | Yes — archetypal nodulating/nif genus | Davis et al. (1988) |
| `nifH` | nif | Sinorhizobium | 104 | Yes — nif genes physically mapped on the *S. meliloti* pSymA megaplasmid | Barnett et al. (2001) |
| `nifH` | nif | Burkholderia | 31 | Yes — "beta-rhizobia"; nodulate *Mimosa* spp. and fix N2 | Moulin et al. (2001) |
| `nifH` | nif | Neorhizobium | 21 | Yes — nodulates *Galega* spp. | Mousavi et al. (2014) |
| `nifH` | nif | Ensifer | 17 | Taxonomic synonym of *Sinorhizobium* — not independent confirmation | Martens et al. (2008) |
| `nifD` | nif | Microvirga | 564 | Yes — nifD used directly in species description | Ardley et al. (2012) |
| `nifD` | nif | Mesorhizobium | 455 | Yes — nod/nif genes carried on mobile ICE elements | Colombi et al. (2021) |
| `nifD` | nif | Sinorhizobium | 247 | Yes | Barnett et al. (2001) |
| `nifD` | nif | Rhizobium | 82 | Yes | Davis et al. (1988) |
| `nifD` | nif | Azorhizobium | 37 | Yes — stem/root-nodulating N2-fixer of *Sesbania rostrata* | Dreyfus et al. (1988) |
| `nifD` | nif | Ensifer | 31 | Taxonomic synonym of *Sinorhizobium* — not independent confirmation | Martens et al. (2008) |
| `nifD` | nif | Bradyrhizobium | 23 | Yes; nifD/nifH are the standard symbiotic marker genes for this genus | Stępkowski et al. (2018) |
| `nifD` | nif | Methylobacterium | 17 | Partial — only *M. nodulans* and a few pink-pigmented strains confirmed nodulating; most *Methylobacterium* do not | Jourand et al. (2004) |
| `nifD` | nif | Burkholderia | 6 | Yes | Moulin et al. (2001) |
| `nifK` | nif | Mesorhizobium | 382 | Yes (via linked *nifHDK* operon; no genus-specific *nifK* paper identified) | Colombi et al. (2021) |
| `nifK` | nif | Sinorhizobium | 354 | Yes — *nifK* is part of the physically mapped *nifHDK* operon on pSymA | Barnett et al. (2001) |
| `nifK` | nif | Microvirga | 176 | Yes (via linked operon; no genus-specific *nifK* paper identified) | Ardley et al. (2012) |
| `nifK` | nif | Rhizobium | 140 | Yes | Davis et al. (1988) |
| `nifJ` | nif | "Candidatus" (= *Candidatus Frankia alpina*) | 325 | **Reference-database artifact, not independent confirmation.** Only 2 nifJ reference sequences exist in the taxon set, both *Frankia*; every nifJ sample matches by default | Pozzi et al. (2020) |
| `nodL` | canonical nod | Sinorhizobium | 211 | Yes — NodL acetylation of Nod factors demonstrated in *S. meliloti* | López-Lara et al. (2001) |
| `nodL` | canonical nod | Rhizobium | 20 | Yes — NodL originally characterized in *R. leguminosarum* | Bloemberg et al. (1994) |
| `nodX` | canonical nod | Rhizobium | 133 | Yes — nodX is the defining *R. leguminosarum* bv. viciae host-specificity gene (Afghanistan pea system) | Hogg et al. (2002); Davis et al. (1988) |
| `nolG`, `nolF`, `noeA`, `noeB` | accessory / Other | Sinorhizobium | 149–229 | **Reference-database artifact, not independent confirmation.** Taxon reference sets for these four genes are 70–97% Sinorhizobium by composition (52/78, 64/74, 55/59, 57/59 respectively), leaving little alternative genus to assign | Barnett et al. (2001) |


- **Ensifer vs. Sinorhizobium**: these are taxonomic synonyms for the same genus (Martens et al. 2008); treating them as independent genera in counts/statistics will understate real Sinorhizobium-genus support and should be merged before any downstream analysis that counts "distinct genera."
- **nifJ and the four accessory genes** (`nolG`, `nolF`, `noeA`, `noeB`) show apparently strong single-genus signal largely because their taxon reference alignments contain few or no alternative genera — this should be reported as a reference-database limitation, not as biological confirmation of genus identity.
- All other genus assignments are corroborated by dedicated primary literature on nitrogen fixation and/or nodulation for that genus.

## References

Ardley, J. K., Parker, M. A., De Meyer, S. E., Trengove, R. D., O'Hara, G. W., Reeve, W. G., Yates, R. J., Dilworth, M. J., Willems, A., & Howieson, J. G. (2012). *Microvirga lupini* sp. nov., *Microvirga lotononidis* sp. nov., and *Microvirga zambiensis* sp. nov. are alphaproteobacterial root nodule bacteria that specifically nodulate and fix nitrogen with geographically and taxonomically separate legume hosts. *International Journal of Systematic and Evolutionary Microbiology*, 62, 2579–2588. https://doi.org/10.1099/ijs.0.035097-0

Barnett, M. J., Fisher, R. F., Jones, T., Komp, C., Abola, A. P., Barloy-Hubler, F., Bowser, L., Capela, D., Galibert, F., Gouzy, J., Gurjal, M., Hong, A., Huizar, L., Hyman, R. W., Kahn, D., Kahn, M. L., Kalman, S., Keating, D. H., Palm, C., Peck, M. C., Surzycki, R., Wells, D. H., Yeh, K.-C., Davis, R. W., Federspiel, N. A., & Long, S. R. (2001). Nucleotide sequence and predicted functions of the entire *Sinorhizobium meliloti* pSymA megaplasmid. *Proceedings of the National Academy of Sciences*, 98(17), 9883–9888. https://doi.org/10.1073/pnas.161294798

Bloemberg, G. V., Thomas-Oates, J. E., Lugtenberg, B. J. J., & Spaink, H. P. (1994). Nodulation protein NodL of *Rhizobium leguminosarum* O-acetylated lipo-oligosaccharides, chitin fragments and N-acetylglucosamine in vitro. *Molecular Microbiology*, 11(5), 793–804. https://doi.org/10.1111/j.1365-2958.1994.tb00356.x

Colombi, E., Perry, B. J., Sullivan, J. T., Bekuma, A. A., Terpolilli, J. J., Ronson, C. W., & Ramsay, J. P. (2021). Comparative analysis of integrative and conjugative mobile genetic elements in the genus *Mesorhizobium*. *Microbial Genomics*, 7(10), 000657. https://doi.org/10.1099/mgen.0.000657

Colombi, E., Hill, Y., Lines, R., Sullivan, J. T., Kohlmeier, M. G., Christophersen, C. T., Ronson, C. W., Terpolilli, J. J., & Ramsay, J. P. (2023). Population genomics of Australian indigenous *Mesorhizobium* reveals diverse nonsymbiotic genospecies capable of nitrogen-fixing symbioses following horizontal gene transfer. *Microbial Genomics*, 9(1), 000918. https://doi.org/10.1099/mgen.0.000918

Davis, E. O., Evans, I. J., & Johnston, A. W. B. (1988). Identification of nodX, a gene that allows *Rhizobium leguminosarum* biovar viciae strain TOM to nodulate Afghanistan peas. *Molecular and General Genetics*, 212(3), 531–535. https://doi.org/10.1007/BF00330860

Dreyfus, B., Garcia, J.-L., & Gillis, M. (1988). Characterization of *Azorhizobium caulinodans* gen. nov., sp. nov., a stem-nodulating nitrogen-fixing bacterium isolated from *Sesbania rostrata*. *International Journal of Systematic Bacteriology*, 38(1), 89–98. https://doi.org/10.1099/00207713-38-1-89

Hogg, B., Davies, A. E., Wilson, K. E., Bisseling, T., & Downie, J. A. (2002). Competitive nodulation blocking of cv. Afghanistan pea is related to high levels of nodulation factors made by some strains of *Rhizobium leguminosarum* bv. viciae. *Molecular Plant-Microbe Interactions*, 15(1), 60–68. https://doi.org/10.1094/MPMI.2002.15.1.60

Jourand, P., Giraud, E., Béna, G., Sy, A., Willems, A., Gillis, M., Dreyfus, B., & de Lajudie, P. (2004). *Methylobacterium nodulans* sp. nov., for a group of aerobic, facultatively methylotrophic, legume root-nodule-forming and nitrogen-fixing bacteria. *International Journal of Systematic and Evolutionary Microbiology*, 54, 2269–2273. https://doi.org/10.1099/ijs.0.02902-0

López-Lara, I. M., Kafetzopoulos, D., Spaink, H. P., & Thomas-Oates, J. E. (2001). Rhizobial NodL O-acetyl transferase and NodS N-methyl transferase functionally interfere in production of modified Nod factors. *Journal of Bacteriology*, 183(11), 3408–3416. https://doi.org/10.1128/JB.183.11.3408-3416.2001

Martens, M., Dawyndt, P., Coopman, R., Gillis, M., De Vos, P., & Willems, A. (2008). Advantages of multilocus sequence analysis for taxonomic studies: a case study using 10 housekeeping genes in the genus *Ensifer* (including former *Sinorhizobium*). *International Journal of Systematic and Evolutionary Microbiology*, 58, 200–214. https://doi.org/10.1099/ijs.0.65392-0

Moulin, L., Munive, A., Dreyfus, B., & Boivin-Masson, C. (2001). Nodulation of legumes by members of the β-subclass of Proteobacteria. *Nature*, 411, 948–950. https://doi.org/10.1038/35082070

Mousavi, S. A., Österman, J., Wahlberg, N., Nesme, X., Lavire, C., Vial, L., Paulin, L., de Lajudie, P., & Lindström, K. (2014). Phylogeny of the Rhizobium–Allorhizobium–Agrobacterium clade supports the delineation of *Neorhizobium* gen. nov. *Systematic and Applied Microbiology*, 37(3), 208–215. https://doi.org/10.1016/j.syapm.2013.12.007

Ormeño-Orrillo, E., & Martínez-Romero, E. (2019). A genomotaxonomy view of the *Bradyrhizobium* genus. *Frontiers in Microbiology*, 10, 1334. https://doi.org/10.3389/fmicb.2019.01334

Pozzi, A. C. M., Herrera-Belaroussi, A., Schwob, G., Bautista-Guerrero, H. H., Bethencourt, L., Fournier, P., Dubost, A., Abrouk, D., Normand, P., & Fernandez, M. P. (2020). Proposal of 'Candidatus Frankia alpina', the uncultured symbiont of *Alnus alnobetula* and *A. incana* that forms spore-containing nitrogen-fixing root nodules. *International Journal of Systematic and Evolutionary Microbiology*, 70(10), 5453–5459. https://doi.org/10.1099/ijsem.0.004433

Stępkowski, T., Banasiewicz, J., Granada, C. E., Andrews, M., & Passaglia, L. M. P. (2018). Phylogeny and phylogeography of rhizobial symbionts nodulating legumes of the tribe Genisteae. *Genes*, 9(3), 163. https://doi.org/10.3390/genes9030163

