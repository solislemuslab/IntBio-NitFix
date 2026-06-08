## `nif_nod_region_coverage_all_samples.tsv`

This table summarizes gene-region coverage for every mapped sample across all matched central `nif`/`nod` target locations. It was generated from the existing BAM files produced by mapping trimmed `symbiosis_sorted` reads to Ryan’s `symbiosis_islands.fasta` reference. Coverage was calculated using `samtools depth`.

The file contains **260,029 lines**: one header line plus **260,028 sample-region rows**. This corresponds to **1,116 samples × 233 matched target locations**.

### Column descriptions

| Column | Meaning |
|---|---|
| `sample` | Sample name, such as `BLAN-1-1-No` or `TALL-13-3-Ro`. |
| `sample_type` | Sample tissue/type parsed from the sample name. For example, `No` = nodule, `Rh` = rhizosphere, `Ro` = root. |
| `blank_control` | Indicates whether the sample is a BLAN control. Values are `yes` or `no`. |
| `gene` | Target gene name, such as `nifA`, `nifH`, or `nodC`. |
| `original_reference` | The original reference sequence in Ryan’s `symbiosis_islands.fasta` where the target region was found. |
| `bed_start` | Start coordinate of the target region on the original reference, using BED-style 0-based coordinates. |
| `bed_end` | End coordinate of the target region on the original reference. |
| `target_id` | Unique target identifier created during the `nif`/`nod` target extraction step. It includes the gene name, accession/reference information, and target number. |
| `strand` | Strand of the annotated gene region: `+` for forward strand or `-` for reverse strand. |
| `target_length` | Length of the target gene region in base pairs. |
| `covered_bases` | Number of bases in the target region that had at least one mapped read covering them. |
| `percent_covered` | Percentage of the target region covered by at least one mapped read. This is calculated as `covered_bases / target_length × 100`. |
| `mean_depth` | Average read depth across the target region. Higher values mean more reads cover the region on average. |
| `max_depth` | Maximum read depth observed at any single base position within the target region. |

### How to interpret the table

Each row represents one sample tested against one `nif`/`nod` target location. For example, if a row has high `percent_covered` and high `mean_depth`, that means the sample has strong read support for that gene region. If `percent_covered` is low or `mean_depth` is near zero, that target is poorly supported in that sample.

This table is important because it determines which genes and samples are suitable for downstream consensus-sequence extraction and phylogenetic tree construction.

The full gene-region coverage table was generated on the cluster but was not uploaded to GitHub **because it is large**. It is stored here:

```bash
/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv
```



# nif/nod Coverage Summary

This folder contains summarized coverage results for central `nif` and `nod` target regions across the August 2025 `symbiosis_sorted` samples.

The full coverage table is large and was not uploaded to GitHub. It is stored on the cluster here:

```text
/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv
```

That full table contains 260,029 lines: one header plus 260,028 sample-region rows, representing 1,116 samples × 233 matched `nif`/`nod` target locations.

## Step Summary

Coverage was calculated from the existing BAM files that were mapped to Ryan's original `symbiosis_islands.fasta` reference. The matched `nif`/`nod` coordinates from `nif_nod_original_reference_regions/` were used to measure coverage for each gene region in each sample.

A target was counted as having good coverage when:

```text
percent_covered >= 80
mean_depth >= 10
```

This means at least 80% of the gene region was covered, with an average read depth of at least 10 reads per base.

## Files in This Folder

| File | Description | Why it is important |
|---|---|---|
| `nif_nod_gene_coverage_by_sample_type_heatmap.svg` | Figure summarizing coverage by gene and sample group. | Main report-ready visualization for this step. |
| `nif_nod_gene_coverage_summary_by_sample_type.tsv` | Gene-level summary by sample group (`BLAN`, `No`, `Rh`, `Ro`). | Shows which genes have good coverage in each sample type. |
| `nif_nod_target_coverage_summary.tsv` | Target-level summary by gene target and sample group. | Useful for selecting specific target references for consensus extraction. |
| `nif_nod_sample_coverage_summary.tsv` | Sample-level summary of how many target regions and genes passed coverage thresholds. | Useful for identifying samples with strong or weak overall coverage. |
| `make_coverage_summary_figure.py` | Local Python script used to generate the SVG heatmap from the gene-level summary table. | Makes the figure reproducible. |

## Figure Explanation

The heatmap summarizes central `nif`/`nod` gene coverage across sample groups:

- `BLAN` = blank controls
- `No` = nodule samples
- `Rh` = rhizosphere samples
- `Ro` = root samples

Each cell shows the percentage of samples in that group with at least one target region passing the good-coverage threshold. Darker color means a larger fraction of samples had good coverage for that gene.

The strongest biological support was observed for several `nif` genes, especially in nodule (`No`) samples. Genes with strong nodule support included `nifH`, `nifU`, `nifD`, `nifK`, `nifE`, and `nifB`. However, BLAN controls also showed target-associated signal, so final target selection requires blank-aware filtering.

## Short Report Text

Coverage summaries showed strong support for several central `nif` genes, especially in nodule (`No`) samples. However, BLAN controls also showed target-associated signal, so coverage alone was not sufficient for final target selection. These results were used to guide the next blank-aware filtering step before consensus-sequence extraction and phylogenetic tree construction.


## Blank-Aware Target Ranking

The file `nif_nod_target_blank_aware_ranking_v2.tsv` ranks each target by comparing nodule (`No`) support with BLAN-control background. Targets were classified as:

- `PROMISING`: stronger nodule signal with lower blank-control background.
- `BIOLOGICAL_SIGNAL_BUT_BLANK_BACKGROUND`: biological signal was present, but BLAN controls also showed substantial signal.
- `LOW_OR_LIMITED_SIGNAL`: limited nodule support.
- `NOT_SUPPORTED`: insufficient support.

The ranking table classified 21 targets as `PROMISING`. Applying a stricter pilot-tree filter requiring at least 75 good nodule samples left 10 promising targets:

| Gene | Target ID | Good nodule samples | Good BLAN samples | Status |
|---|---|---:|---:|---|
| `nifB` | `nifB|NC_008278|NC_008278_-_Symbiosis_Island|ref15` | 76 | 2 | `PROMISING` |
| `nifB` | `nifB|NC_009937|NC_009937_-_Symbiosis_Island_2|ref16` | 108 | 1 | `PROMISING` |
| `nifD` | `nifD|NC_008278|NC_008278_-_Symbiosis_Island|ref28` | 137 | 2 | `PROMISING` |
| `nifD` | `nifD|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island|ref36` | 133 | 1 | `PROMISING` |
| `nifE` | `nifE|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island|ref49` | 82 | 0 | `PROMISING` |
| `nifH` | `nifH|NC_002678|NC_002678_-_Symbiosis_Island_3|ref52` | 83 | 0 | `PROMISING` |
| `nifH` | `nifH|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island|ref63` | 159 | 1 | `PROMISING` |
| `nifK` | `nifK|NC_009937|NC_009937_-_Symbiosis_Island_4|ref71` | 137 | 0 | `PROMISING` |
| `nifT` | `nifT|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island_2|ref109` | 78 | 0 | `PROMISING` |
| `nifX` | `nifX|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island|ref157` | 93 | 1 | `PROMISING` |

All 10 strict-filter targets were `nif` genes. No `nod` target passed this stricter pilot filter. This does not mean `nod` genes are absent or unimportant; rather, under the current strict coverage and blank-background thresholds, `nod` targets did not have enough clean support for the first pilot tree. Because `nod` genes are part of the project goal, they should be evaluated separately using a `nod`-focused target-selection strategy.

## Short Report Text

Coverage summaries showed strong support for several central `nif` genes, especially in nodule (`No`) samples. However, BLAN controls also showed target-associated signal, so coverage alone was not sufficient for final target selection. Blank-aware ranking classified 21 targets as `PROMISING`; with a stricter filter of at least 75 good nodule samples, 10 promising `nif` targets remained for the first pilot consensus/tree analysis. No `nod` target passed this stricter pilot filter, so `nod` genes should be evaluated separately with a `nod`-focused strategy.





