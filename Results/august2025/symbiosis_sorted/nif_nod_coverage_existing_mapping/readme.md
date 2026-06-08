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

