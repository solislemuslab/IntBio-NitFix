
| Figure file | Definition |
|---|---|
| `fastp_qc_summary.svg` | Summary of fastp read-trimming quality across all 2,907 V2 samples. Shows read retention, base retention, Q30 before trimming, Q30 after trimming, and Q20 after trimming. |
| `fastp_quality_after_trimming_all_samples_lightpurple_mean_blue.svg` | Per-base read-quality profile after trimming. Light purple lines are individual samples; the blue line is the mean across all samples. |
| `fastp_quality_before_trimming_all_samples_lightpurple_mean_blue.svg` | Per-base read-quality profile before trimming. Light purple lines are individual samples; the blue line is the mean across all samples. |
| `mixed-IUPAC_label.svg` | Mixed-IUPAC nifH tree with tip labels visible. Contains 4,212 sample-target nifH consensus sequences; useful for zoomed inspection but crowded for full-tree display. |
| `mixed-IUPAC_ref.svg` | Mixed-IUPAC nifH tree colored by nifH target reference (`ref52`, `ref56`, `ref60`, `ref62`, `ref63`). Shows how mixed-aware sequences cluster by target reference. |
| `mixed-IUPAC_sample_type.svg` | Mixed-IUPAC nifH tree colored by sample type (`No`, `Rh`, `Ro`). Shows whether nodule, rhizosphere, and root sequences cluster separately or are mixed across nifH lineages. |
| `single_dominant_ref.svg` | Strict single-dominant nifH maximum-likelihood tree colored by nifH target reference. Contains 1,802 high-confidence sample-target consensus sequences. |
| `single_dominant_sample_type.svg` | Strict single-dominant nifH maximum-likelihood tree colored by sample type (`No`, `Rh`, `Ro`). This is the conservative primary tree view. |
| `symbiosis_mapping_qc_histograms.svg` | Mapping QC histogram across all 2,907 samples after mapping trimmed reads to `symbiosis_islands.fasta`. Shows mapped-read percentage and properly paired-read percentage. |
| `v2_gene_coverage_heatmap_pct80_depth10.svg` | Gene-coverage heatmap using the good-coverage threshold: percent covered ≥80% and mean depth ≥10x. Summarizes which nif/nod genes show strong coverage across sample groups. |
| `v2_threshold_sensitivity_summary.svg` | Threshold sensitivity summary showing how retained samples/targets change under different coverage and depth filters. Used to evaluate whether conclusions depend strongly on the chosen thresholds. |
