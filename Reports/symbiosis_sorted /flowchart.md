## symbiosis_sorted Analysis

Ryan’s notes explain why I used `symbiosis_sorted` for this pipeline:

> “It takes a lot of compute to sort these loci out, so I thought I would send along pre-sorted subsets. These are in 16S_sorted, ITS_sorted, and symbiosis_sorted.”

> “I have attached a genomic reference in case you are interested in the functional genes in symbiosis_sorted.”

> “I am still finalizing functional gene data, which are subject to a different pipeline. Claudia and I already discussed using this to derive phylogenetic trees...”
>
> symbiosis_sorted contains the functional-gene data.

> Ryan provided symbiosis_islands.fasta as the genomic reference for the functional genes.

> Ryan said the functional-gene data could be used to “derive phylogenetic trees.”

> Based on this, I used `symbiosis_sorted` for the `nif`/`nod` functional-gene analysis.


16S_sorted
  -> bacterial community / taxonomic composition
  -> Which bacteria are present?
  -> Useful for microbiome composition and interaction networks

ITS_sorted
  -> fungal community / taxonomic composition
  -> Which fungi are present?
  -> Useful if the project includes fungal/community context

symbiosis_sorted
  -> functional symbiosis genes
  -> Which nif/nod/symbiosis genes are present?
  -> **Useful for nif/nod phylogenetic trees**

---
**BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"**
```bash
module load SolisLemus-BioPhylo/2026.04.20
module load samtools-1.9
```
### 1. Copy Raw Reads Into Processed-Data

**Purpose:** Put the raw symbiosis reads into the working folder.

**Input:**

```text
raw symbiosis FASTQ files
*_R1.fq.gz
*_R2.fq.gz
```
**Output:**

```text
$BASE/raw_symbiosis_full/*_R1.fq.gz
$BASE/raw_symbiosis_full/*_R2.fq.gz
```

**Result:**

```text
1,116 R1 files
1,116 R2 files
```

**Output in Git:**

```text
Raw FASTQ files were not uploaded to Git because they are large.
```

The raw sequencing reads are about:
150 bp


---

### 2. Trim Reads With fastp

**Purpose:** Clean adapters, low-quality bases, poly-G/poly-X artifacts, and short reads. **Followed Ryan email**.

**Input:**

```text
$BASE/raw_symbiosis_full/*_R1.fq.gz
$BASE/raw_symbiosis_full/*_R2.fq.gz
```

**Code to run:**

```bash
bash "$BASE/Rscripts/trim_symbiosis_full_ryan.sh"
```

**Output:**

```text
$BASE/symbiosis_trimmed_ryan_full/*_P1.fastq.gz
$BASE/symbiosis_trimmed_ryan_full/*_P2.fastq.gz
$BASE/symbiosis_trimmed_ryan_full/*_U1.fastq.gz
$BASE/symbiosis_trimmed_ryan_full/*_U2.fastq.gz
$BASE/fastp_logs_ryan_full/*.fastp.json
```

**Result:**

```text
1,116 paired samples trimmed
0 failed samples
These are the outputs from the fastp trimming step: P1 and P2 are the cleaned paired reads used for mapping, U1 and U2 are cleaned reads whose mate was removed during trimming, the .fastp.json files contain per-sample trimming/QC statistics, and the manifest records whether each sample finished successfully.
```
---

### 3. Check Trimming Quality

**Purpose:** Confirm that reads are clean enough for mapping.

**Input:**

```text
$BASE/fastp_logs_ryan_full/*.fastp.json
```

**Code to run:**

```bash
bash "$BASE/Rscripts/summarize_fastp_qc.sh" \
  | tee "$BASE/trimmed_fastp_QC_checking/fastp_qc_summary_run.log"

bash "$BASE/Rscripts/make_fastp_quality_before_all_samples.sh" \
  | tee "$BASE/trimmed_fastp_QC_checking/fastp_quality_before_all_samples_run.log"

bash "$BASE/Rscripts/make_fastp_quality_profile_figures.sh" \
  | tee "$BASE/trimmed_fastp_QC_checking/fastp_quality_profile_figures_run.log"
```

**Output:**

```text
$BASE/trimmed_fastp_QC_checking/fastp_logs_ryan_full_summary.tsv
$BASE/trimmed_fastp_QC_checking/fastp_qc_summary_run.log
$BASE/trimmed_fastp_QC_checking/fastp_quality_before_trimming_all_samples_lightpurple_mean_blue.svg
$BASE/trimmed_fastp_QC_checking/fastp_quality_after_trimming_all_samples_lightpurple_mean_blue.svg
```

**Result:**

```text
Quality was good.
Mean read retention: 95.84%
Mean Q30 improved from 95.39% to 97.49%.
```

**Output in Git:**

https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/QC%20Check

---

### 4. Map Reads to Ryan’s Symbiosis Reference

**Purpose:** Find where trimmed reads align on Ryan’s symbiosis reference.

**Input:**

```text
$BASE/symbiosis_trimmed_ryan_full/*_P1.fastq.gz
$BASE/symbiosis_trimmed_ryan_full/*_P2.fastq.gz
$BASE/symbiosis_islands.fasta
```

**Ryan’s reference FASTA contains:**
85 reference sequences/accessions

Their lengths vary a lot:

Shortest reference sequence: 111 bp

Longest reference sequence: 34,967 bp


Interpretation:
symbiosis_islands.fasta is not one genome.
It is a collection of 85 reference sequences.
Some are short gene/CDS sequences.
Some are longer symbiosis-island regions.


**Code to run:**

```bash
bash "$BASE/Rscripts/map_symbiosis_full.sh"
```

**Output:**

```text
$BASE/symbiosis_mapped_full/*.bam
$BASE/symbiosis_mapped_full/*.bam.bai
$BASE/symbiosis_mapping_logs_full/*.flagstat.txt
$BASE/symbiosis_mapping_full_manifest.tsv
```

**Result:**

```text
1,116 BAM files created
1,116 BAM index files created
1,116 flagstat reports created
0 failed samples
```

**Output in Git:**

```text
Large BAM files were not uploaded to Git.
The mapping script produced sorted BAM files, BAM index files, per-sample `samtools flagstat` reports, and a mapping manifest summarizing sample-level completion status.
```

---

### 5. Check Mapping Quality

**Purpose:** Confirm that the reads mapped well to `symbiosis_islands.fasta`.

**Input:**

```text
$BASE/symbiosis_mapping_logs_full/*.flagstat.txt
$BASE/symbiosis_mapping_full_manifest.tsv
```

**Code to run:**

```bash
bash "$BASE/Rscripts/summarize_symbiosis_mapping_qc.sh" \
  | tee "$BASE/symbiosis_mapping_QC_checking/symbiosis_mapping_qc_summary_run.log"
```

**Output:**

```text
$BASE/symbiosis_mapping_QC_checking/symbiosis_mapping_completion_summary.txt
$BASE/symbiosis_mapping_QC_checking/symbiosis_mapping_flagstat_summary.tsv
$BASE/symbiosis_mapping_QC_checking/symbiosis_mapping_full_manifest.tsv
$BASE/symbiosis_mapping_QC_checking/symbiosis_mapping_qc_summary_run.log
$BASE/symbiosis_mapping_QC_checking/symbiosis_mapping_qc_figure.svg
```

**Result:**

```text
Mean mapped reads: 93.3%
Mean properly paired reads: 85.0%
Mapping quality was good across the dataset.
```

**Output in Git:**

https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/mapping_QC_checking

---

### 6. Extract nif/nod Target Regions

**Purpose:** Find the exact central `nif` and `nod` gene regions from Ryan’s annotation.
Ryan’s broader project goal includes three bacterial phylogenies: **a core genome tree, a `nif` gene tree, and a `nod` gene tree**. However, this step only extracts central **`nif` and `nod`** target regions from Ryan’s `symbiosis_islands.gb` annotation. The core genome tree is a separate future analysis and is not generated from this `symbiosis_islands.gb` target-extraction step.

**Input:**

```text
$BASE/symbiosis_islands.gb
$BASE/symbiosis_islands_gene_list.xlsx
```

**Code to run:**

```bash
python3 "$BASE/Rscripts/extract_all_nif_nod_targets.py"
```

**Output:**

```text
$BASE/nif_nod_target_reference/all_central_nif_nod_target_records.tsv
$BASE/nif_nod_target_reference/central_nif_nod_gene_summary.tsv
$BASE/nif_nod_target_reference/all_central_nif_nod_targets.fasta
$BASE/nif_nod_target_reference/all_central_nif_targets.fasta
$BASE/nif_nod_target_reference/all_central_nod_targets.fasta
$BASE/nif_nod_target_reference/central_genes_not_found_in_genbank.txt
```

**Result:**

```text
Ryan’s gene list contained 18 central nif genes and 17 central nod genes.

From these:
- 16 of 18 nif genes were found in the GenBank annotation.
- 2 of 18 nif genes were not found: nifM and nifY.
- 13 of 17 nod genes were found in the GenBank annotation.
- 4 of 17 nod genes were not found: nodE, nodF, nodP, and nodT.

The script extracted 231 unique central nif/nod target regions:
- 169 nif target regions
- 62 nod target regions
```

**Output in Git:**

https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_target_reference

---

### 7. Check BLAN Controls

**Purpose:** Test whether blank controls contain target-associated signal.

**Input:**

```text
$BASE/symbiosis_trimmed_ryan_full/BLAN*_P1.fastq.gz
$BASE/symbiosis_trimmed_ryan_full/BLAN*_P2.fastq.gz
$BASE/nif_nod_target_reference/all_central_nif_nod_targets.fasta
```

**Code to run:**

```bash
bash "$BASE/Rscripts/map_all_blank_nif_nod_targets.sh"
```

**Output:**

```text
$BASE/nif_nod_blank_mapping_summary.tsv
$BASE/nif_nod_blank_mapping_logs/
$BASE/nif_nod_blank_mapping/
```

**Result:**

```text
All 12 BLAN controls showed target-associated signal.
Mapped-read percentages ranged from 47.35% to 92.08%.
BLAN samples were kept for QC but excluded from final biological trees.
```

**Output in Git:**


https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_blank_mapping

---

### 8. Match nif/nod Targets Back to Ryan’s Original Reference

**Purpose:** Connect extracted nif/nod gene regions back to the coordinates in `symbiosis_islands.fasta`.

**Input:**

```text
$BASE/symbiosis_islands.fasta
$BASE/nif_nod_target_reference/all_central_nif_nod_targets.fasta
```

**Code to run:**

```bash
python3 "$BASE/Rscripts/match_nif_nod_targets_to_original_reference.py"
```

**Output:**

```text
$BASE/nif_nod_original_reference_regions/nif_nod_matches_in_original_reference.tsv
$BASE/nif_nod_original_reference_regions/nif_nod_matches_in_original_reference.bed
$BASE/nif_nod_original_reference_regions/nif_nod_match_summary_by_gene.tsv
$BASE/nif_nod_original_reference_regions/nif_nod_targets_without_exact_match.tsv
```

**Result:**

```text
All 231 extracted targets matched Ryan’s original FASTA.
233 exact reference locations were found because nifQ and nifW had multiple exact locations.
```

**Output in Git:**

https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_original_reference_regions

---

### 9. Measure Coverage for Each Gene in Each Sample

**Purpose:** Ask which samples have enough reads covering each nif/nod target location.

**Input:**

```text
$BASE/symbiosis_mapped_full/*.bam
$BASE/nif_nod_original_reference_regions/nif_nod_matches_in_original_reference.tsv
```

**Code to run:**

```bash
bash "$BASE/Rscripts/calculate_nif_nod_region_coverage.sh"
```

**Output:**

```text
$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv
```

**Result:**

```text
260,028 sample-by-target rows generated.
This means:
1,116 samples × 233 nif/nod target locations = 260,028 rows.
```

**Output in Git:**

```text
The full file was too large for Git.
Cluster path:
$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv
Summary files and figures are here:
```
https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_coverage_existing_mapping



---

### 10. Summarize Coverage by Sample Type

**Purpose:** Compare coverage across BLAN, No, Rh, and Ro groups.

**Input:**

```text
$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv
```

**Code to run:**

```bash
python3 "$BASE/Rscripts/summarize_nif_nod_coverage.py"
```

**Output:**

```text
$BASE/nif_nod_coverage_existing_mapping/nif_nod_gene_coverage_summary_by_sample_type.tsv
$BASE/nif_nod_coverage_existing_mapping/nif_nod_target_coverage_summary.tsv
$BASE/nif_nod_coverage_existing_mapping/nif_nod_sample_coverage_summary.tsv
$BASE/nif_nod_coverage_existing_mapping/nif_nod_gene_coverage_by_sample_type_heatmap.svg
```

**Result:**

```text
percent_covered >= 80
mean_depth >= 10

Nodule samples showed strongest support for several nif genes.
Best-supported genes included nifH, nifU, nifD, nifK, nifE, and nifB.
BLAN controls also showed target-associated signal, so blank-aware filtering was needed.
```

**Output in Git:**

https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_coverage_existing_mapping


---

### 11. Rank Targets With Blank-Aware Filtering

**Purpose:** Choose targets with strong nodule signal and lower BLAN-control background. **Choose the best target regions for making trees.**
 After Step 10, we knew that some genes had good coverage in real samples, but we also saw that some BLAN controls had signal too. So we needed to be careful. We had many possible nif/nod target regions.
But not all targets are good for tree building. A good target should have:
```text
No = biological nodule signal (Nodules are the biologically most relevant tissue for nitrogen-fixation genes like nifH, nifK, nifD, etc., because that is where the symbiotic nitrogen-fixing bacteria are expected to be strongest.)
BLAN = background/control signal
 
Does this target have strong signal in real nodule samples?
AND
Does this target have low signal in BLAN controls?
enough samples with good coverage
```

**Input:**

```text
$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv
```

**Code to run:**

```bash
python3 "$BASE/Rscripts/rank_nif_nod_targets_blank_aware_v2.py"
```

**Output:**
```text
no_samples-Total number of nodule (No) samples checked for this target. Usually 366.
blank_samples-Total number of BLAN/control samples checked for this target. Usually 12.
no_good_samples-Number of No samples that had good coverage for this target.
blank_good_samples-Number of BLAN samples that had good coverage for this target.
no_good_fraction-Fraction of No samples with good coverage.
blank_good_fraction-Fraction of BLAN samples with good coverage.
```



```text
$BASE/nif_nod_coverage_existing_mapping/nif_nod_target_blank_aware_ranking_v2.tsv
```

### Blank-Aware Target Ranking Result

 A target was classified as `PROMISING` only if it met all of the following criteria:

```text
no_good_samples >= 50
blank_good_fraction <= 0.25
depth_ratio_no_vs_blank >= 2
good_fraction_difference >= 0.10
```
**Result**
After applying the  blank-aware criteria, the 231 matched target locations were classified as:
```text
PROMISING: 19
BIOLOGICAL_SIGNAL_BUT_BLANK_BACKGROUND: 21
LOW_OR_LIMITED_SIGNAL: 118
NOT_SUPPORTED: 73
```
```text
A target was classified as BIOLOGICAL_SIGNAL_BUT_BLANK_BACKGROUND if:
no_good_samples >= 50
depth_ratio_no_vs_blank >= 1
but it did not pass all PROMISING requirements.
Meaning: the target had biological signal in nodule samples, but BLAN background was too high or the nodule-vs-BLAN difference was not strong enough.
LOW_OR_LIMITED_SIGNAL
A target was classified as LOW_OR_LIMITED_SIGNAL if:
no_good_samples > 0
but it did not pass the stronger rules above.
Meaning: at least one nodule sample had good coverage, but the target did not have enough support for tree building.
NOT_SUPPORTED
A target was classified as NOT_SUPPORTED if:
no_good_samples = 0
Meaning: no nodule samples had good coverage for that target.
```
**Output in Git:**


https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_coverage_existing_mapping


---

### 12. Extract Pilot nifH Consensus Sequences

**Purpose:** Build one representative `nifH` sequence per good nodule sample.

**Input:**

```text
$BASE/symbiosis_islands.fasta
$BASE/symbiosis_mapped_full/*.bam
$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv
$BASE/nif_nod_original_reference_regions/nif_nod_matches_in_original_reference.tsv
```

**Code to run:**

```bash
bash "$BASE/Rscripts/extract_pilot_nifH_ref63_consensus_v2.sh"
```

**Output:**

```text
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_good_nodule_samples.txt
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_all_samples.fasta
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_qc.tsv
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_Nle5.fasta
```

**Result:**

```text
159 full-length 894 bp consensus sequences were recovered.
104 high-quality sequences remained after filtering for N_percent <= 5.
```

**Output in Git:**

https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nifH_ref63_pilot_consensus_v2

---

### 13. Align nifH Sequences With MAFFT

**Purpose:** Line up the 104 high-quality `nifH` sequences base-by-base.
MAFFT is a multiple-sequence-alignment tool. We used it to align the 104 filtered nifH consensus sequences so IQ-TREE could compare the same base positions across samples and build a phylogenetic tree.

**Input:**

```text
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_Nle5.fasta
```

**Code to run:**

```bash
mafft --auto "$FILTERED" > "$OUTDIR/nifH_ref63_consensus_Nle5.mafft.fasta"
```

**Output:**

```text
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_Nle5.mafft.fasta
```

**Result:**

```text
104 sequences aligned.
Alignment length: 894 bp.
```

**Output in Git:**


https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nifH_ref63_pilot_consensus_v2

---

### 14. Build Pilot nifH Tree With IQ-TREE

**Purpose:** Build a pilot maximum-likelihood `nifH` tree.

**Input:**

```text
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_Nle5.mafft.fasta
```

**Code to run:**

```bash
iqtree -s "$ALIGN" -m MFP -B 1000 -T AUTO \
  --prefix "$OUTDIR/nifH_ref63_Nle5_iqtree"
```

**Output:**

```text
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_Nle5_iqtree.treefile
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_Nle5_iqtree.iqtree
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_Nle5_iqtree.log
```

**Result:**

```text
Pilot maximum-likelihood nifH tree created.
IQ-TREE selected GTR+F+I+G4 as the best-fit model for this alignment.
```

**Output in Git:**

https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nifH_ref63_pilot_consensus_v2


---

### 15. Visualize Tree in iTOL

**Purpose:** Make a report-ready tree figure.

**Input:**

```text
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_Nle5_iqtree.treefile
```

**Code/tool to run:**
https://itol.embl.de
```text
Upload a new tree -> choose nifH_ref63_Nle5_iqtree.treefile
```

**Output:**

```text
nifH_ref63_iTOL_pilot_tree.svg
or
nifH_ref63_iTOL_pilot_tree.png
```

**Result:**

```text
Pilot nifH tree was visualized in iTOL.
```

**Output in Git:**

WILL Add the exported iTOL figure to:

https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nifH_ref63_pilot_consensus_v2



**Conclusion**
The pilot nifH tree shows that Ryan’s symbiosis_sorted functional-gene data can successfully produce a sample-level nitrogen-fixation gene phylogeny. From the nodule samples, we recovered 104 high-quality nifH consensus sequences, aligned them across 894 bp, and built a maximum-likelihood tree. The samples show visible phylogenetic structure, meaning their nifH sequences are not all identical and contain enough variation for downstream biological interpretation. At this stage, the tree supports the use of nifH for functional-gene phylogenetic analysis, but the biological meaning of the clusters still needs metadata, such as host plant, site, or sample code.



**Next Steps**

Add metadata colors to the nifH tree, especially host plant/site/sample-code groups, to see whether clusters match biological patterns.

Check bootstrap support values to identify which branches are well supported.

Repeat consensus extraction and tree building for other strong targets such as nifK, nifD, and nifB.

Compare the separate gene trees to see whether different nitrogen-fixation genes show similar sample relationships.

If multiple genes give consistent results, consider building a concatenated multi-gene tree using shared high-quality samples.

---

