## Simple Pipeline Flowchart With Inputs, Code, and Outputs

Ryan’s notes explain why I used `symbiosis_sorted` for this pipeline:

> “It takes a lot of compute to sort these loci out, so I thought I would send along pre-sorted subsets. These are in 16S_sorted, ITS_sorted, and symbiosis_sorted.”

> “I have attached a genomic reference in case you are interested in the functional genes in symbiosis_sorted.”

> “I am still finalizing functional gene data, which are subject to a different pipeline. Claudia and I already discussed using this to derive phylogenetic trees...”

Based on this, I used `symbiosis_sorted` for the `nif`/`nod` functional-gene analysis.


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

---

### 2. Trim Reads With fastp

**Purpose:** Clean adapters, low-quality bases, poly-G/poly-X artifacts, and short reads. Followed Ryan email.

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
$BASE/symbiosis_trimmed_ryan_full_manifest.tsv
```

**Result:**

```text
1,116 paired samples trimmed
0 failed samples
```

**Output in Git:**

https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/QC%20Check

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
Mapping QC summaries were uploaded here:
```
https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/mapping_QC_checking

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
231 unique central nif/nod target regions extracted
169 nif regions
62 nod regions
6 genes were not found in the GenBank annotation
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
Nodule samples showed strongest support for several nif genes.
Best-supported genes included nifH, nifU, nifD, nifK, nifE, and nifB.
BLAN controls also showed target-associated signal, so blank-aware filtering was needed.
```

**Output in Git:**

https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_coverage_existing_mapping


---

### 11. Rank Targets With Blank-Aware Filtering

**Purpose:** Choose targets with strong nodule signal and lower BLAN-control background.

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
$BASE/nif_nod_coverage_existing_mapping/nif_nod_target_blank_aware_ranking_v2.tsv
```

**Result:**

```text
21 target regions were classified as PROMISING.
With a stricter pilot filter, 10 targets remained.
All 10 strict pilot targets were nif genes.
No nod target passed the strict pilot filter.
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

```text
https://itol.embl.de/
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


---

### Next Step

## What This Pipeline Produced So Far

- [x] Cleaned `symbiosis_sorted` reads using Ryan’s `fastp` settings.

- [x] Trimming QC summaries and figures.

- [x] Mapped reads against Ryan’s `symbiosis_islands.fasta` reference.

- [x] Mapping QC summaries and figures.

- [x] A central `nif`/`nod` target-reference set from `symbiosis_islands.gb`.

- [x] A table showing which `nif`/`nod` genes were found or missing in the GenBank annotation.

- [x] BLAN-control mapping results for QC.

- [x] A matched-coordinate table connecting extracted `nif`/`nod` targets back to Ryan’s original `symbiosis_islands.fasta`.

- [x] A region-level coverage table for all samples and matched `nif`/`nod` target locations.

- [x] Coverage summaries by sample type: BLAN, nodule (`No`), rhizosphere (`Rh`), and root (`Ro`).

- [x] A blank-aware target ranking table.

- [x] A list of `PROMISING` target regions for functional-gene tree building.

- [x] A pilot `nifH` consensus FASTA file.

- [x] A filtered high-quality `nifH` consensus FASTA file.

- [x] A MAFFT alignment of 104 high-quality `nifH` sequences.

- [x] A pilot `nifH` maximum-likelihood tree built with IQ-TREE.

- [x] A report-ready pilot `nifH` tree visualization in iTOL.

## What Still Needs To Be Done

- Build additional trees for the other strong `nif` targets.

- Decide whether to use one best `nif` gene tree or combine several `nif` genes.

- Re-check `nod` targets using a nod-focused strategy.

- Build a `nod` gene tree if suitable `nod` consensus sequences can be recovered.

- Build or obtain the core genome tree.

- Compare the core genome tree, `nif` tree, and `nod` tree.

- Compare bacterial trees with plant phylogeny.

- Compare bacterial trees with plant-microbe interaction networks.

- Use 16S data for bacterial community composition and interaction-network analysis.

- Use ITS data if fungal/community context is needed.

- Combine tree results, network results, and environmental metadata to test fuzzy specificity.

### References:

### References supporting Step 11: Blank-aware target ranking

| Pipeline decision | Paper | Exact sentence from the paper | How it helps this step |
|---|---|---|---|
| Focus first on nodule samples (`No`) | Sprent et al. / review on legume-rhizobia specificity, **“Specificity in Legume-Rhizobia Symbioses”** | “Most legume species can fix atmospheric nitrogen (N2) via symbiotic bacteria (general term ‘rhizobia’) in root nodules...” | This supports using nodule samples as the main biological group for choosing the first `nif`/`nod` phylogenetic target, because nodules are the tissue where nitrogen-fixing symbionts are expected to be strongest. |
| Use BLAN controls before selecting targets | Salter et al. 2014, **“Reagent and laboratory contamination can critically impact sequence-based microbiome analyses”** | “Concurrent sequencing of negative control samples is strongly advised.” | This supports keeping BLAN controls in the QC step instead of immediately removing them, because blanks help identify background or contamination signal. |
| Compare biological samples with negative controls | Davis et al. 2018, **“Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data”** | “Sequences from contaminating taxa are likely to have higher prevalence in control samples than in true samples.” | This supports our blank-aware ranking logic: a target is more trustworthy when it is stronger in biological samples than in BLAN controls. |
| Use negative controls to identify problematic sequence features | decontam documentation / Davis and Callahan | “The prevalence ... of each sequence feature in true positive samples is compared to the prevalence in negative controls to identify contaminants.” | This supports comparing target-region support in real samples versus BLAN controls before selecting targets for tree building. |


### References supporting Step 12: Pilot `nifH` consensus extraction

| Pipeline decision | Paper/tool reference | Exact sentence or relevant statement | How it helps this step |
|---|---|---|---|
| Use `nifH` as a nitrogen-fixation marker | Gaby & Buckley 2012, **“A comprehensive evaluation of PCR primers to amplify the nifH gene of nitrogenase”** | “The nifH gene encodes the iron protein subunit of nitrogenase, the enzyme responsible for biological nitrogen fixation.” | This supports choosing `nifH` as a biologically meaningful pilot gene for nitrogen-fixation phylogeny. |
| Use `nifH` for diversity/phylogenetic analysis of nitrogen-fixing organisms | Gaby & Buckley 2012 | “The nifH gene is the most widely used molecular marker for studying the ecology and evolution of nitrogen-fixing microorganisms.” | This supports building a pilot tree from `nifH` sequences. |
| Generate consensus sequence from mapped reads | BCFtools consensus documentation | “The consensus command applies VCF variants to a reference sequence to create consensus sequence.” | This supports the computational step where mapped reads are converted into one representative sequence per sample. |
| Use SAMtools/BCFtools for variant/consensus workflows | Danecek et al. 2021, **“Twelve years of SAMtools and BCFtools”** | SAMtools and BCFtools are described as tools for “high-throughput sequencing data processing.” | This supports using `samtools`, `bcftools mpileup`, variant calling, and consensus generation from BAM files. |
