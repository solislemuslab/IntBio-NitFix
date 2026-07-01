## symbiosis_sorted Analysis V2

**Version 2** 

Ryan clarified that `BLAN` is a real NEON field site https://www.neonscience.org/field-sites/blan , not a blank/negative control. The negative controls should be `MC-1` and `MC-2`. Ryan also raised two important next-analysis concerns: the current pilot tree produces **one dominant consensus sequence per sample**, but real samples may contain multiple nitrogen-fixing organisms; and reference-guided mapping/consensus may miss divergent nif/nod sequences. Version 2 keeps the successful read-cleaning and mapping foundation, then revises the downstream plan to support negative-control checking, broader sample inclusion, and future multi-copy-aware tree construction.


**Sample**

Sample IDs contain both site and sample-type information. For example, in BLAN-3-1-No, BLAN is the NEON field-site code, 3-1 identifies the individual sample or replicate, and No is the sample type. Therefore, BLAN-3-1-No is not a blank control; it is a real nodule sample collected from the BLAN site . In this dataset, samples ending in -No are nodule samples, samples ending in -Rh are rhizosphere samples, and samples ending in -Ro are root samples. Negative controls should instead be identified using Ryan’s MC-1 and MC-2 labels.

sample_type: No, Rh, Ro
```
370 No
376 Rh
370 Ro
```
site+sample type (number of samples)
```
site	sample_type	count
BLAN	No	4
BLAN	Rh	4
BLAN	Ro	4
CLBJ	No	6
CLBJ	Rh	6
```

number of samples per site:
    ```
     12 BLAN
     18 CLBJ
     15 CPER
     24 DCFS
     79 DEJU
     66 DELA
     27 DSNY
      9 GRSM
     26 HEAL
    156 JERC
      1 KONA
     86 LAJA
      4 LENO
      1 MLBS
      6 NIWO
     21 NOGP
     12 OAES
      9 ORNL
     72 OSBS
     15 RMNP
     38 SCBI
     20 SERC
     15 SJER
    222 TALL
      8 TEAK
    102 TOOL
     34 WOOD
     18 YELL
```


**BASE**

```bash
BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"
module load SolisLemus-BioPhylo/2026.04.20
module load samtools-1.9
```

---

### 1. Copy Raw Reads Into Processed-Data

**Purpose:** Put the raw `symbiosis_sorted` FASTQ files into the working folder.

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

**Interpretation:** These are the original sequencing reads for the symbiosis functional-gene dataset. Each read is approximately 150 bp before downstream mapping and consensus construction.

---

### 2. Trim Reads With fastp

**Purpose:** Clean the raw paired-end reads before mapping.

**Input:**

```text
$BASE/raw_symbiosis_full/*_R1.fq.gz
$BASE/raw_symbiosis_full/*_R2.fq.gz
```

**Method:** `fastp`

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
$BASE/fastp_logs_ryan_full/*.fastp.html
```

**Result:**

```text
1,116 paired samples trimmed
0 failed samples
```

**Interpretation:** `P1` and `P2` are the cleaned paired reads used for mapping. `U1` and `U2` are reads whose mate was removed during trimming. The `fastp` JSON/HTML files store per-sample trimming and quality statistics.

---

### 3. Check Trimming Quality

**Purpose:** Confirm that trimming retained most reads and improved read quality.

**Input:**

```text
$BASE/fastp_logs_ryan_full/*.fastp.json
```

**Method:** custom QC scripts summarizing `fastp` reports

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
$BASE/trimmed_fastp_QC_checking/fastp_quality_before_trimming_all_samples_lightpurple_mean_blue.svg
$BASE/trimmed_fastp_QC_checking/fastp_quality_after_trimming_all_samples_lightpurple_mean_blue.svg
$BASE/trimmed_fastp_QC_checking/fastp_quality_all_samples_gray_mean_blue.svg
```

**Result:**

```text
Mean read retention: 95.84%
Mean Q30 improved from 95.39% to 97.49%
```

**Interpretation:** Trimming removed lower-quality sequence while retaining most of the data. This supports using the cleaned reads for reference mapping.

---

### 4. Map Reads to  Symbiosis Reference provided by Ryan

**Purpose:** Align cleaned reads from each sample to symbiosis functional-gene reference provided by Ryan.

**Input:**

```text
$BASE/symbiosis_trimmed_ryan_full/*_P1.fastq.gz
$BASE/symbiosis_trimmed_ryan_full/*_P2.fastq.gz
$BASE/symbiosis_islands.fasta
```

**Method:** `BWA-MEM` for mapping, `samtools` for sorting/indexing

**Ryan's reference FASTA contains:**

```text
85 reference sequences/accessions
Minimum length: 111 bp
Maximum length: 34,967 bp
Average length: ~3,957 bp
```

**Important interpretation:** `symbiosis_islands.fasta` is not a whole genome. It is a collection of 85 selected symbiosis/gene reference sequences. Some entries are short CDS/gene sequences, while others are longer symbiosis-island regions.

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

**Interpretation:** Each BAM file stores the read-level alignment information for one sample: which reads mapped, which reference section they mapped to, and where they mapped on that reference.

---

### 5. Check Mapping Quality

**Purpose:** Confirm that reads mapped well to `symbiosis_islands.fasta`.

**Input:**

```text
$BASE/symbiosis_mapping_logs_full/*.flagstat.txt
$BASE/symbiosis_mapping_full_manifest.tsv
```

**Method:** `samtools flagstat` summaries and custom mapping-QC script

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

**Interpretation:** The mapping results support using the BAM files for downstream nif/nod target coverage analysis. Version 2 keeps this mapping foundation but revises the downstream control interpretation and consensus/tree strategy based on Ryan's feedback.

---

### 6. Extract Central nif/nod Target Regions From GenBank

**Purpose:** Define the nif/nod gene target regions available in GenBank annotation provided by Ryan. These target regions are the reference features that can later be evaluated for coverage, consensus construction, and tree building.

**Input:**

```text
$BASE/symbiosis_islands.gb
$BASE/symbiosis_islands_gene_list.xlsx
```

**Method:** custom GenBank target extraction with `extract_all_nif_nod_targets.py`

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
Ryan's gene list contained:
- 18 central nif genes
- 17 central nod genes

Found in the GenBank annotation:
- 16 of 18 nif genes
- 13 of 17 nod genes

Not found in the GenBank annotation:
- nifM
- nifY
- nodE
- nodF
- nodP
- nodT

Extracted target regions:
- 231 total unique central nif/nod target regions
- 169 nif target regions
- 62 nod target regions
```

**Interpretation:** This step defines which nif/nod genes are available in Ryan's annotated reference set. It does not yet decide which targets are best for tree building. The extracted targets can support several downstream strategies: one tree per gene, one concatenated nif tree and one concatenated nod tree, or separate trees by sample group (`No`, `Rh`, `Ro`) versus one combined tree across biological samples.

**Version 2 change from the original plan:** Step 6 remains conceptually valid, but the downstream interpretation changes. We should no longer use `BLAN` as a blank-control group because Ryan clarified that BLAN is a real NEON field site. Instead, `MC-1` and `MC-2` should be used as negative controls in later QC steps. This means the target extraction remains the same, but target ranking and filtering should be revised in later steps to use true negative controls and to treat BLAN as a biological site.

**Next analysis goal after Step 6:** Build a revised phylogenetic pipeline that can answer Ryan's main questions:

```text
1. Which nif/nod targets are well supported across biological samples?
2. Do MC-1 and MC-2 show low/no nif/nod signal as true negative controls?
3. How many samples can be recovered if coverage/N filters are adjusted?
4. Can we build one tree per gene, one combined nif tree, and/or one combined nod tree?
5. Can the consensus step be adapted to detect multiple nifH/nif/nod sequence types per sample?
6. How sensitive is reference-guided mapping to divergent sequences?
```

