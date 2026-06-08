# Symbiosis Gene Analysis Pipeline Report

## Step 1. Identify the Correct Dataset for Functional-Gene Phylogenies

### Purpose

The first step was to choose the correct dataset for building `nif` and `nod` gene trees. Ryan's notes show that functional genes are in `symbiosis_sorted`, while 16S data are for bacterial community/taxonomic analysis.

### Connection to Ryan's Email

Ryan explicitly says:

- "Basically, we use a novel method for sequence capture of three separate locus classes simultaneously: 16S (full length), ITS (full length plus flanking 16S and 26S and 5.8S), and a symbiosis island dataset."

- "It takes a lot of compute to sort these loci out, so I thought I would send along pre-sorted subsets. These are in 16S_sorted, ITS_sorted, and symbiosis_sorted."

- "I have attached a genomic reference in case you are interested in the functional genes in symbiosis_sorted."

- "For 16S, if you need a reference I recommend SILVA. For ITS, I recommend UNITE, preferably the all_eukaryote version so it will be easier to detect host DNA."

- "I am still finalizing functional gene data, which are subject to a different pipeline. Claudia and I already discussed using this to derive phylogenetic trees but I wanted to raise it if others are interested in functional genes."

### Support From Relevant Literature

The two reference papers support the scientific motivation for this workflow.

Parker (2015) studied host use by *Bradyrhizobium* root-nodule bacteria across legume clades using bacterial phylogenies, including housekeeping genes and the symbiosis-island gene `nifD`. This supports using symbiosis-related functional genes to study host-symbiont specificity across phylogenetic scales.

Braga et al. (2018) showed that combining phylogenetic information with interaction-network analysis can help test hypotheses about host-associated diversification and host-use patterns. This supports the broader project goal of comparing bacterial phylogenies with plant-microbe interaction networks to study fuzzy specificity.

Together, these papers support the overall strategy of recovering bacterial functional-gene sequences, building gene trees, and later comparing those trees with host identity, plant phylogeny, geography, and interaction networks.

### Interpretation

Ryan separates the sequencing data into three different locus classes:

```text
16S_sorted         -> bacterial 16S/community/taxonomic analysis
ITS_sorted         -> fungal ITS/community/taxonomic analysis
symbiosis_sorted   -> symbiosis-island / functional-gene analysis
```

Therefore, the `nif` and `nod` phylogenetic analysis should use:

```text
symbiosis_sorted
symbiosis_islands.fasta
symbiosis_islands.gb
symbiosis_islands_gene_list.xlsx
```

not the 16S/SILVA workflow.

### Why I Did Not Continue With 16S/SILVA

I initially tried the 16S/SILVA workflow because Ryan recommended SILVA for 16S and the project notes included bacterial phylogeny goals. However, the 16S/SILVA consensus output was not suitable for `nif`/`nod` phylogenetic analysis.

The issue was that SILVA is a broad 16S taxonomic reference. Mixed-community reads mapped across many related 16S records, producing many mostly `N`-filled consensus sequences rather than one reliable sequence per sample.

Concise explanation:

```text
The 16S/SILVA workflow was stopped because SILVA is a broad 16S taxonomic reference; mixed-community reads mapped across many related 16S records and produced mostly N-filled consensus sequences, not reliable nif/nod gene sequences.
```

### Correct Dataset Used

The corrected pipeline uses Ryan's functional-gene dataset:

```bash
/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted/
```

Important files:

```text
raw_symbiosis_full/              # copied raw paired reads
symbiosis_islands.fasta          # Ryan's functional-gene mapping reference
symbiosis_islands.gb             # GenBank annotation file
symbiosis_islands_gene_list.xlsx # gene summary file
```

### Result of Step 1

We concluded that:

- 16S data should be used for bacterial community/taxonomic analyses.
- `symbiosis_sorted` is the correct dataset for functional-gene phylogenetic analysis.
- `symbiosis_islands.fasta` is the correct reference for mapping `symbiosis_sorted` reads.
- `symbiosis_islands.gb` and `symbiosis_islands_gene_list.xlsx` should be used to identify annotated `nif` and `nod` targets.
- The 16S/SILVA work should be described only as an earlier attempted workflow, not as part of the final corrected pipeline.

### Step 1 Final Decision

```text
Use symbiosis_sorted and Ryan's symbiosis-island reference/annotation files for the functional-gene pipeline.

Do not use the 16S/SILVA consensus output for nif/nod phylogenetic tree construction.
```


```text
Raw symbiosis_sorted reads
  ↓
fastp trimming using Ryan's settings
  ↓
Map trimmed reads to Ryan's symbiosis_islands.fasta
  ↓
Use symbiosis_islands.gb + gene list to identify nif/nod targets
  ↓
Calculate nif/nod gene-region coverage from mapped BAMs
  ↓
Compare biological samples with BLAN controls
  ↓
Select promising targets using blank-aware filtering
  ↓
Extract pilot nifH consensus sequences
  ↓
Filter low-quality consensus sequences
  ↓
MAFFT multiple sequence alignment
  ↓
IQ-TREE pilot nifH phylogeny
  ↓
iTOL tree visualization
```

## 2. Starting Point

Working folder:

```bash
/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted/
```
BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"

Initial organization:

```text
symbiosis_sorted/
├─ raw_symbiosis_full/                 # copied raw paired reads
├─ symbiosis_islands.fasta             # Ryan's symbiosis reference for mapping
├─ symbiosis_islands.gb                # gene annotations and coordinates
├─ symbiosis_islands_gene_list.xlsx    # summary of nif/nod genes
└─ Rscripts/                           # scripts used in this pipeline
```

Dataset size:

```text
1,116 paired samples total
12 BLAN controls
1,104 biological samples
```

See the data:
Each FASTQ read uses 4 lines:

```text
@read_id/R1orR2
sequence
+
quality (The many I characters usually indicate high base quality. The 9 and * characters indicate lower quality bases than I.)
```

```bash
zcat "$BASE/raw_symbiosis_full/TALL-13-3-Ro_R2.fq.gz" | head -20
```

```text
@LH00516:359:23333FLT3:3:1101:1536:25596/2
CCTTGGAGATGTTGTT...
+
IIIIIIII9IIIIIIII...
```


# 3. Step-by-Step Pipeline

## 3.1 Read Quality Control and Trimming

### Purpose

Clean raw paired-end reads before mapping by removing adapters, low-quality bases, poly-G/poly-X artifacts, low-complexity reads, and reads shorter than 36 bp.

### Ryan Connection

Ryan recommended using `fastp` and provided the trimming settings in the project notes. This was important because Ryan noted that NovaSeq data can contain problematic poly-G read-through artifacts. We therefore used Ryan's `fastp` settings for the full August 2025 `symbiosis_sorted` dataset.

### Input

```bash
$BASE/raw_symbiosis_full/*_R1.fq.gz
$BASE/raw_symbiosis_full/*_R2.fq.gz
```

where:

```bash
BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"
```

### Script Used for Trimming

```bash
bash "$BASE/Rscripts/trim_symbiosis_full_ryan.sh"
```

This script runs `fastp` on each paired sample using Ryan's settings. For each sample, it creates trimmed paired reads, unpaired reads, HTML/JSON QC reports, and a trimming manifest.

### Folder Organization for This Step

```text
symbiosis_sorted/
├─ raw_symbiosis_full/
│  ├─ *_R1.fq.gz                         # raw forward reads
│  └─ *_R2.fq.gz                         # raw reverse reads
│
├─ symbiosis_trimmed_ryan_full/
│  ├─ *_P1.fastq.gz                      # trimmed paired Read 1
│  ├─ *_P2.fastq.gz                      # trimmed paired Read 2
│  ├─ *_U1.fastq.gz                      # unpaired Read 1 after trimming
│  └─ *_U2.fastq.gz                      # unpaired Read 2 after trimming
│
├─ fastp_logs_ryan_full/
│  ├─ *.fastp.html                       # per-sample fastp HTML QC report
│  ├─ *.fastp.json                       # per-sample fastp JSON QC report
│  └─ *.fastp.stderr.log                 # fastp error/output log per sample
│
├─ trimmed_fastp_QC_checking/
│  ├─ fastp_logs_ryan_full_summary.tsv
│  ├─ fastp_qc_summary_run.log
│  ├─ fastp_quality_before_trimming_all_samples_lightpurple_mean_blue.svg
│  └─ fastp_quality_after_trimming_all_samples_lightpurple_mean_blue.svg
│
└─ Rscripts/
   ├─ trim_symbiosis_full_ryan.sh (Main trimming script)
   ├─ summarize_fastp_qc.sh 
   ├─ make_fastp_quality_profile_figures.sh
   └─ make_fastp_quality_before_all_samples.sh
```

### Trimming Output Example

For sample `TALL-13-3-Ro`, `fastp` produced four output files:

```text
TALL-13-3-Ro_P1.fastq.gz
TALL-13-3-Ro_P2.fastq.gz
TALL-13-3-Ro_U1.fastq.gz
TALL-13-3-Ro_U2.fastq.gz
```

The `P1/P2` files are the trimmed paired reads and were used for downstream paired-end mapping. The `U1/U2` files are reads whose mate was removed during trimming; these were kept as trimming outputs but were not used for the main paired-end mapping step.

Example commands to view the trimmed reads:

```bash
zcat "$BASE/symbiosis_trimmed_ryan_full/TALL-13-3-Ro_P1.fastq.gz" | head -12
zcat "$BASE/symbiosis_trimmed_ryan_full/TALL-13-3-Ro_P2.fastq.gz" | head -12
zcat "$BASE/symbiosis_trimmed_ryan_full/TALL-13-3-Ro_U1.fastq.gz" | head -12
zcat "$BASE/symbiosis_trimmed_ryan_full/TALL-13-3-Ro_U2.fastq.gz" | head -12
```

Each FASTQ record has the standard four-line format:

```text
@read_id
DNA sequence
+
base-quality scores
```

### Example Read Counts

For sample `TALL-13-3-Ro`, the raw files contained:

```bash
zcat "$BASE/raw_symbiosis_full/TALL-13-3-Ro_R1.fq.gz" | awk 'END {print NR/4}'
zcat "$BASE/raw_symbiosis_full/TALL-13-3-Ro_R2.fq.gz" | awk 'END {print NR/4}'
```

Result:

```text
R1 raw reads: 123,581
R2 raw reads: 123,581
```

After trimming, the paired reads were:

```bash
zcat "$BASE/symbiosis_trimmed_ryan_full/TALL-13-3-Ro_P1.fastq.gz" | awk 'END {print NR/4}'
zcat "$BASE/symbiosis_trimmed_ryan_full/TALL-13-3-Ro_P2.fastq.gz" | awk 'END {print NR/4}'
```

Result:

```text
P1 trimmed paired reads: 118,873
P2 trimmed paired reads: 118,873
```

Unpaired reads after trimming:

```bash
zcat "$BASE/symbiosis_trimmed_ryan_full/TALL-13-3-Ro_U1.fastq.gz" | awk 'END {print NR/4}'
zcat "$BASE/symbiosis_trimmed_ryan_full/TALL-13-3-Ro_U2.fastq.gz" | awk 'END {print NR/4}'
```

Result:

```text
U1 unpaired reads: 2,681
U2 unpaired reads: 1,075
```

### Trimming Completion Check

The trimming manifest records whether each sample completed successfully:

```bash
tail -n +2 "$BASE/symbiosis_trimmed_ryan_full_manifest.tsv" | cut -f2 | sort | uniq -c
```

Result summary:

```text
1,116 paired samples completed
0 failed samples
0 missing R2 samples
```

Output files generated:

```text
1,116 P1 files
1,116 P2 files
1,116 fastp JSON reports
```

BLAN control samples were retained for later QC.

### Quality-Control Summary

Quality-control metrics were summarized from the `fastp` JSON reports using:
### 1. Summarize fastp JSON reports into a QC table and log

```bash
bash "$BASE/Rscripts/summarize_fastp_qc.sh" \
  | tee "$BASE/trimmed_fastp_QC_checking/fastp_qc_summary_run.log"
```

The QC summary table contains per-sample read counts, read-retention percentage, Q20/Q30 rates, and GC content before and after trimming.

### Additional Quality-Profile Figures

Two additional scripts were run to make clearer report-ready per-base quality figures:
### 2. Make before-trimming all-sample quality figure
```bash
bash "$BASE/Rscripts/make_fastp_quality_before_all_samples.sh" \
  | tee "$BASE/trimmed_fastp_QC_checking/fastp_quality_before_all_samples_run.log"
```

### 3. Make after-trimming and mean-only quality figures

```
bash "$BASE/Rscripts/make_fastp_quality_profile_figures.sh" \
  | tee "$BASE/trimmed_fastp_QC_checking/fastp_quality_profile_figures_run.log"
``` 
The result of this step is in this path:
[https://github.com/solislemuslab/IntBio-NitFix/main/Results/august2025/symbiosis_sorted](https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/QC%20Check)


### 3.2. Mapping Reads to Ryan's Symbiosis-Island Reference

**Purpose:** Align cleaned functional-gene reads to the reference provided by Ryan.

**Ryan connection:** Ryan stated that `symbiosis_sorted` contains functional-gene data and provided `symbiosis_islands.fasta` as the genomic reference for these reads. Mapping to this reference answers which reads support the provided symbiosis-island and functional-gene regions.

**Reference used:** `$BASE/symbiosis_islands.fasta`

**Input:** `$BASE/symbiosis_trimmed_ryan_full/`, `$BASE/symbiosis_islands.fasta`

**Output:** `$BASE/symbiosis_mapped_full/*.bam`, `$BASE/symbiosis_mapped_full/*.bam.bai`, `$BASE/symbiosis_mapping_logs_full/*.flagstat.txt`, `$BASE/symbiosis_mapping_full_manifest.tsv`

**Script/command used:**

```bash
bash $BASE/Rscripts/map_symbiosis_full.sh
```

**What the script does:** Indexes the FASTA if needed, maps reads with `bwa mem`, sorts/indexes BAM files with `samtools`, and writes `flagstat` reports.

**Samples:** 1,116.

**Runtime:** Completed within the same day; much faster than the earlier 16S/SILVA mapping.

**Result:** All samples mapped successfully: 1,116 BAM files, 1,116 indexes, 1,116 flagstat reports, 0 failed samples, and 0 zero-size BAMs.

**Example of checking mapping summary:** 
```bash
module load samtools-1.9
BASE="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/august2025/symbiosis_sorted"
BAM="$BASE/symbiosis_mapped_full/TALL-13-3-Ro.bam"

ls -lh "$BAM" "$BAM.bai"
samtools flagstat "$BAM"
```
```text
239767 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
2021 + 0 supplementary
0 + 0 duplicates
229655 + 0 mapped (95.78% : N/A)
237746 + 0 paired in sequencing
118873 + 0 read1
118873 + 0 read2
214626 + 0 properly paired (90.28% : N/A)
218190 + 0 with itself and mate mapped
9444 + 0 singletons (3.97% : N/A)
3550 + 0 with mate mapped to a different chr
3367 + 0 with mate mapped to a different chr (mapQ>=5)

Total reads in BAM: 239,767
Mapped reads: 229,655
Mapped percent: 95.78%
Paired reads: 237,746
Properly paired reads: 214,626
Properly paired percent: 90.28%
Singletons: 9,444

Important connection to trimming:
P1 reads after trimming: 118,873
P2 reads after trimming: 118,873
Total paired reads entering mapping: 237,746
```

For example, sample `TALL-13-3-Ro` had 118,873 trimmed paired reads in each direction (`P1` and `P2`), giving 237,746 paired reads entering the mapping step. After mapping to `symbiosis_islands.fasta`, `samtools flagstat` reported 229,655 mapped reads (95.78%) and 214,626 properly paired reads (90.28%). This confirms that the trimmed reads mapped successfully to Ryan’s symbiosis-island reference.



**See which reference regions got the most mapped reads:**

The mapping reference for this step was symbiosis_islands.fasta, which contains multiple symbiosis-island and functional-gene DNA reference sequences. Each sequence in this FASTA file has a header name, such as NC_009937_Symbiosis_Island_4_Azorhizobium_caulinodans_ORS_571, followed by the actual DNA sequence. During mapping, reads are aligned to the DNA sequences, but the BAM file reports the matching reference by its FASTA header name. For the example sample TALL-13-3-Ro, samtools idxstats was used to summarize how many reads mapped to each reference region. The output columns are reference name, reference length in base pairs, number of mapped reads, and number of unmapped reads. The results were sorted by mapped-read count, showing that the strongest supported region for this sample was NC_009937_Symbiosis_Island_4_Azorhizobium_caulinodans_ORS_571, which is 8,188 bp long and had 83,325 mapped reads. This confirms that reads from TALL-13-3-Ro mapped strongly to symbiosis-island reference regions. These values are raw read counts for this one sample, so they are useful for inspection but are not normalized by reference length.

```bash
samtools idxstats "$BAM" | sort -k3,3nr | head -20
```
```text
NC_009937_Symbiosis_Island_4_Azorhizobium_caulinodans_ORS_571	8188	83325	1269
```
**Example BAM Alignment Inspection for `TALL-13-3-Ro`**
```bash
samtools view "$BAM" | head -5
```

**Example Coverage Inspection on the Top Mapped Reference**
```text
Coverage was checked for the top mapped reference in sample `TALL-13-3-Ro` using `samtools depth`. The top reference was `NC_009937_Symbiosis_Island_4_Azorhizobium_caulinodans_ORS_571`. In the output, the columns are reference name, base position, and read depth. For example, position 68 had depth 1, meaning one read covered that base, and position 69 had depth 2, meaning two reads covered that base. This confirms that mapped reads cover specific positions along the reference sequence.
```

```bash
TOPREF=$(samtools idxstats "$BAM" | sort -k3,3nr | head -1 | cut -f1)

echo "$TOPREF"

samtools depth -r "$TOPREF" "$BAM" | head -20
```

```text
samtools depth -r "$TOPREF" "$BAM" | head -20
NC_009937_Symbiosis_Island_4_Azorhizobium_caulinodans_ORS_571	68	1
NC_009937_Symbiosis_Island_4_Azorhizobium_caulinodans_ORS_571	69	2
```


Mapping quality was evaluated across all 1,116 samples using samtools flagstat summaries, including mapped-read percentage and properly paired-read percentage after mapping to Ryan’s symbiosis_islands.fasta reference. The mapping QC results are located in
https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/mapping_QC_checking
, including the summary table, completion check, run log, and mapping-quality figure.


### 3.3 Extract Central `nif`/`nod` Targets

**Purpose:** Identify which central `nif` and `nod` genes are present in Ryan's GenBank annotation.
Mapping used the full reference to place reads broadly on Ryan’s symbiosis-island sequences, but downstream **phylogenetic analysis needs gene-specific targets**, so we used the GenBank annotation to identify and extract the exact nif/nod gene regions.
```text
which reads belong to nif genes
which reads belong to nod genes
which gene regions have enough coverage
which consensus sequences should be used for nif/nod trees
```
**Ryan connection:** Ryan supplied `symbiosis_islands.gb` and the gene list, so the annotation file was used to identify available functional-gene targets.

**Input:** `$BASE/symbiosis_islands.gb`, `$BASE/symbiosis_islands_gene_list.xlsx`

**Gene-list summary-symbiosis_islands_gene_list.xlsx:**  
The `symbiosis_islands_gene_list.xlsx` file summarizes the annotated genes in Ryan’s symbiosis-island reference. According to the summary sheet, the reference contains 73 unique genes, including 18 unique central nitrogen-fixation genes, 17 unique central nodulation genes, and 38 accessory genes. For this step, we focused on the central `nif` and `nod` genes because the project goal is to build functional-gene phylogenies for nitrogen-fixation and nodulation genes.

The central genes listed in the spreadsheet were:

```text
nif genes: nifA, nifB, nifD, nifE, nifH, nifJ, nifK, nifM, nifN, nifQ, nifS, nifT, nifU, nifV, nifW, nifX, nifY, nifZ

nod genes: nodA, nodB, nodC, nodD, nodE, nodF, nodH, nodI, nodJ, nodL, nodP, nodQ, nodS, nodT, nodU, nodX, nodZ

**symbiosis_islands.gb**
```text
where each gene starts
where each gene ends
which strand it is on
what gene name it has
```
So the script used the GenBank annotation to extract the DNA regions corresponding to central nif and nod genes.
**Output:** `$BASE/nif_nod_target_reference/`

**Script/command used:**

```bash
python3 $BASE/Rscripts/extract_all_nif_nod_targets.py
```

**What the script does:** Parses the GenBank file, extracts annotated central `nif`/`nod` CDS sequences, removes exact duplicates, and writes FASTA/summary tables.

**Runtime:** Less than 1 minute.

**Result:** Ryan’s `symbiosis_islands_gene_list.xlsx` identified 18 central `nif` genes and 17 central `nod` genes, and `symbiosis_islands.gb` was used to extract the matching annotated gene sequences.

From the 369 annotated gene/CDS regions in `symbiosis_islands.gb`, we extracted **231 unique central `nif`/`nod` target regions**. These included **169 regions annotated as central `nif` genes** and **62 regions annotated as central `nod` genes**. These numbers are reference gene-region counts, not sample counts. Six genes from the spreadsheet were not found in the GenBank annotation (`nifM`, `nifY`, `nodE`, `nodF`, `nodP`, and `nodT`), so they were recorded as missing and excluded from downstream target analysis.

results are located in https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_target_reference


### 3.4 Blank-Control Evaluation

**Purpose:** Check whether BLAN control samples contain reads that map to the extracted central `nif`/`nod` target reference.

**Ryan connection:** BLAN samples are control samples, not biological plant samples. They were retained during trimming and mapping so that background signal, possible contamination, or technical carryover could be evaluated before interpreting biological `nif`/`nod` signal.

**Input:**

```bash
$BASE/symbiosis_trimmed_ryan_full/BLAN*_P1.fastq.gz
$BASE/symbiosis_trimmed_ryan_full/BLAN*_P2.fastq.gz
$BASE/nif_nod_target_reference/all_central_nif_nod_targets.fasta
```

**Output:** `$BASE/nif_nod_blank_mapping/`, `$BASE/nif_nod_blank_mapping_logs/`, `$BASE/nif_nod_blank_mapping_summary.tsv`

**Script/command used:**

```bash
bash $BASE/Rscripts/map_all_blank_nif_nod_targets.sh
```

**What the script does:** 
The script maps the 12 BLAN control samples to the extracted central nif/nod target FASTA, generates BAM files and mapping logs, and summarizes total reads, mapped reads, mapped percentage, properly paired reads, and the top mapped target for each BLAN sample.

**Samples:** 12 BLAN controls.

**Runtime:** Fast; completed in the same session.

**Result:** All 12 BLAN controls showed substantial mapping to the extracted nif/nod target reference, with mapped-read percentages ranging from 47.35% to 92.08%. This indicates that the BLAN controls contain target-associated background signal and must be considered during QC. Therefore, BLAN samples were retained for quality-control evaluation but excluded from biological consensus-sequence extraction and phylogenetic tree construction.

Result of nif_nod_blank_mapping_summary.tsv locaed at:
https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_blank_mapping

### Step 5. Match Targets Back to Ryan's Original FASTA

**Purpose:** Confirm that extracted `nif`/`nod` targets are present in the original reference used for mapping.

**Ryan connection:** This keeps downstream coverage analysis tied to Ryan's original `symbiosis_islands.fasta`, rather than replacing the main mapping reference.

**Input:** `$BASE/symbiosis_islands.fasta`, `$BASE/nif_nod_target_reference/all_central_nif_nod_targets.fasta`

**Output:** `$BASE/nif_nod_original_reference_regions/`

**Script/command used:**

```bash
python3 $BASE/Rscripts/match_nif_nod_targets_to_original_reference.py
```

**What the script does:** Searches target sequences against Ryan's FASTA and records exact reference coordinates in TSV/BED format.

**Runtime:** Less than 1 minute.

**Result:** All 231 extracted target sequences matched Ryan's original FASTA. A total of 233 exact locations were found because `nifQ` and `nifW` had multiple exact reference locations.

Results are located at:

https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_original_reference_regions

### 3.6. Gene-Level Coverage From Existing BAMs

**Purpose:** Measure coverage of each matched `nif`/`nod` region in each mapped sample.

**Ryan connection:** This uses Ryan's original mapping reference and annotation files to determine which functional genes are supported enough for phylogenetic analysis.

**Input:** `$BASE/symbiosis_mapped_full/*.bam`, `$BASE/nif_nod_original_reference_regions/nif_nod_matches_in_original_reference.tsv`

**Output:** `$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv`

**Script/command used:**

```bash
bash $BASE/Rscripts/calculate_nif_nod_region_coverage.sh
```

**What the script does:** Uses `samtools depth` to calculate percent covered, mean depth, max depth, and covered bases for every sample/region.

**Samples/regions:** 1,116 samples x 233 target regions.

**Runtime:** Completed in the same session

**Result:** Final table had 260,029 lines: 1 header plus 260,028 sample-region rows.


### 3.7. Coverage Summary by Sample Type

**Purpose:** Summarize coverage by gene and sample group before consensus extraction.

**Input:** `$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv`

**Output:** `$BASE/nif_nod_coverage_existing_mapping/nif_nod_gene_coverage_summary_by_sample_type.tsv`, `$BASE/nif_nod_coverage_existing_mapping/nif_nod_target_coverage_summary.tsv`, `$BASE/nif_nod_coverage_existing_mapping/nif_nod_sample_coverage_summary.tsv`

**Script/command used:**

```bash
python3 $BASE/Rscripts/summarize_nif_nod_coverage.py
```

**What the script does:** Summarizes good coverage using thresholds `percent_covered >= 80` and `mean_depth >= 10` across BLAN, No, Rh, and Ro groups.

**Runtime:** Less than 1 minute.

**Result:** Several `nif` genes had strong nodule support, including `nifH`, `nifU`, `nifD`, `nifK`, `nifE`, and `nifB`. Many BLAN samples also showed strong signal, so blank-aware filtering was required.




### Step 8. Blank-Aware Target Ranking

**Purpose:** Identify targets with strong nodule signal and lower blank-control background.

**Input:** `$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv`

**Output:** `$BASE/nif_nod_coverage_existing_mapping/nif_nod_target_blank_aware_ranking_v2.tsv`

**Script/command used:**

```bash
python3 $BASE/Rscripts/rank_nif_nod_targets_blank_aware_v2.py
```

**What the script does:** Classifies each target as `PROMISING`, `BIOLOGICAL_SIGNAL_BUT_BLANK_BACKGROUND`, `LOW_OR_LIMITED_SIGNAL`, or `NOT_SUPPORTED`.

**Runtime:** Less than 1 minute.

**Result:** 21 targets were `PROMISING`. With a stricter cutoff of at least 75 good nodule samples, 10 promising targets remained; all were `nif` genes. No `nod` target passed this stricter filter.

### Step 9. Pilot Consensus Extraction for Best `nifH` Target

**Purpose:** Test whether a strong target can produce usable sample-level consensus sequences.

**Target:** `nifH|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island|ref63`

**Input:** `$BASE/symbiosis_islands.fasta`, `$BASE/symbiosis_mapped_full/*.bam`, coverage table from Step 6

**Output:** `$BASE/nifH_ref63_pilot_consensus_v2/`

**Script/command used:**

```bash
bash $BASE/Rscripts/extract_pilot_nifH_ref63_consensus_v2.sh
```

**What the script does:** Selects good nodule samples, calls variants with `bcftools`, applies variants with `bcftools consensus`, masks zero-depth sites as `N`, and writes FASTA/QC outputs.

**Runtime:** Fast; completed in the same session.

**Result:** 159 full-length 894 bp consensus sequences were recovered. After filtering for `N_percent <= 5`, 104 high-quality sequences remained.

### Step 10. Alignment, Tree Construction, and Visualization

**Purpose:** Build a pilot `nifH` phylogenetic tree.

**Input:** `$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_Nle5.fasta`

**Output:** `$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_Nle5.mafft.fasta`, `$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_Nle5_iqtree.treefile`

**Commands used:**

```bash
mafft --auto "$FILTERED" > "$OUTDIR/nifH_ref63_consensus_Nle5.mafft.fasta"

iqtree -s "$ALIGN" -m MFP -B 1000 -T AUTO \
  --prefix "$OUTDIR/nifH_ref63_Nle5_iqtree"
```

**What the tools do:** MAFFT aligns the filtered consensus sequences. IQ-TREE performs model selection and maximum-likelihood tree inference with ultrafast bootstrap support.

**Runtime:** MAFFT was fast. IQ-TREE completed the pilot tree in the same terminal session.

**Result:** The alignment contained 104 sequences and was 894 bp long. IQ-TREE produced a pilot tree with best-fit model `GTR+F+I+G4`. The tree was uploaded and visualized in iTOL.

## 4. Important Results to Copy to Laptop

| Result type | File/folder | Why important |
|---|---|---|
| Gene summary | `$BASE/nif_nod_target_reference/central_nif_nod_gene_summary.tsv` | Shows available `nif`/`nod` genes |
| Mapping manifest | `$BASE/symbiosis_mapping_full_manifest.tsv` | Confirms all samples mapped |
| Coverage table | `$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv` | Main region-level coverage result |
| Coverage summaries | `$BASE/nif_nod_coverage_existing_mapping/nif_nod_gene_coverage_summary_by_sample_type.tsv` | Compares BLAN/No/Rh/Ro |
| Target ranking | `$BASE/nif_nod_coverage_existing_mapping/nif_nod_target_blank_aware_ranking_v2.tsv` | Lists promising targets |
| Filtered consensus FASTA | `$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_Nle5.fasta` | Final pilot sequences |
| MAFFT alignment | `$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_Nle5.mafft.fasta` | Alignment used for tree |
| IQ-TREE tree | `$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_Nle5_iqtree.treefile` | Main pilot phylogeny |
| IQ-TREE report | `$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_Nle5_iqtree.iqtree` | Model and run details |
| Tree figure | iTOL export | Report-ready visualization |

## 5. Final Folder Organization

```text
symbiosis_sorted/
├─ raw_symbiosis_full/
├─ symbiosis_trimmed_ryan_full/
├─ fastp_logs_ryan_full/
├─ symbiosis_islands.fasta
├─ symbiosis_islands.gb
├─ symbiosis_islands_gene_list.xlsx
├─ symbiosis_mapped_full/
├─ symbiosis_mapping_logs_full/
├─ nif_nod_target_reference/
├─ nif_nod_blank_mapping/
├─ nif_nod_original_reference_regions/
├─ nif_nod_coverage_existing_mapping/
├─ nifH_ref63_pilot_consensus_v2/
│  ├─ sequences/
│  ├─ logs/
│  ├─ nifH_ref63_consensus_Nle5.fasta
│  ├─ nifH_ref63_consensus_Nle5.mafft.fasta
│  └─ nifH_ref63_Nle5_iqtree.treefile
└─ Rscripts/
```

## 6. Final Notes

This pipeline successfully processed Ryan's `symbiosis_sorted` functional-gene data through trimming, mapping, annotation-based target identification, coverage analysis, blank-aware filtering, pilot consensus extraction, alignment, and pilot tree construction. The most important final outputs are the blank-aware target ranking, the filtered `nifH` consensus FASTA, the MAFFT alignment, and the IQ-TREE `.treefile`.

The main caution is that BLAN controls showed strong target-associated signal, so final gene selection must continue to use blank-aware filtering. The next recommended step is to repeat consensus extraction and tree construction for the other strong `nif` targets, then compare trees before deciding whether to build individual or concatenated `nif` phylogenies.
