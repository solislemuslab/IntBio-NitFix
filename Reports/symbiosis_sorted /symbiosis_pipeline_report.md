# Symbiosis Gene Analysis Pipeline Report

## Summary of Ryan’s Main Idea

Ryan’s dataset was designed to study plant-microbe specificity in nitrogen-fixing plant systems across many NEON sites. The samples include rhizosphere, root, and nodule material from many host plants, so the data can be used to compare microbial communities across tissue types, sites, environments, and plant lineages. Ryan provided three main kinds of sequencing data because they answer different parts of the project: `16S_sorted` describes bacterial community composition, `ITS_sorted` describes fungal/community composition, and `symbiosis_sorted` contains functional symbiosis-gene reads. 

The broader project goal is to build and compare three bacterial phylogenies: a core genome tree, a `nif` gene tree, and a `nod` gene tree. These trees can then be compared with plant phylogeny, plant identity, microbial community composition, environmental covariates, and plant-microbe interaction networks. The main biological question is whether related plants associate with related bacteria, but not necessarily with the exact same bacterial strains. This is the idea of fuzzy specificity: host-symbiont associations may be conserved at broader bacterial clade levels rather than at strict strain or species levels. 


In this pipeline, I focused on the `symbiosis_sorted` data because it is the correct dataset for recovering `nif`/`nod` functional-gene sequences and building functional-gene trees. The 16S and ITS datasets remain important for community composition, diversity metrics, and interaction-network analyses, but they are not the correct input for the `nif`/`nod` phylogenetic pipeline.


## Summary of the Parker Paper

Parker (2015), “The Spread of *Bradyrhizobium* Lineages Across Host Legume Clades: from *Abarema* to *Zygia*,” studied how *Bradyrhizobium* root-nodule bacteria are distributed across many legume host lineages. The paper used bacterial phylogenetic information, including housekeeping genes and the symbiosis-island gene `nifD`, to ask whether related legumes tend to host related bacterial symbionts. The main result was that many bacterial lineages were found across broad and diverse legume host groups, showing that host use is not strictly one plant lineage to one bacterial lineage. However, the paper also found evidence that some related host plant clades contained more similar bacterial symbionts than expected by chance. This supports Ryan’s idea of fuzzy specificity: plant-microbe associations may be structured at broader clade levels, even when they are not perfectly specific at the species or strain level. This paper motivates our goal of building bacterial gene trees and comparing them with plant phylogeny and interaction networks.

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

**Ryan connection:** Ryan defined the biological sample types as `rh = rhizosphere`, `ro = root`, and `no = nodule`. 
Samples beginning with `BLAN` were treated as blank controls based on their names and were retained during trimming and mapping so that background signal, possible contamination, or technical carryover could be evaluated before interpreting biological `nif`/`nod` signal.

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

### 3.5. Match Targets Back to Ryan's Original FASTA

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

**Samples/regions:** 1,116 samples x 233 are exact matched locations in Ryan’s original reference.

**Runtime:** fast

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


Summary tables and the coverage heatmap for this step are available in the GitHub results folder:  
https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_coverage_existing_mapping

### 3.8. Blank-Aware Target Ranking

**Purpose:** Identify targets with strong nodule signal and lower blank-control background.

**Input:** `$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv`

**Output:** `$BASE/nif_nod_coverage_existing_mapping/nif_nod_target_blank_aware_ranking_v2.tsv`

**Script/command used:**

```bash
python3 $BASE/Rscripts/rank_nif_nod_targets_blank_aware_v2.py
```

**What the script does:** Classifies each target as `PROMISING`, `BIOLOGICAL_SIGNAL_BUT_BLANK_BACKGROUND`, `LOW_OR_LIMITED_SIGNAL`, or `NOT_SUPPORTED`.

**Runtime:** Less than 1 minute.


**Result:** The blank-aware ranking classified 21 targets as `PROMISING`. When an additional stricter filter was applied requiring at least 75 good nodule (`No`) samples, 10 promising targets remained, all annotated as `nif` genes. No `nod` target passed this stricter pilot filter. This does not mean `nod` genes are absent or unimportant; rather, under the current strict coverage and blank-background thresholds, `nod` targets did not have enough clean support for the first pilot tree. Because `nod` genes are central to the project goals, they will be evaluated separately using a `nod`-focused target-selection strategy.

nif_nod_target_blank_aware_ranking_v2.tsv is available in the GitHub results folder:  
https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nif_nod_coverage_existing_mapping

### 3.9 Pilot Consensus Extraction for Best `nifH` Target

**Purpose:** Test whether one strong blank-aware target can produce usable sample-level consensus sequences for phylogenetic analysis.

**Target:** `nifH|NZ_CP141049|NZ_CP141049_-_Symbiosis_Island|ref63`

**Why this target was selected:** This `nifH` target was one of the strongest blank-aware candidates. It had good coverage in 159 nodule (`No`) samples, only 1 BLAN control with good coverage, high nodule median depth, and a strong nodule-vs-blank depth ratio. Because `nifH` is also a key nitrogen-fixation gene, it was selected as the first pilot target.

**Input:**

```bash
$BASE/symbiosis_islands.fasta
$BASE/symbiosis_mapped_full/*.bam
$BASE/nif_nod_coverage_existing_mapping/nif_nod_region_coverage_all_samples.tsv
$BASE/nif_nod_original_reference_regions/nif_nod_matches_in_original_reference.tsv

**Output:**
$BASE/nifH_ref63_pilot_consensus_v2/

**Script/command used:**

```bash
bash $BASE/Rscripts/extract_pilot_nifH_ref63_consensus_v2.sh
```

**What the script does:** Selects good nodule samples, calls variants with `bcftools`, applies variants with `bcftools consensus`, masks zero-depth sites as `N`, and writes FASTA/QC outputs.


**Runtime:** Fast;

**Result:** The pilot extraction recovered 159 full-length 894 bp consensus sequences for the selected `nifH` target. After filtering for `N_percent <= 5`, 104 high-quality consensus sequences remained for alignment and pilot tree construction.

Results are available at:
https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nifH_ref63_pilot_consensus_v2

### 3.10 Alignment, Tree Construction, and Visualization

**Purpose:** Build a pilot `nifH` phylogenetic tree from the filtered high-quality consensus sequences.

**Input:**

```bash
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_Nle5.fasta


**Output:**

$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_consensus_Nle5.mafft.fasta
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_Nle5_iqtree.treefile
$BASE/nifH_ref63_pilot_consensus_v2/nifH_ref63_Nle5_iqtree.iqtree

**Commands used:**

```bash
mafft --auto "$FILTERED" > "$OUTDIR/nifH_ref63_consensus_Nle5.mafft.fasta"

iqtree -s "$ALIGN" -m MFP -B 1000 -T AUTO \
  --prefix "$OUTDIR/nifH_ref63_Nle5_iqtree"
```

**What the tools do:** MAFFT aligns the filtered consensus sequences so homologous nucleotide positions are compared across samples. IQ-TREE performs model selection and maximum-likelihood tree inference with ultrafast bootstrap support.

**Runtime:** MAFFT was very fast. IQ-TREE was fast.

**Result:** The filtered `nifH ref63` consensus set contained 104 high-quality sequences after applying the `N_percent <= 5` filter. MAFFT produced an alignment with 104 sequences and a length of 894 bp. IQ-TREE produced a pilot maximum-likelihood tree and selected `GTR+F+I+G4` as the best-fit model. The tree was uploaded and visualized in iTOL.

Results are available at:  
https://github.com/solislemuslab/IntBio-NitFix/tree/main/Results/august2025/symbiosis_sorted/nifH_ref63_pilot_consensus_v2


Visualization: Go to:

https://itol.embl.de/

Then:
```text
Click Upload a new tree.
Choose the file:
nifH_ref63_Nle5_iqtree.treefile
```


## 4. 
Do a nod-focused target review
Purpose: find the best nod candidate(s) for a separate nod tree


## 5. Final Notes

This pipeline successfully processed Ryan's `symbiosis_sorted` functional-gene data through trimming, mapping, annotation-based target identification, coverage analysis, blank-aware filtering, pilot consensus extraction, alignment, and pilot tree construction. The most important final outputs are the blank-aware target ranking, the filtered `nifH` consensus FASTA, the MAFFT alignment, and the IQ-TREE `.treefile`.

The main caution is that BLAN controls showed strong target-associated signal, so final gene selection must continue to use blank-aware filtering. The next recommended step is to repeat consensus extraction and tree construction for the other strong `nif` targets, then compare trees before deciding whether to build individual or concatenated `nif` phylogenies.





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



**Literature Review for step 12**
Functional nitrogen-fixation genes, especially nifH, have been widely used to study the diversity and relationships of nitrogen-fixing bacteria. The nifH gene encodes part of the nitrogenase enzyme complex and is commonly used as a marker for diazotrophs because it is conserved enough to align across organisms but variable enough to reveal genetic differences among lineages. For this reason, many studies build single-gene nifH phylogenies where the tree tips are sample-derived sequences, environmental clones, isolates, or consensus sequences, rather than different genes.
Ryan’s notes directly support this approach for our project. He identified symbiosis_sorted as the dataset containing the functional-gene reads, provided symbiosis_islands.fasta as the reference for these reads, and noted that the functional-gene data could be used to derive phylogenetic trees. Therefore, our workflow follows the intended structure of the dataset: trim the symbiosis_sorted reads, map them to Ryan’s symbiosis reference, identify well-covered nif/nod target regions, and build a pilot functional-gene tree from a strong nifH target.
The 2012 mangrove rhizosphere study provides a close example of this logic. Liu et al. analyzed nitrogen-fixing bacterial diversity using nifH sequences recovered from mangrove rhizosphere soil. Their method was to amplify and sequence nifH, align the nifH nucleotide sequences with ClustalW, and build neighbor-joining phylogenetic trees in MEGA with 1000 bootstrap replicates. The tips in their tree were not different genes; they were nifH sequences from environmental clones, cultured isolates, and GenBank references. Their result was that nifH sequences clustered into multiple groups, and sequences from specific sites tended to cluster together. This supports our interpretation that a single-gene nifH tree can be used to study nitrogen-fixing bacterial diversity and sample-level structure.
Older papers are still important because they established nifH as a standard marker gene. For example, Zehr et al. used amplified nifH sequences from environmental samples to detect previously unrecognized nitrogen-fixing microorganisms in ocean systems. Vinuesa et al. used nifH together with other genes in rhizobial phylogenetic/systematic analysis. These older studies are foundational, not outdated; they explain why nifH became a standard marker for nitrogen-fixation work.
Recent studies show that this field is still active, but the methods have expanded. Many modern papers now combine nifH marker-gene information with metagenomics, genomes, phylogenomics, imaging, or symbiosis biology. For example, Coale et al. 2024 reported a nitrogen-fixing organelle in a marine alga, showing that nitrogen-fixing symbioses remain an active and important research area. Kantor et al. 2024 used metagenomics to study genetic diversity among UCYN-A sublineages and their algal hosts. These recent studies show the modern direction of the field: nifH remains useful, but genome-scale data can provide stronger organism-level resolution when available.
Our tree should therefore be described carefully. It is not a whole-genome species tree and it is not a tree of different genes. It is a pilot single-locus nifH gene tree. Each tip represents one nodule sample’s nifH consensus sequence from the selected nifH ref63 target region. This tree shows that Ryan’s symbiosis_sorted data can recover high-quality sample-level nitrogen-fixation gene sequences and that those sequences contain enough variation to form phylogenetic structure. The next biological interpretation step is to color the tips by metadata, such as host plant, site, or sample code, and test whether well-supported clades correspond to biological groups.
References To Include
Ryan’s project notes: support use of symbiosis_sorted, symbiosis_islands.fasta, and functional-gene phylogenetic trees.
Liu et al. 2012: nifH phylogenetic diversity from mangrove rhizosphere soil; tree tips are nifH clone/isolate/reference sequences.
Zehr et al. 1998: environmental nifH amplification used to detect nitrogen-fixing microorganisms.
Vinuesa et al. 2005: rhizobial phylogenies using nifH with other genes.
Coale et al. 2024: nitrogen-fixing organelle in a marine alga.
Kantor et al. 2024: metagenomics of UCYN-A sublineage diversity and algal hosts.

Liu et al. 2012
https://doi.org/10.1139/W2012-016

Zehr et al. 1998
https://doi.org/10.1128/AEM.64.9.3444-3450.1998

Vinuesa et al. 2005
https://en.wikipedia.org/wiki/Neorhizobium_huautlense

Coale et al. 2024
https://en.wikipedia.org/wiki/Braarudosphaera_bigelowii

Kantor et al. 2024
https://en.wikipedia.org/wiki/Atelocyanobacterium_thalassa
