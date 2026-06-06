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
  -> fastp trimming with Ryan's settings
  -> map reads to Ryan's symbiosis_islands.fasta
  -> use symbiosis_islands.gb to identify nif/nod targets
  -> calculate gene-region coverage from mapped BAMs
  -> compare biological samples with BLAN controls
  -> extract pilot nifH consensus sequences
  -> filter low-quality sequences
  -> MAFFT alignment
  -> IQ-TREE pilot tree
  -> iTOL visualization
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

## 3. Step-by-Step Pipeline

### Step 1. Read Quality Control and Trimming

**Purpose:** Clean raw reads before mapping by removing adapters, low-quality bases, poly-G/poly-X artifacts, and short reads.

**Ryan connection:** Ryan recommended `fastp` and provided the trimming settings because NovaSeq data can contain problematic poly-G read-through artifacts.

**Input:** `$BASE/raw_symbiosis_full/*_R1.fq.gz`, `$BASE/raw_symbiosis_full/*_R2.fq.gz`

**Output:** `$BASE/symbiosis_trimmed_ryan_full/`, `$BASE/fastp_logs_ryan_full/`, `$BASE/symbiosis_trimmed_ryan_full_manifest.tsv`

**Script/command used:**

```bash
bash $BASE/Rscripts/trim_symbiosis_full_ryan.sh
```

**What the script does:** Runs `fastp` on each paired sample using Ryan's settings and writes trimmed paired/unpaired reads plus QC logs.

**Samples:** 1,116 paired samples.

**Runtime:** Fast; completed in the same work session.

**Result:** Trimming completed successfully: 1,116 P1 files, 1,116 P2 files, 1,116 JSON reports, and 0 failed samples. BLAN controls were retained for QC.

### Step 2. Mapping Reads to Ryan's Symbiosis-Island Reference

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

### Step 3. Extract Central `nif`/`nod` Targets

**Purpose:** Identify which central `nif` and `nod` genes are present in Ryan's GenBank annotation.

**Ryan connection:** Ryan supplied `symbiosis_islands.gb` and the gene list, so the annotation file was used to identify available functional-gene targets.

**Input:** `$BASE/symbiosis_islands.gb`, `$BASE/symbiosis_islands_gene_list.xlsx`

**Output:** `$BASE/nif_nod_target_reference/`

**Script/command used:**

```bash
python3 $BASE/Rscripts/extract_all_nif_nod_targets.py
```

**What the script does:** Parses the GenBank file, extracts annotated central `nif`/`nod` CDS sequences, removes exact duplicates, and writes FASTA/summary tables.

**Runtime:** Less than 1 minute.

**Result:** 233 annotated CDS records were extracted and 231 unique sequences retained. Of 35 genes listed in the summary sheet, 29 were found. Missing genes were `nifM`, `nifY`, `nodE`, `nodF`, `nodP`, and `nodT`.

### Step 4. Blank-Control Evaluation

**Purpose:** Check whether BLAN controls contain target-associated signal.

**Ryan connection:** BLAN samples are not biological samples, but controls are needed before interpreting functional-gene signal.

**Input:** BLAN reads from `$BASE/symbiosis_trimmed_ryan_full/`, target FASTA from `$BASE/nif_nod_target_reference/all_central_nif_nod_targets.fasta`

**Output:** `$BASE/nif_nod_blank_mapping/`, `$BASE/nif_nod_blank_mapping_logs/`, `$BASE/nif_nod_blank_mapping_summary.tsv`

**Script/command used:**

```bash
bash $BASE/Rscripts/map_all_blank_nif_nod_targets.sh
```

**What the script does:** Maps all 12 BLAN controls to the extracted `nif`/`nod` target FASTA and summarizes mapped reads/top targets.

**Samples:** 12 BLAN controls.

**Runtime:** Fast; completed in the same session.

**Result:** All BLAN controls showed strong target mapping, ranging from 47.35% to 92.08% mapped reads. This means BLAN profiles must be used for QC, and BLAN samples must be excluded from final trees.

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

### Step 6. Gene-Level Coverage From Existing BAMs

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

**Runtime:** Completed in the same session; exact time not recorded.

**Result:** Final table had 260,029 lines: 1 header plus 260,028 sample-region rows.

### Step 7. Coverage Summary by Sample Type

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
