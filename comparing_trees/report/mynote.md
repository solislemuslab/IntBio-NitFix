Email 11 July:

I checked the sample IDs in the metadata file (intbio_metadata_draft3.csv) against the sample IDs used in analysis.

The corrected metadata file contains 2,898 unique sample IDs, while the  analysis contains 2,907 unique samples. After harmonizing the OAES naming format, the final comparison gave 24 samples present in the analysis but not in the metadata, and 15 samples present in the metadata but not in the analysis. 

One formatting issue I corrected was for the OAES_19 and OAES_24 samples. In the raw sequencing files and analysis outputs, these samples use an underscore in the site label, for example OAES_19-2-Rh and OAES_24-3-No.  I renamed these in the metadata to match the sequencing/analysis IDs, so they are no longer counted as missing.

The six renamed metadata IDs were:

OAES-19-2-No  ->  OAES_19-2-No
OAES-19-2-Rh  ->  OAES_19-2-Rh
OAES-19-2-Ro  ->  OAES_19-2-Ro
OAES-24-3-No  ->  OAES_24-3-No
OAES-24-3-Rh  ->  OAES_24-3-Rh
OAES-24-3-Ro  ->  OAES_24-3-Ro
After this correction, the samples present in the analysis but not in the metadata are:

```text
N13-2-D-No
N13-2-D-Rh
N13-2-D-Ro
N34-2-D-No
N34-2-D-Rh
N34-2-D-Ro
N35-4-D-No
N35-4-D-Rh
N35-4-D-Ro
N36-4-D-No
N36-4-D-Rh
N36-4-D-Ro
N38-2-D-No
N38-2-D-Rh
N38-2-D-Ro
N39-2-D-No
N39-2-D-Rh
N39-2-D-Ro
N49-3-D-No
N49-3-D-Rh
N49-3-D-Ro
UKFS-52-1-No
UKFS-52-1-Rh
UKFS-52-1-Ro
```

The samples present in the metadata but not in the analysis are:

```text
CLBJ-24-Ro
CLBJ-29-2-No
JERC-4-1-Ro
JERC-6-3-Ro
JERC-8-2-Ro
JERC-10-3-Ro
JERC-RHRE-18-2-No
OAES-5-No
TALL-2-1-No
TALL-2-1-Rh
TALL-2-1-Ro
TALL-91-2-No
UKFS-2-1-No
UKFS-2-1-Rh
UKFS-2-1-Ro
```


##########################################################################################

The per-gene taxon reference alignments were checked against the gene-list table. 
The only gene listed in symbiosis_islands_gene_list.xlsx that was not present in taxon_fastas_by_gene was nopT. The selected 10 tree genes all had matching taxon-reference alignments.

###################################################################################

# Strict and Mixed-IUPAC Consensus Calling

This note explains how sample-gene pairs were selected, how consensus sequences were constructed, why strict and mixed-IUPAC trees have different numbers of samples, and what thresholds were used.

This part is easy to confuse because there are two different decisions:

1. **Should this sample-gene pair be used at all?**
2. **If yes, what nucleotide should we write at each position?**

These two steps use different rules.

---

## 1. Sample-Gene Filtering

First, for each sample and each gene, reads were mapped to the consensus50 reference gene sequence. Then we measured whether enough of that gene was covered by mapped reads.

A sample-gene pair was kept only if it passed both filters:

```text
percent covered >= 80%
mean depth >= 10
```

### Percent Covered

Percent covered means the fraction of the gene length covered by mapped reads.

For example, if `nifH` is 997 bp long and reads cover 850 positions:

```text
percent covered = 850 / 997 = 85.3%
```

This sample passes the coverage filter because:

```text
85.3% >= 80%
```

This does **not** mean that 80% of reads match exactly. It means that at least 80% of gene positions have mapped read coverage.

### Mean Depth

Mean depth means the average number of reads covering positions across the gene.

For example:

```text
mean depth = 35x
```

This passes because:

```text
35 >= 10
```

So after this first step, we keep only sample-gene pairs with enough coverage and depth.

---

## 2. Position-by-Position Consensus Calling

After a sample-gene pair passes coverage/depth filtering, we examine the gene position by position.

At each position, we count the number of high-quality reads supporting:

```text
A, C, G, T
```

Then we decide whether the site is:

- clearly single-dominant
- missing or uncertain
- mixed

---

## 3. Mixed-Site Rule

A position was called mixed only when all of these were true:

```text
site depth >= 20
minor allele count >= 5
minor allele fraction >= 0.20
```

This means we do not call a site mixed just because one or two reads disagree with the dominant base. The minor base must have enough read support.

For example, if a site has 100 reads, a minor allele must be supported by at least 5 reads and represent at least 20% of the reads to be considered mixed.

---

## 4. Example of One Mixed Position

Suppose one position has 100 reads:

| Base | Read count |
|---|---:|
| A | 80 |
| G | 20 |
| C | 0 |
| T | 0 |

Total depth:

```text
100
```

Minor allele:

```text
G
```

Minor allele count:

```text
20
```

Minor allele fraction:

```text
20 / 100 = 0.20
```

This position passes the mixed-site rule:

```text
depth 100 >= 20
minor allele count 20 >= 5
minor allele fraction 0.20 >= 0.20
```

So this site is considered mixed.

In the **strict consensus**, we write the dominant base:

```text
A
```

In the **mixed-IUPAC consensus**, we write the ambiguity code:

```text
R
```

because:

```text
R = A/G
```

So for the same position:

| Consensus type | Output base |
|---|---|
| Strict consensus | A |
| Mixed-IUPAC consensus | R |

---

## 5. Why Strict Consensus Does Not Contain IUPAC Letters

The strict consensus FASTA contains only:

```text
A, C, G, T, N
```

It does **not** contain IUPAC ambiguity letters such as `R`, `Y`, `M`, or `K`.

Even if a position has mixed evidence, the strict consensus writes the dominant base.

For example:

| Base | Read count |
|---|---:|
| A | 80 |
| G | 20 |

Strict output:

```text
A
```

Mixed-IUPAC output:

```text
R
```

So the strict sequence itself does not show the mixed letter. The mixed-site count is calculated from the read pileup before writing the final strict FASTA.

---

## 6. Why We Used the “More Than 10 Mixed Positions” Rule

After checking every position in a gene, we count how many positions are mixed.

A sample is considered **single-dominant** when:

```text
mixed positions <= 10
```

A sample is considered **mixed possible multitemplate** when:

```text
mixed positions > 10
```

This rule is applied after consensus calling evidence is evaluated across the gene. In other words, if the consensus-building step finds more than 10 positions that pass the mixed-site criteria:

```text
site depth >= 20
minor allele count >= 5
minor allele fraction >= 0.20
```

then we classify that sample-gene sequence as mixed instead of strict.

A few mixed-looking positions can happen because of sequencing noise, mapping uncertainty, local alignment issues, or low-frequency variation. Therefore, we allowed up to 10 mixed positions and still treated the sample as mostly single-dominant.

But if a sample has more than 10 mixed positions, that suggests stronger mixed-template, multicopy, or multiple-organism signal. In that case, we do not trust the sample as one clean dominant sequence for the strict tree.

So the classification rule is:

```text
mixed positions <= 10  -> pass_single_dominant
mixed positions > 10   -> pass_mixed_possible_multitemplate
```

Important: the `>10 mixed positions` rule does not mean that the strict sequence has more than 10 IUPAC letters. The strict sequence never writes IUPAC letters. The count comes from the read evidence at each position before the strict FASTA is written.

---

## 7. Why Strict and Mixed-IUPAC Trees Have Different Numbers of Samples

The strict tree includes only:

```text
pass_single_dominant
```

The mixed-IUPAC tree includes:

```text
pass_single_dominant + pass_mixed_possible_multitemplate
```

So the mixed-IUPAC tree is larger because it keeps samples that passed coverage/depth filters but had more than 10 mixed positions.

Example for `nifH`:

| Category | Count |
|---|---:|
| fail_high_N | 34 |
| pass_single_dominant | 387 |
| pass_mixed_possible_multitemplate | 1,151 |

Therefore:

```text
Strict nifH tree = 387 sequences
Mixed-IUPAC nifH tree = 387 + 1,151 = 1,538 sequences
```

The strict tree is smaller because it excludes samples with strong mixed signal.

---

## 8. Example: Sample Included in Mixed-IUPAC Tree but Not Strict Tree

Suppose one sample has good `nifH` coverage:

```text
percent covered = 92%
mean depth = 35x
```

So it passes the first filter.

Now suppose we inspect all 997 positions of `nifH`, and 25 positions pass the mixed-site rule.

```text
mixed positions = 25
```

Because:

```text
25 > 10
```

this sample is not included in the strict tree.

But it is included in the mixed-IUPAC tree.

| Tree type | Included? | Reason |
|---|---|---|
| Strict tree | No | More than 10 mixed positions |
| Mixed-IUPAC tree | Yes | Passed coverage/depth; mixed positions represented using IUPAC codes |

In the mixed-IUPAC sequence, those 25 mixed positions would appear as IUPAC ambiguity letters such as:

```text
R, Y, M, K, W, S, B, D, H, V
```

---

## 9. IUPAC Codes Used

| Supported bases | IUPAC code |
|---|---|
| A + G | R |
| C + T | Y |
| G + C | S |
| A + T | W |
| G + T | K |
| A + C | M |
| C + G + T | B |
| A + G + T | D |
| A + C + T | H |
| A + C + G | V |
| A + C + G + T | N |

Note: In this pipeline, `N` is also used when a position is missing, uncertain, or has insufficient reliable support.

---

## 10. Threshold Summary

The thresholds used here are conservative pipeline choices, but they follow the same logic used in existing consensus-calling tools.

| Threshold in our pipeline | Support |
|---|---|
| Base quality `>=20` | iVar uses default base quality `20`; ViralConsensus also uses default base quality `20`. |
| Mean depth `>=10` | iVar uses default minimum depth `10`; ViralConsensus also uses default minimum depth `10`. |
| Minor allele fraction `>=0.20` | iVar lists `0.20` as a commonly used frequency threshold for bases supported by at least 20% of reads. |
| IUPAC ambiguity codes | Supported by iVar, Day and McMorris 1992, and ANDES. ViralConsensus does not output IUPAC codes. |
| UFBoot convergence `>=0.99` | IQ-TREE uses bootstrap correlation coefficient `0.99` as the default UFBoot convergence criterion. |

Other thresholds, including percent covered `>=80%`, site depth `>=20`, minor allele count `>=5`, and `>10` mixed positions, are pipeline-specific conservative filters and should be reported explicitly.

### References

- iVar manual: https://andersen-lab.github.io/ivar/html/manualpage.html  
- ViralConsensus: https://pmc.ncbi.nlm.nih.gov/articles/PMC10212278/  
- Day and McMorris 1992: https://doi.org/10.1016/S0022-5193(05)80692-7  
- ANDES: https://pmc.ncbi.nlm.nih.gov/articles/PMC2921379/  
- IQ-TREE command reference: https://www.iqtree.org/doc/Command-Reference



## All steps (my note)

## Per-Gene Tree Pipeline Explanation

Before this per-gene pipeline, all trimmed reads were already mapped to the consensus50 functional-gene reference. That mapping step produced one BAM file per sample. The per-gene pipeline starts from those BAM files.

Although each gene folder contains a file named `01_extract_gene_consensus_consensus50.py`, we do not run `STEP=01` directly with the wrapper script. That Python file contains the consensus-calling logic and is called inside `STEP=02`.

### Step 02. Extract Per-Gene Consensus Sequences

For each selected gene, the pipeline looks inside each sample BAM file and extracts the reads that mapped to that gene. Then it uses the pileup to summarize the mapped reads into one consensus sequence per sample-gene pair.

A pileup means that, for every position in the gene, the code checks which bases were observed in the mapped reads.

Example:

```text
Gene position 25:
Reads covering this position: A, A, A, A, G
Consensus base: A
```
So instead of keeping thousands of reads for one sample, we create one representative DNA sequence for that sample and gene.
The important sentence:
“the consensus-calling logic was designed to avoid copying reference ambiguity into sample sequences”
means this:
The consensus50 reference contains some ambiguous bases, such as N, R, Y, S, etc. These ambiguous bases come from the reference-building process. They represent uncertainty or variation in the reference, not necessarily in our sample. Therefore, when we build a sample consensus sequence, we do not want to automatically copy those ambiguous letters from the reference into the sample.
Instead, the sample consensus is based on the actual reads from that sample.
```text
Example 1:
Reference base: R
Meaning of R: A or G
Reads from sample: A, A, A, A, A
Sample consensus: A
```
We call A, not R, because the sample reads support A.
Example 2:
```text
Reference base: N
Meaning of N: unknown/ambiguous reference position
Reads from sample: C, C, C, C, C
Sample consensus: C
```
We call C, not N, because the sample reads give clear evidence.
Example 3:
```text
Reference base: R
Reads from sample: A, A, A, G, G
Sample consensus in strict file: A
Sample consensus in mixed-IUPAC file: R
```
Here, the sample itself has evidence for both A and G, so the mixed-IUPAC sequence records that mixed signal as R.
This is different from copying R from the reference. The R in the sample sequence is written only because the sample reads support both bases.
Thresholds Used In Consensus Calling
A sample-gene pair was considered for consensus calling when:
```text
percent covered >= 80%
mean depth >= 10
For each position, mixed signal was counted when:
site depth >= 20
minor allele count >= 5
minor allele fraction >= 0.20
```
Then the sample was classified as:
pass_single_dominant
if it had 10 or fewer mixed positions.
pass_mixed_possible_multitemplate
if it had more than 10 mixed positions.
The strict FASTA contains only samples classified as pass_single_dominant.
The mixed-IUPAC FASTA contains both:
pass_single_dominant
pass_mixed_possible_multitemplate
That is why the mixed-IUPAC tree has more tips than the strict tree.
What Each Pipeline Step Does
STEP 02: Consensus Extraction
This step creates sample-specific consensus sequences for each gene.
It calls the Python script:
01_extract_gene_consensus_consensus50.py
inside each gene folder.
Outputs include:
```text
<gene>_consensus_qc.tsv
<gene>_consensus50_strict_single_dominant_pct80_depth10_Nle20.fasta
<gene>_consensus50_iupac_all_pass_pct80_depth10_Nle20.fasta
```
STEP 03: MAFFT Alignment
This step aligns the consensus sequences so that homologous positions are compared correctly.
Inputs:
strict FASTA
mixed-IUPAC FASTA
Outputs:
```text
<gene>_consensus50_strict_single_dominant.mafft.fasta
<gene>_consensus50_iupac_all_pass.mafft.fasta
STEP 04: Strict Tree Construction
```
This step builds the strict single-dominant tree with IQ-TREE.
Input:
<gene>_consensus50_strict_single_dominant.mafft.fasta
This tree includes only samples with clean single-dominant consensus sequences.
```text
IQ-TREE uses:
-st DNA
-m MFP
-B 1000
-alrt 1000
-nm 5000
-T AUTO
```
Output:
```text
<gene>_consensus50_strict_single_dominant_nm5000.treefile
<gene>_consensus50_strict_single_dominant_nm5000.contree
<gene>_consensus50_strict_single_dominant_nm5000.iqtree
<gene>_consensus50_strict_single_dominant_nm5000.log
```
STEP 05: Mixed-IUPAC Tree Construction
This step builds the mixed-IUPAC tree with IQ-TREE.
Input:
<gene>_consensus50_iupac_all_pass.mafft.fasta
This tree includes:
pass_single_dominant samples
pass_mixed_possible_multitemplate samples
It is useful as a broader exploratory/sensitivity tree, but it should be interpreted carefully because mixed-IUPAC calls may reflect multiple templates in the same sample.
Outputs:
```text
<gene>_consensus50_iupac_all_pass_nm5000.treefile
<gene>_consensus50_iupac_all_pass_nm5000.contree
<gene>_consensus50_iupac_all_pass_nm5000.iqtree
<gene>_consensus50_iupac_all_pass_nm5000.log
```
STEP 06: BLAST/Taxon Assignment
This step compares each sample consensus sequence to known reference sequences for the same gene.
It asks:
Which known reference sequence is the closest match to this sample consensus sequence?
Outputs include:
<gene>_best_reference_hit_with_taxon.tsv
This table contains:
```text
sample ID
best reference ID
closest BLAST taxon
percent identity
query coverage
e-value
bitscore
```
Important interpretation:
The BLAST label is the closest known reference hit. It is not proof of exact species identity.
STEP 07: Check Tree Outputs
```text
This step checks and summarizes the tree results.
It reports things such as:
number of tree tips
alignment length
best-fit model
tree length
bootstrap convergence
warnings
```
This is the summary/QC step used for the report.



Mapping reads to consensus50 reference -> BAM files

Then per-gene pipeline:
STEP=02 -> runs 02 shell script -> internally calls 01 Python script
STEP=03 -> MAFFT alignment
STEP=04 -> strict tree
STEP=05 -> mixed-IUPAC tree
STEP=06 -> BLAST/taxon assignment
STEP=07 -> check outputs

```text
Mapped BAM files
   |
   | STEP 02
   v
Per-gene consensus sequences
   |
   | STEP 03
   v
MAFFT alignments
   |
   | STEP 04 and STEP 05
   v
Strict tree and mixed-IUPAC tree
   |
   | STEP 06
   v
BLAST/taxon annotation tables
   |
   | STEP 07
   v
Tree summary/QC tables for report
```













