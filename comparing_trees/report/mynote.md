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

| Step | Threshold | Meaning |
|---|---:|---|
| Sample-gene coverage | percent covered >= 80% | At least 80% of the gene length has mapped read coverage |
| Sample-gene depth | mean depth >= 10 | Average coverage across the gene is at least 10x |
| Mixed-site depth | site depth >= 20 | A position must have enough reads before testing mixed signal |
| Mixed-site minor count | minor allele count >= 5 | The minor base must be supported by at least 5 reads |
| Mixed-site minor fraction | minor allele fraction >= 0.20 | The minor base must represent at least 20% of reads at that site |
| Strict-vs-mixed classification | mixed positions > 10 | Samples with more than 10 mixed positions are excluded from strict tree and retained in mixed-IUPAC tree |

---

## 11. Threshold Interpretation and Literature Support

These exact thresholds are pipeline parameters chosen to be conservative and reproducible. They are not universal biological constants.

However, the general strategy is well supported:

- Consensus calling from mapped reads commonly uses base-quality filtering, minimum depth, and allele-frequency thresholds.
- IUPAC ambiguity codes are commonly used when more than one nucleotide has enough support at a position.
- A 20% allele-frequency threshold is a common consensus/ambiguity threshold in tools such as iVar.
- Minimum depth thresholds such as 10x or 20x are commonly used to avoid making consensus calls from very low read support.
- The final `>10 mixed positions` rule is a conservative sample-level filter used here to separate mostly single-dominant sequences from sequences with stronger mixed-template signal.

Therefore, the thresholds should be reported explicitly as part of the method rather than described as universal standards.

### References Supporting the General Approach

- iVar documentation supports consensus generation from `samtools mpileup`, use of minimum base quality, minimum depth, frequency thresholds, and IUPAC ambiguity codes. It lists `q=20`, minimum depth `10`, and frequency thresholds including `0.2` as commonly used settings:  
  https://andersen-lab.github.io/ivar/html/manualpage.html

- ViralConsensus describes consensus generation from mapped reads using base-quality filtering, minimum depth, and frequency thresholds; its default minimum base quality is 20 and default minimum depth is 10:  
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10212278/

- Threshold-based consensus methods for molecular sequences are a recognized class of consensus methods, where ambiguity codes can be returned depending on nucleotide support thresholds:  
  https://doi.org/10.1016/S0022-5193(05)80692-7

- ANDES discusses threshold-driven consensus generation and polymorphism detection, emphasizing that consensus sequence construction depends on user-defined thresholds for underlying variation:  
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2921379/

- IQ-TREE UFBoot convergence uses a bootstrap correlation coefficient threshold of 0.99 by default, supporting our use of `>=0.99` as the convergence criterion for tree support stability:  
  https://www.iqtree.org/doc/Command-Reference















