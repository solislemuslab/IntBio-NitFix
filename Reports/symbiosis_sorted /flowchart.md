## Why I Used `symbiosis_sorted`

Ryan’s notes explain that the sequencing data were already separated into locus-specific folders. He wrote:

> “It takes a lot of compute to sort these loci out, so I thought I would send along pre-sorted subsets. These are in 16S_sorted, ITS_sorted, and symbiosis_sorted.”

He also explained that the symbiosis data are the functional-gene data for tree building:

> “I have attached a genomic reference in case you are interested in the functional genes in symbiosis_sorted.”

And for phylogenetic analysis, he wrote:

> “I am still finalizing functional gene data, which are subject to a different pipeline. Claudia and I already discussed using this to derive phylogenetic trees...”

Based on these instructions, I used `symbiosis_sorted` for the `nif`/`nod` functional-gene analysis. I did not continue with the 16S/SILVA approach because 16S is mainly for taxonomic/community analysis, while `nif` and `nod` are functional symbiosis genes.

```text
Ryan’s instructions
  |
  |  “pre-sorted subsets... 16S_sorted, ITS_sorted, and symbiosis_sorted”
  |  “functional genes in symbiosis_sorted”
  |  “using this to derive phylogenetic trees”
  v
I used symbiosis_sorted for nif/nod analysis
  |
  v
1. Copy raw reads into processed-data
   Input: raw symbiosis FASTQ files
   Result: 1,116 R1 and 1,116 R2 files copied
  |
  v
2. Trim reads with fastp
   Purpose: clean adapters, low-quality bases, poly-G/poly-X artifacts
   Result: 1,116 paired samples trimmed; 0 failed
  |
  v
3. Check trimming quality
   Purpose: confirm reads are clean enough for mapping
   Result: quality was good; mean read retention 95.84%;
           mean Q30 improved from 95.39% to 97.49%
  |
  v
4. Map reads to Ryan’s symbiosis reference
   Reference: symbiosis_islands.fasta
   Purpose: find where reads align on Ryan’s symbiosis reference
   Result: 1,116 BAM files created successfully
  |
  v
5. Check mapping quality
   Purpose: confirm reads mapped well
   Result: mean mapped reads 93.3%;
           mean properly paired reads 85.0%
  |
  v
6. Extract nif/nod target regions
   Input: symbiosis_islands.gb + gene list
   Purpose: find the exact nif/nod gene regions
   Result: 231 unique central nif/nod regions extracted
           169 nif regions + 62 nod regions
  |
  v
7. Check BLAN controls
   Purpose: test whether blank controls contain target signal
   Result: BLAN controls had target-associated signal,
           so they are used for QC but excluded from final trees
  |
  v
8. Match nif/nod targets back to Ryan’s reference
   Purpose: connect extracted nif/nod genes to original mapping coordinates
   Result: all 231 targets matched;
           233 exact reference locations found
  |
  v
9. Measure coverage for each gene in each sample
   Purpose: ask which samples have enough reads covering each nif/nod gene
   Result: 260,028 sample-region coverage rows generated [1,116 samples × 233 nif/nod target locations = 260,028 rows]
   So each row means:
   one sample + one nif/nod target locatio
   ex:Sample CLBJ-40-1-No + nifH target ref63 = one row
  |
  v
```
