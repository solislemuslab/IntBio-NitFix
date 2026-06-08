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
   one sample + one nif/nod target location
   ex:Sample CLBJ-40-1-No + nifH target ref63 = one row
  |
  v

10. Summarize coverage by sample type
    Groups:
      BLAN = blank control samples, used to detect background/contamination signal - they are not real samples
      No   = nodule samples, the main biological samples for nitrogen-fixing symbionts
      Rh   = rhizosphere samples, soil/root-zone microbial community
      Ro   = root samples, microbes associated with root tissue

    Purpose: compare coverage across controls and biological sample types

   Result: nodule samples showed the strongest support for several nif genes.
        The best-supported genes included nifH, nifU, nifD, nifK, nifE, and nifB.
        For example, nifH, nifU, and nifD had good coverage in more than 96%
        of nodule samples. BLAN controls also showed target-associated signal,
        so blank-aware filtering was needed before choosing final targets.
|
v

11. Rank targets with blank-aware filtering

Why I did this:
I had many possible nif/nod target regions.
Before building a tree, I needed to choose a target that was strongly supported
in real biological samples but not strongly supported in blank controls.

Why I kept BLAN samples in this step:
BLAN samples are controls.
They help show background signal, contamination, or technical carryover.
If I removed BLAN samples too early, I would not know whether a target also
appeared strongly in controls.

What I compared:
  - No samples = nodule samples, where nif/nod signal is expected to be strongest
  - BLAN samples = blank controls, where strong signal is not expected

What I looked for:
  - many nodule samples with good coverage
  - few BLAN samples with good coverage
  - higher read depth in nodules than in BLAN controls

Result:
  21 target regions were classified as PROMISING.

Then I applied a stricter pilot filter:
  - target must be PROMISING
  - at least 75 nodule samples must have good coverage

After this stricter filter:
  - 10 targets remained
  - all 10 were nif genes
  - no nod target passed this strict pilot filter

Meaning:
This step helped me choose the cleanest target for the first pilot tree.
It does not mean nod genes are unimportant.
It means nod genes need a separate nod-focused review.



```



References:

### References supporting Step 11: Blank-aware target ranking

| Pipeline decision | Paper | Exact sentence from the paper | How it helps this step |
|---|---|---|---|
| Focus first on nodule samples (`No`) | Sprent et al. / review on legume-rhizobia specificity, **“Specificity in Legume-Rhizobia Symbioses”** | “Most legume species can fix atmospheric nitrogen (N2) via symbiotic bacteria (general term ‘rhizobia’) in root nodules...” | This supports using nodule samples as the main biological group for choosing the first `nif`/`nod` phylogenetic target, because nodules are the tissue where nitrogen-fixing symbionts are expected to be strongest. |
| Use BLAN controls before selecting targets | Salter et al. 2014, **“Reagent and laboratory contamination can critically impact sequence-based microbiome analyses”** | “Concurrent sequencing of negative control samples is strongly advised.” | This supports keeping BLAN controls in the QC step instead of immediately removing them, because blanks help identify background or contamination signal. |
| Compare biological samples with negative controls | Davis et al. 2018, **“Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data”** | “Sequences from contaminating taxa are likely to have higher prevalence in control samples than in true samples.” | This supports our blank-aware ranking logic: a target is more trustworthy when it is stronger in biological samples than in BLAN controls. |
| Use negative controls to identify problematic sequence features | decontam documentation / Davis and Callahan | “The prevalence ... of each sequence feature in true positive samples is compared to the prevalence in negative controls to identify contaminants.” | This supports comparing target-region support in real samples versus BLAN controls before selecting targets for tree building. |
