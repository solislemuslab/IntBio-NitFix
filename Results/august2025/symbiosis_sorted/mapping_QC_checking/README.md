This figure summarizes mapping quality across all 1,116 samples after mapping the trimmed reads to Ryan’s symbiosis_islands.fasta reference.

Panel A shows the distribution of mapped-read percentages. This metric asks whether reads mapped anywhere to the reference. Most samples have high mapping rates, clustered around 90–100%, with a mean of 93.3% and median of 93.6%. This means that, for most samples, the majority of trimmed reads successfully aligned to one or more sequences in the symbiosis-island reference.

Panel B shows the distribution of properly paired-read percentages. This is a stricter paired-end mapping metric. Because these data are paired-end, each DNA fragment was sequenced from two directions: Read 1 and Read 2. A read pair is counted as “properly paired” only when both reads from the same DNA fragment map to the reference in the expected orientation and distance. Therefore, Panel B does not only ask whether reads mapped, but whether the paired-end reads mapped together in a biologically and technically consistent way. The mean properly paired rate was 85.0%, with a median of 85.8%. This value is expected to be lower than the overall mapped-read percentage because some reads can map individually while their mate is missing, clipped, mapped to a different reference region, or not positioned as an expected pair.

Overall, the figure shows that mapping was successful across the dataset. Most samples had high mapped-read percentages and strong properly paired-read percentages, indicating that Ryan’s symbiosis_islands.fasta reference was appropriate for these cleaned symbiosis_sorted reads. The 12 BLAN controls are included in the summary because they were retained for QC evaluation.

## my example

Imagine one sample has 100 read pairs:

100 R1 reads
100 R2 reads

So there are 200 total reads from this sample.

The reference is not reads; it is a set of known DNA sequences from `symbiosis_islands.fasta`.

During mapping, the computer asks:

“Where does each sample read best match on the reference DNA?”

If 140 out of 200 reads match somewhere on the reference, then:

140 / 200 = 70%

So the mapped-read percentage is 70%.

**For properly paired:**
Because the data are paired-end, each DNA fragment has two reads: R1 and R2. A read pair is “properly paired” when both R1 and R2 map to the same reference region in the expected direction and distance from each other.

For example, if 150 out of 200 reads are part of properly paired alignments, then:

150 / 200 = 75%

So the properly paired percentage is 75%.

**summary**
Mapping means matching sample reads to known reference DNA sequences. If a sample has 200 reads total and 140 reads match the reference, the mapped-read percentage is 70%. Properly paired means the R1 and R2 reads from the same DNA fragment both map to the reference in the expected orientation and distance. This tells us that the paired-end reads are mapping consistently, not randomly.
