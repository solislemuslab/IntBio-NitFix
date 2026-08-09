# Per-gene functional tree pipelines

This directory replaces the earlier combined multigene script with one complete folder per gene.
Each gene folder contains all code needed for that gene only:

1. extract strict and mixed-IUPAC consensus sequences
2. align both FASTA files with MAFFT
3. build strict tree with IQ-TREE nm5000
4. build mixed-IUPAC tree with IQ-TREE nm5000
5. BLAST sample consensus sequences against the corresponding taxon reference alignment
6. check all output files and IQ-TREE warnings

Genes included:

```text
nifH,nifD,nifK,nifJ,nodL,nolG,nolF,noeA,noeB,nodX
```

Copy this whole folder to the cluster:

```bash
LOCAL="/Users/rosa/Desktop/ALLWork/Madison/Project/ryan-nitfix/intbio-nitfix-git/symbiosis_sorted_all_sample_consensus_sequences_50percent/code/gene_tree_pipelines"
REMOTE="raghdam@solislemus-001.discovery.wisc.edu:/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent/Rscripts_v2/"
scp -r "$LOCAL" "$REMOTE"
```

Run one step for all genes in parallel, for example:

```bash
OUT="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent"
GENES="nifH,nifD,nifK,nifJ,nodL,nolG,nolF,noeA,noeB,nodX"

STEP=02 JOBS=2 GENES="$GENES" bash "$OUT/Rscripts_v2/gene_tree_pipelines/run_selected_step_for_genes.sh"
STEP=03 JOBS=2 GENES="$GENES" bash "$OUT/Rscripts_v2/gene_tree_pipelines/run_selected_step_for_genes.sh"
STEP=04 JOBS=2 GENES="$GENES" bash "$OUT/Rscripts_v2/gene_tree_pipelines/run_selected_step_for_genes.sh"
STEP=05 JOBS=2 GENES="$GENES" bash "$OUT/Rscripts_v2/gene_tree_pipelines/run_selected_step_for_genes.sh"
STEP=06 JOBS=2 GENES="$GENES" bash "$OUT/Rscripts_v2/gene_tree_pipelines/run_selected_step_for_genes.sh"
STEP=07 JOBS=2 GENES="$GENES" bash "$OUT/Rscripts_v2/gene_tree_pipelines/run_selected_step_for_genes.sh"
```

For maximum safety, run tree steps in tmux.
