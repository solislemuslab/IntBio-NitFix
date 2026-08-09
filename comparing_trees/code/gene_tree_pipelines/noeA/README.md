# noeA consensus50 tree pipeline

This folder contains the complete per-gene code for `noeA`.

Run order on the cluster:

```bash
OUT="/mnt/dv/wid/projects6/SolisLemus-Intbio-raw/processed-data/symbiosis_sorted_v2_consensus_sequences_50percent"

bash "$OUT/Rscripts_v2/gene_tree_pipelines/noeA/02_run_noeA_consensus_consensus50.sh"
bash "$OUT/Rscripts_v2/gene_tree_pipelines/noeA/03_align_noeA_consensus50_with_mafft.sh"
bash "$OUT/Rscripts_v2/gene_tree_pipelines/noeA/04_build_noeA_strict_tree_nm5000_consensus50.sh"
bash "$OUT/Rscripts_v2/gene_tree_pipelines/noeA/05_build_noeA_iupac_tree_nm5000_consensus50.sh"
bash "$OUT/Rscripts_v2/gene_tree_pipelines/noeA/06_blast_noeA_consensus_to_taxon_reference.sh"
bash "$OUT/Rscripts_v2/gene_tree_pipelines/noeA/07_check_noeA_tree_outputs.sh"
```

Default thresholds:

- sample-gene coverage selection: percent covered >= 80 and mean depth >= 10
- consensus sequence inclusion: N percent <= 20
- pileup filters: base quality >= 20 and mapping quality >= 10
- mixed-site rule: depth >= 20, minor allele count >= 5, minor allele fraction >= 0.20
- sequence flagged mixed when mixed positions > 10
- tree model/search: IQ-TREE `-st DNA -m MFP -B 1000 -alrt 1000 -nm 5000 -T AUTO`

Outputs are written under:

```text
$OUT/gene_trees/noeA/
```
