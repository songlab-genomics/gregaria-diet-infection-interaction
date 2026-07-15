# Code

Save command-line scripts and shared R code here.

## Snakemake RNA-seq pipeline copy

This folder now contains a working copy of the Snakemake RNA-seq pipeline from
`../locust-time-course-RNAseq/code/`. The copy is intended as the starting point
for rerunning the fat body RNA-seq data from raw reads on the cluster.

Copied components:

- `Snakefile` plus specialized Snakefiles for RNA-seq variants,
  metatranscriptomics QC, and eggNOG annotations.
- `rules/` with sample discovery, fastp trimming, STAR/Salmon mapping,
  featureCounts, RNA-seq variant calling, eggNOG, and metatranscriptomics QC.
- SLURM launchers and `cluster.json` from the time-course project.
- Python helper scripts used by the annotation and metatranscriptomics rules.

Before running this for the fat body project, update the copied configuration to
point to the Mehreen raw reads, the intended S. gregaria reference genome/GTF,
the Metarhizium reference resources, and a new run-specific output directory.
