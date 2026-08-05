# Code

For the independent audit of DEG-rich unplaced scaffolds, see
[`SCAFFOLD_ORIGIN_AUDIT.md`](SCAFFOLD_ORIGIN_AUDIT.md). That target starts from
the frozen 77-scaffold candidate set and cannot rerun RNA-seq preprocessing.

For the focused BLASTn/DIAMOND follow-up of the 948 unresolved genes in the
infected diet 33 versus diet 50 contrast, see
[`SCAFFOLD_ORIGIN_HOMOLOGY.md`](SCAFFOLD_ORIGIN_HOMOLOGY.md). It is a separate
run and cannot resubmit RNA-seq or the first scaffold audit.

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

## Host-pathogen manuscript rerun

The additive dual RNA-seq entry point is:

```text
Snakefile.host_pathogen_dual_rnaseq
```

It maps every freshly trimmed library twice: first to *S. gregaria* alone and
then competitively to *S. gregaria* plus the *M. robertsii* proxy reference.
The default HPC target stops after fastp, STAR, featureCounts, mapping
comparisons, and Kraken/Bracken. Host and fungal DESeq2 analyses are run locally
after transferring the compact count and taxonomy tables back from Sol. See
[`HOST_PATHOGEN_DUAL_RNASEQ.md`](HOST_PATHOGEN_DUAL_RNASEQ.md) for exact HPC
paths, fungal-reference download commands, preflight checks, and launch steps.
