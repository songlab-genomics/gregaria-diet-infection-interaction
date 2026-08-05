# Unplaced-scaffold origin audit

## Purpose

This independent workflow tests the 77 unplaced scaffolds carrying the 986
unified-model fat-body DEG genes. It does not rerun fastp, STAR,
featureCounts, Kraken on RNA-seq reads, or DESeq2.

The candidate definition is frozen in `config/scaffold_origin/` with:

- the 77 scaffold accessions;
- the 986 unique candidate DEG genes and their contrast memberships;
- the source run/model and NCBI assembly placement;
- the official NCBI FCS report used for the reference assembly.

The old 80-scaffold list remains untouched but is not used by this target.

## What the workflow tests

1. Extract the complete DNA sequence of each candidate scaffold from the
   *S. gregaria* RefSeq FASTA.
2. Calculate sequence length, GC content, ambiguous-base content, all annotated
   genes, candidate DEG genes, and gene density.
3. Align each complete scaffold to the *S. piceifrons*, *S. americana*, and
   *S. serialis cubense* genome assemblies to measure cross-locust conservation.
4. Classify each complete scaffold with the installed Kraken PlusPF database
   and reconstruct its NCBI taxonomic lineage.
5. Combine cross-species alignment, Kraken, NCBI assembly, annotation, and FCS
   evidence in one table.
6. Retain every scaffold pending independent validation. Kraken alone never
   triggers automatic removal.

The official NCBI FCS report did not flag these scaffolds, but its screen used a
2023 database. That is useful evidence, not proof of locust origin. The next
strong tests are nucleotide/protein homology, conservation in other
*Schistocerca* genomes, and genomic DNA coverage/Hi-C support. Parallel qPCR
would have been useful, but no matching tissue remains.

## Start a new HPC run

Sync the updated project folder to:

```text
/scratch/mtecher/gregaria-diet-infection-interaction
```

Then:

```bash
cd /scratch/mtecher/gregaria-diet-infection-interaction/code
module purge
module load mamba
source activate myENV

export ORIGIN_RUN_ID=scaffold_origin_$(date +%Y%m%d_%H%M%S)
unset ORIGIN_DIR
export ORIGIN_RESUME_EXISTING=false

bash pilot/check_scaffold_origin_audit_hpc.sh
snakemake -s Snakefile.scaffold_origin_audit \
  scaffold_origin_audit_hpc --dry-run --rerun-incomplete
sbatch snakemake.scaffold_origin_audit.slurm
```

If preflight reports that minimap2 is missing, submit
`scripts/install_metatranscriptomics_qc_env.slurm`, wait for that installer to
finish, and rerun preflight. It updates the existing Kraken/Bracken environment
in place without deleting it.

The dry-run must list only rules beginning with `origin_` plus the aggregate
target. If it lists fastp, STAR, or featureCounts, stop and check the Snakefile
name.

## Check progress

```bash
bash pilot/check_scaffold_origin_audit_progress.sh "$ORIGIN_RUN_ID"
```

This command is read-only.

## Main outputs

```text
output/runs/<ORIGIN_RUN_ID>/
├── 00-input/       frozen inputs and SHA-256 provenance
├── 01-sequences/   candidate scaffold FASTA and sequence metrics
├── 02-annotation/  all genes and gene-density summary
├── 03-taxonomy/    Kraken calls and reconstructed lineage
├── 04-cross-species/ minimap2 PAF files and conservation summary
├── 05-evidence/    integrated evidence and follow-up table
└── 06-local-transfer/scaffold_origin_audit_for_local.tar.gz
```

The main interpretation table is:

```text
05-evidence/scaffold_origin_evidence.tsv
```

Its `automatic_filter_decision` is deliberately
`retain_pending_independent_validation` for every scaffold.
