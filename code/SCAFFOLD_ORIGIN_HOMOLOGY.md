# Focused BLASTn and DIAMOND scaffold-origin audit

## Purpose

This workflow tests the 948 unresolved DEG genes on 59 unplaced scaffolds that
contribute to the infected diet 33 versus diet 50 contrast.

It is independent of the completed RNA-seq and first scaffold-audit runs. It
does not invoke fastp, STAR, featureCounts, DESeq2, Kraken, or minimap2.

The target set is frozen in `config/scaffold_homology/`:

- `infected_diet33_vs50_unresolved_scaffolds_20260731.txt`
- `infected_diet33_vs50_unresolved_genes_20260731.tsv`
- `infected_diet33_vs50_unresolved_provenance_20260731.tsv`

## Searches

1. Extract one strand-oriented genomic DNA query for each of the 948 genes.
2. Extract every RefSeq protein isoform linked to those genes in the GFF.
3. Search gene DNA with BLASTn against NCBI `core_nt`.
4. Search proteins with DIAMOND against NCBI ClusteredNR, whose standalone
   database basename is `nr_cluster_seq`.
5. Exclude *S. gregaria* taxid 7010 in both search commands and again during
   summarization. The post-search filter is intentional: some database records
   can retain a self taxid despite the search-level exclusion.
6. After removing self-hits, retain hits within 5% of each query's best external
   bit score and summarize their
   NCBI taxonomy by query, gene, and scaffold.

NCBI describes `core_nt` as a broad, non-redundant nucleotide database intended
for sequence characterization. ClusteredNR preserves broad `nr` taxonomic
coverage while using representatives of highly similar protein clusters. NCBI
publishes it separately from the standard BLAST database collection under the
local basename `nr_cluster_seq`; `clustered_nr` is not a valid local basename.

References:

- <https://ncbiinsights.ncbi.nlm.nih.gov/2023/09/26/blast-clusterednr-database-available-download/>
- <https://ftp.ncbi.nlm.nih.gov/blast/db/experimental/>
- <https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options>

The output categories are evidence labels, not automatic filters:

- locust/orthopteran support;
- broader arthropod support;
- other metazoan support;
- mixed/ambiguous;
- non-metazoan candidate;
- no informative hit.

## One-time HPC resources

After syncing the project folder, install the dedicated environment:

```bash
cd /scratch/mtecher/gregaria-diet-infection-interaction/code
sbatch scripts/install_scaffold_homology_env.slurm
```

Download the RefSeq protein FASTA. The script resumes an interrupted `wget` and
leaves an existing valid FASTA unchanged:

```bash
bash scripts/download_gregaria_protein_fasta_hpc.sh
```

Check scratch capacity and download the NCBI databases in a long `public`
partition job:

```bash
df -h /scratch/mtecher
sbatch scripts/download_scaffold_homology_databases.slurm
```

The current ClusteredNR release is a very large 160-volume database. NCBI
specifies at least 128 GB RAM and substantial scratch space for the standalone
database. This workflow keeps downloaded archives as well as extracted database
files, so check that several hundred GB remain before submission.

Monitor the database job:

```bash
squeue --me
tail -f slurm-homology-db-<JOBID>.out
bash pilot/check_scaffold_homology_databases_hpc.sh
```

The download script preserves and resumes completed files. It downloads
`core_nt` with `update_blastdb.pl`, obtains the experimental `nr_cluster_seq`
volume list from NCBI metadata, validates every archive checksum, and writes
`download_complete.tsv` only after `blastdbcmd` can read both databases. It does
not remove downloaded archives.

## Start a new search run

```bash
cd /scratch/mtecher/gregaria-diet-infection-interaction/code
module purge
module load mamba
source activate myENV

export HOMOLOGY_RUN_ID=scaffold_homology_$(date +%Y%m%d_%H%M%S)
unset HOMOLOGY_DIR
export HOMOLOGY_RESUME_EXISTING=false

bash pilot/check_scaffold_origin_homology_hpc.sh

snakemake -s Snakefile.scaffold_origin_homology \
  scaffold_origin_homology_hpc --dry-run --rerun-incomplete

sbatch --export=ALL snakemake.scaffold_origin_homology.slurm
```

The dry-run must list only `homology_` rules and the aggregate
`scaffold_origin_homology_hpc` target.

## Check progress

```bash
bash pilot/check_scaffold_origin_homology_progress.sh "$HOMOLOGY_RUN_ID"
```

This status command is read-only.

To resume the exact same run after confirming no controller or workers remain:

```bash
export HOMOLOGY_RUN_ID=<exact-existing-run-id>
export HOMOLOGY_RESUME_EXISTING=true
sbatch --export=ALL snakemake.scaffold_origin_homology.slurm
```

## Output structure

```text
output/runs/<HOMOLOGY_RUN_ID>/
├── 00-input/                  frozen target set and checksums
├── 01-sequences/              scaffold, gene-DNA, and protein queries
├── 02-blastn-core-nt/         nucleotide hits excluding S. gregaria
├── 03-diamond-nr-cluster-seq/ protein hits excluding S. gregaria
├── 04-evidence/               query, gene, and scaffold origin summaries
└── 05-local-transfer/         checksummed local-analysis archive
```

The main result is:

```text
04-evidence/gene_homology_summary.tsv
```

`04-evidence/non_metazoan_candidate_genes.tsv` is a review list. Its genes must
not be removed solely on the basis of a best-hit label; inspect alignment
coverage, identity, competing top hits, expression, and scaffold-level
consistency first.
