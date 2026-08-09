# Fat-body host-pathogen dual RNA-seq

## Metadata correction and missing library 1044

The original RNA-seq datasheet contains 45 libraries. Sample `1024` is an
infected diet-33 library, not diet 83, and sample `1044` is an infected diet-83
library that was absent from the first 44-library run. The corrected inputs are:

- `../data/metadata/mehreen_metadata_corrected_45.txt` for local R analyses;
- `config/fatbody_dual_rnaseq_samples_corrected_45.tsv` for the final design;
- `config/fatbody_dual_rnaseq_sample_1044.tsv` for the one-library HPC add-on.

The add-on run writes to a new run directory and reuses only the immutable
reference products and STAR indexes from `host_pathogen_dual_20260727_000952`.
It does not modify the completed 44-library run. If the uploaded delivery still
contains separate L1/L2 files, merge it first with:

```bash
export SOURCE_DIR=/path/to/uploaded/mehreen_1044
sbatch --export=ALL pilot/merge_fatbody_1044_lanes.slurm
```

After the two protected `mehreen_1044_MERGE_*.fq.gz` files exist in the raw-read
directory, launch only sample 1044:

```bash
module purge
module load mamba
source activate myENV
bash pilot/submit_fatbody_1044_addon_hpc.sh
```

After its transfer bundle is copied locally, combine its count column with the
frozen 44-sample matrices in a new 45-sample run folder. Then rerun DESeq2 and
all report pages using the corrected metadata; do not reinterpret the existing
44-sample pages as corrected results.

This workflow reruns the fat-body libraries from the original FASTQs and maps
the same freshly trimmed reads in two independent ways:

1. host-only STAR mapping to *S. gregaria*;
2. competitive STAR mapping to *S. gregaria* plus the *M. robertsii* proxy.

It does not reuse the colleague-generated BAMs or count matrices. Every launch
writes to a new timestamped folder.

## HPC locations

The original paired reads are located at:

```text
/data/songlab/sequencing_data/RNAseq/mehreen
```

For example:

```text
/data/songlab/sequencing_data/RNAseq/mehreen/mehreen_1002_MERGE_1.fq.gz
/data/songlab/sequencing_data/RNAseq/mehreen/mehreen_1002_MERGE_2.fq.gz
```

The HPC repository checkout and new results are expected under:

```text
/scratch/mtecher/gregaria-diet-infection-interaction
```

The preferred locust reference remains:

```text
/data/songlab/maeva/gregaria-timecourse/reference
```

## Download the fungal reference

The biological isolate was *Metarhizium robertsii* DWR2009, also known as
ARSEF 10343. Because an assembly under that exact isolate name was not found,
the workflow records ARSEF 23 (`GCF_000187425.2`) as a mapping proxy.

Download it into the project scratch reference folder:

```bash
cd /scratch/mtecher/gregaria-diet-infection-interaction/code
bash scripts/download_metarhizium_arsef23_hpc.sh
```

The standardized files will be:

```text
/scratch/mtecher/gregaria-diet-infection-interaction/reference/metarhizium/GCF_000187425.2/GCF_000187425.2_genomic.fna
/scratch/mtecher/gregaria-diet-infection-interaction/reference/metarhizium/GCF_000187425.2/GCF_000187425.2_genomic.gtf
```

The helper uses the equivalent NCBI Datasets command:

```bash
datasets download genome accession GCF_000187425.2 \
  --include genome,gtf,protein,rna,seq-report \
  --filename GCF_000187425.2_ncbi_dataset.zip \
  --no-progressbar
```

## Raw FASTQ input and trimming

The default starts from the original paired FASTQs:

```text
/data/songlab/sequencing_data/RNAseq/mehreen/mehreen_1002_MERGE_1.fq.gz
/data/songlab/sequencing_data/RNAseq/mehreen/mehreen_1002_MERGE_2.fq.gz
```

No trimmed reads need to exist before launch. Snakemake runs fastp independently
for every sample and writes the trimmed reads and QC reports inside the new run:

```text
01-trimmed-fastp/1002_1.trimmed.fastq.gz
01-trimmed-fastp/1002_2.trimmed.fastq.gz
01-trimmed-fastp/qc/1002.fastp.json
01-trimmed-fastp/qc/1002.fastp.html
```

The trimming settings match the time-course workflow: paired adapter detection,
2-bp front clipping on both mates, and a minimum retained length of 50 bp. This
single fresh set of trimmed reads is shared by both mapping branches.

## Two mapping branches

The host-only branch builds a run-local STAR index from
`GCF_023897955.1_iqSchGreg1.2` and maps every library only to the locust
reference. This provides the conventional host expression analysis and is the
primary comparison with the earlier study.

The competitive branch builds a separate STAR index from prefixed locust and
fungal references. It produces independent host and fungal count matrices and
reveals genes whose apparent host signal is reduced when reads can map to the
fungal proxy.

Both mappings use the same two-pass STAR settings, unique-read MAPQ convention,
paired featureCounts settings, and reverse-stranded library assumption.

For each host branch and for the competitive fungal branch, featureCounts uses
one active counting definition:

```text
-t transcript,exon -g gene_id
```

This keeps transcript-span and exon features together under their gene ID. The
pipeline does not generate a parallel exon-only matrix.

Before DESeq2, loci listed in
`data/excluded_loci/gregaria_rrna_list.txt` are removed. All biological samples
are retained, DESeq2 count replacement is disabled, and Cook's filtering is
disabled for the planned contrasts.

The optional `trimmed` input mode remains available for a future rerun that
intentionally reuses a complete fastp dataset:

```bash
export DUAL_INPUT_MODE=trimmed
export DUAL_TRIMMED_DIR=/path/to/completed/trimmed-fastp
```

## Kraken and Bracken

Install the Kraken2 and Bracken commands in a small environment on scratch:

```bash
cd /scratch/mtecher/gregaria-diet-infection-interaction/code
sbatch scripts/install_metatranscriptomics_qc_env.slurm
```

The installer creates or safely updates:

```text
/scratch/mtecher/conda-envs/metatranscriptomics-qc
```

It also writes `explicit-spec.txt` inside the environment to record the exact
installed package builds. The pipeline preflight activates this same
environment and verifies `kraken2` and `bracken` before worker submission.

The default Kraken database is:

```text
/scratch/mtecher/kraken2/pluspf_20260626
```

This is the versioned June 2026 PlusPF database, which contains RefSeq fungal
and protozoan genomes in addition to the standard Kraken content. The workflow
classifies:

1. all trimmed reads, to reproduce the colleague's family composition view;
2. reads unmapped to the host-only locust reference, enriching the pathogen and
   other non-host signal;
3. reads unmapped to both locust and *M. robertsii*, for more specific
   discovery of other symbionts or contaminants.

Clavicipitaceae is retained as a focal family in every composition plot even
when it is not among the globally most abundant families.

Kraken reports remain plain text for Bracken and summary plotting. The much
larger per-read classification files are written as `*.kraken2.output.gz` to
limit scratch use without discarding the read-level assignments.

Install the pinned database release on scratch:

```bash
cd /scratch/mtecher/gregaria-diet-infection-interaction/code
sbatch scripts/download_kraken_pluspf_20260626.slurm
```

The download uses one 24-hour job on the `public` partition. The archive is
downloaded with resumable `wget`, and `pigz` uses eight allocated cores for
extraction. Marker files preserve completed download, extraction, and checksum
stages if the job needs to be resubmitted. The download itself remains
network-limited and does not become faster with additional CPU cores.

Bracken is enabled by default for family-level estimated abundances. The
database needs a file such as:

```text
database150mers.kmer_distrib
```

If that file is unavailable, the pipeline can still make the same two stacked
barplot types directly from Kraken family-level clade counts:

```bash
export DUAL_RUN_BRACKEN=false
```

## HPC and local analysis boundary

The default Sol target is:

```text
host_pathogen_dual_hpc_preprocessing
```

It runs fastp, both STAR mappings, host and fungal featureCounts, Python mapping
comparisons, Kraken2, and optional Bracken. It does not run R, DESeq2,
enrichment, R-based plotting, or DEG-dependent unplaced-scaffold selection.
Those analyses are performed locally after transferring the compact result
tables back. The broader `host_pathogen_dual_all` target remains defined for
development, but it is not the default HPC submission target.

## Preflight and launch

Activate the Snakemake controller environment and run the preflight:

```bash
cd /scratch/mtecher/gregaria-diet-infection-interaction/code
module purge
module load mamba
source activate myENV
export DUAL_INPUT_MODE=raw
export DUAL_RAW_DIR=/data/songlab/sequencing_data/RNAseq/mehreen
export DUAL_KRAKEN_SOURCES=all_trimmed,host_unmapped,competitive_unmapped
bash pilot/check_host_pathogen_dual_hpc.sh

snakemake \
  -s Snakefile.host_pathogen_dual_rnaseq \
  host_pathogen_dual_hpc_preprocessing \
  --dry-run \
  --rerun-incomplete \
  --rerun-triggers mtime
```

Submit the HPC preprocessing workflow:

```bash
export DUAL_RUN_ID=host_pathogen_dual_$(date +%Y%m%d_%H%M%S)
unset DUAL_DIR
export DUAL_RESUME_EXISTING=false
export DUAL_INPUT_MODE=raw
export DUAL_RAW_DIR=/data/songlab/sequencing_data/RNAseq/mehreen
export DUAL_KRAKEN_SOURCES=all_trimmed,host_unmapped,competitive_unmapped
sbatch snakemake.host_pathogen_dual.slurm
```

The launcher refuses submissions without an explicit `DUAL_RUN_ID`. This
prevents a bare `sbatch` command from silently creating a second run when the
intent was to resume an existing one.

The launcher defaults to `host_pathogen_dual_hpc_preprocessing` and performs a
Snakemake dry run before submitting worker jobs. It does not require R or
DESeq2 packages in `myENV`. It permits up to 200 submitted or active Snakemake
worker jobs; as jobs finish, Snakemake submits later jobs until the DAG is
complete. SLURM may run fewer according to available resources and account
limits. The controller requests `6-23:59:00` on `public` and uses
`--keep-going`, so independent samples continue after an individual rule
failure. Each new launch receives a new output folder:

```text
/scratch/mtecher/gregaria-diet-infection-interaction/output/runs/host_pathogen_dual_YYYYMMDD_HHMMSS
```

Scheduler profiles use `htc` only for jobs with walltimes of four hours or
less. Jobs requiring more than four hours use `public`, with requested
walltimes kept below the seven-day limit.

### Resume an interrupted run

Do not start a second controller while the first controller or any of its
workers is active. After confirming with `squeue --me` that the workflow has
stopped, resume the same run ID as follows:

```bash
export DUAL_RUN_ID=host_pathogen_dual_YYYYMMDD_HHMMSS
export DUAL_DIR=/scratch/mtecher/gregaria-diet-infection-interaction/output/runs/${DUAL_RUN_ID}
export DUAL_RESUME_EXISTING=true

sbatch --export=ALL snakemake.host_pathogen_dual.slurm
```

The launcher uses `--rerun-incomplete --rerun-triggers mtime`. Snakemake keeps
complete outputs and schedules only missing outputs or jobs marked incomplete.
Never delete the existing run folder merely to resume it.

### Monitor a run without interrupting it

The progress helper only counts non-empty outputs and queries SLURM. It does
not invoke Snakemake, acquire its lock, submit jobs, unlock the run, or alter
files:

```bash
cd /scratch/mtecher/gregaria-diet-infection-interaction/code
bash pilot/check_host_pathogen_dual_progress.sh \
  host_pathogen_dual_20260727_000952
```

Refresh the same report every minute:

```bash
watch -n 60 \
  'bash pilot/check_host_pathogen_dual_progress.sh host_pathogen_dual_20260727_000952'
```

Follow the controller and worker logs separately:

```bash
squeue -u "$USER" -o "%.18i %.32j %.10T %.10M %.30R"
tail -f slurm-fatbody-dual-JOBID.out
find /scratch/mtecher/gregaria-diet-infection-interaction/logs/slurm/host_pathogen_dual_20260727_000952 \
  -type f -name "*.err" -size +0c -print
```

While the controller is running, the stage counter is the safe way to inspect
what remains. To obtain Snakemake's exact remaining rule DAG, first confirm
that neither the controller nor any worker is active. Then export the same run
ID and run the dry-run command shown above. Do not run a second Snakemake
process against a locked, active run.

## Main HPC outputs

```text
00-reference/
  combined_host_metarhizium.fna
  combined_sequence_id_map.tsv
  STAR-host-only/
  STAR-competitive/
01-trimmed-fastp/
  paired trimmed FASTQs
  qc/*.fastp.html
  qc/*.fastp.json
02-host-only-star/
  fresh host-only BAMs, CSI indexes, STAR logs, and host-unmapped FASTQs
02-competitive-star/
  fresh competitive BAMs, CSI indexes, STAR logs, and competitive-unmapped FASTQs
03-featurecounts/
  host-only transcript+exon tables
  competitive host transcript+exon tables
  competitive fungal transcript+exon tables
04-count-matrices/
  host-only/
  competitive-host/
  competitive-fungus/
05-mapping-comparison/
  host_only_mapping_summary.tsv
  competitive_mapping_summary.tsv
  host_only_vs_competitive_sample_counts.tsv
  host_only_vs_competitive_gene_counts.tsv
08-taxonomy/
  all_trimmed/
  host_unmapped/
  competitive_unmapped/
11-local-transfer/
  transfer_manifest.tsv
  fatbody_hpc_results_for_local.tar.gz
```

## Transfer back for local analysis

After the HPC target finishes, the simplest Globus transfer is this single
archive:

```text
11-local-transfer/fatbody_hpc_results_for_local.tar.gz
```

The accompanying `transfer_manifest.tsv` records the source path, archive path,
file size, and SHA-256 checksum for every included file. The archive contains:

```text
run_provenance.tsv
00-reference/combined_sequence_id_map.tsv
00-reference/combined_gene_id_map.tsv
01-trimmed-fastp/qc/
02-host-only-star/*_Log.final.out
02-competitive-star/*_Log.final.out
03-featurecounts/*/*.summary
04-count-matrices/
05-mapping-comparison/
08-taxonomy/family_abundance_manifest.tsv
08-taxonomy/*/*.kraken2.report
08-taxonomy/*/*.bracken.family.tsv
```

After Globus transfers the archive, extract it inside a new local run folder:

```bash
mkdir -p output/runs/host_pathogen_dual_YYYYMMDD_HHMMSS
tar -xzf fatbody_hpc_results_for_local.tar.gz \
  -C output/runs/host_pathogen_dual_YYYYMMDD_HHMMSS
```

The BAMs, trimmed FASTQs, STAR indexes, and compressed per-read Kraken outputs
can remain on scratch unless a later analysis specifically needs them. DESeq2
requires the count matrices, sample metadata, rRNA exclusion list, and analysis
scripts; it does not require the BAMs or FASTQs.

Local downstream outputs will be written under `06-host-deseq2/`,
`07-metarhizium-deseq2/`, `09-unplaced-origin/`, and
`10-final-evidence-catalogue/`.

The host DESeq2 analyses use the same six-group model for both mappings. Planned
results include infection within each diet, diet contrasts among controls, diet
contrasts among infected animals, the average infection response, and three
diet-by-infection interaction contrasts. Positive log2 fold change always
refers to the numerator recorded in the output table.

The final reconciliation labels significant genes as supported by both
mappings in the same direction, supported by both in opposite directions,
host-only only, or competitive-host only. Every fresh host-only DEG on an
`NW_` scaffold enters the unplaced-scaffold origin screen automatically; the
reconciliation table then shows whether its differential-expression evidence
is retained under competitive mapping.

Fungal differential expression is attempted only among infected libraries that
pass configurable fungal-read and detected-gene thresholds. Control libraries
are retained for background assessment but are not treated as zero-expression
fungal replicates.

The unplaced-scaffold Kraken table is a screening result, not final proof of
origin. Final host/non-host labels should combine competitive alignment,
sequence taxonomy, arthropod/fungal protein homology, repeats, and mapping
uniqueness.

No tissue remains for qPCR validation. Fungal RNA abundance must therefore be
reported as a tissue-specific transcriptional proxy for fungal burden, not as
an absolute measurement of fungal biomass.
