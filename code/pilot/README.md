# Fat Body Pilot Run

This pilot keeps the copied time-course Snakemake rules unchanged and adapts the
Mehreen fat body raw-read layout with symlinks.

The copied time-course rules expect paired reads in this layout:

```text
{DATADIR}/00-reads-{tissue}/{sample}_MERGE_1.fastq.gz
{DATADIR}/00-reads-{tissue}/{sample}_MERGE_2.fastq.gz
```

The raw fat body files are currently named like:

```text
/data/songlab/sequencing_data/RNAseq/mehreen/mehreen_1002_MERGE_1.fq.gz
```

The pilot launcher creates run-specific symlinks with the expected names, then
runs `expression_all`, which covers fastp trimming, STAR alignment, and
featureCounts for only the samples in `fatbody_pilot_samples.tsv`.

Run the preflight from the `code/` folder on the cluster:

```bash
cd /scratch/mtecher/gregaria-diet-infection-interaction/code
bash pilot/check_fatbody_pilot_hpc.sh
```

If the dry-run succeeds, submit the pilot:

```bash
cd /scratch/mtecher/gregaria-diet-infection-interaction/code
sbatch snakemake.fatbody_pilot.slurm
```

Optional arguments:

```bash
sbatch snakemake.fatbody_pilot.slurm expression_all 12 myENV
```

Important inputs to verify before launching:

- Raw reads: `/data/songlab/sequencing_data/RNAseq/mehreen`
- Preferred reference folder: `/data/songlab/maeva/gregaria-timecourse/reference`
- Fallback reference folders: `/scratch/mtecher/gregaria-timecourse/reference` and copied `locust-time-course-RNAseq` reference folders
- Output root: `/scratch/mtecher/gregaria-diet-infection-interaction`

Cluster-worker logs are written under:

```text
/scratch/mtecher/gregaria-diet-infection-interaction/logs/slurm/{RUN_ID}/
```

The first pilot failed at `trimming_fastp`, before STAR or featureCounts. The
trim rule now checks the active `fastp` executable and lets you override the tool
without editing the Snakefile:

```bash
FASTP_CONDA_ENV=myENV sbatch snakemake.fatbody_pilot.slurm
FASTP_MODULE=fastp-0.23.4-gcc-12.1.0 FASTP_CONDA_ENV= sbatch snakemake.fatbody_pilot.slurm
FASTP_CMD=/path/to/fastp sbatch snakemake.fatbody_pilot.slurm
```

The copied time-course `cluster.json` uses `--export=NONE` for ordinary jobs.
For this pilot, the launcher overrides that setting because worker jobs must see
the run-specific scratch paths (`DATADIR`, `WORKDIR`, `REFDIR`, and related
variables). If workers fall back to `/data/songlab/sequencing_data/RNAseq/timecourse/gregaria`
or `/data/songlab/maeva/gregaria-timecourse/data`, the export override is not
being used.
