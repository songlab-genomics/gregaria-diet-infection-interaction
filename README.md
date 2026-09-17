# Effect of Diet on Metarhizium Infection in S. Gregaria

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr

## Data and reproducibility

GitHub contains the analysis code, rendered website, and the analysis-ready
inputs needed by the workflowR reports. These include corrected metadata,
rRNA exclusions, host count matrices, DESeq2 result tables, eggNOG annotations,
phase-reference datasets, taxonomy summaries, and scaffold-origin audit data.

The NCBI GFF is committed as a compressed file to keep the repository compact.
After cloning, prepare its uncompressed working copy once:

```bash
bash scripts/setup_reference_data.sh
```

The expanded GFF is ignored by Git, while every R Markdown page continues to
use the readable project-relative path `data/reference/...genomic.gff`.
Bundled files can be checked against the committed manifest with:

```bash
shasum -a 256 -c data/SHA256SUMS
```

Raw FASTQ files, BAM files, complete Kraken classifications, and historical
timestamped result trees are intentionally not stored in Git. New report runs
are written locally to `output/rmd_runs/`; publication-facing HTML and figure
assets are committed under `docs/`.

## Zenodo code archive

Create the Zenodo code, data, and website zip from a clean commit rather than zipping
the working directory directly:

```bash
bash scripts/build_zenodo_release.sh v1.0.0
```

The script uses `git archive`, validates the resulting zip, and reports its
SHA-256 checksum. It includes the committed analysis-ready data, code, metadata,
and rendered site. Local R history files, `.git`, and generated output folders
remain excluded.

## GitHub Pages deployment

The rendered workflowR website is committed under `docs/`. In the repository
Pages settings, use **Deploy from a branch**, with branch `main` and folder
`/docs`. The included Pages workflow is retained as a manual fallback only; use
it only after intentionally changing the repository Pages source to
**GitHub Actions**. Keeping both deployment modes active on every push can
create competing Pages deployments.
