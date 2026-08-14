# Gregaria Diet-Infection Data Package

This directory contains the large inputs and full analysis outputs associated
with the `songlab-genomics/gregaria-diet-infection-interaction` GitHub
repository. The GitHub repository contains the runnable workflowR and
Snakemake code, corrected metadata, rRNA exclusion list, and rendered website.

## Package layout

- `data/GO_Annotations/`: eggNOG and supporting functional annotations.
- `data/raw_read_counts/`: original and merged gene-count inputs.
- `data/reference/`: the NCBI *Schistocerca gregaria* genome annotation.
- `data/external/`: phase-reference and manually curated candidate inputs.
- `data/scaffold_origin/`: taxonomy lookup inputs used in the scaffold audit.
- `data/legacy_csv/`: retained historical comparison tables.
- `output/runs/`: transferred HPC products and mapping/scaffold audits.
- `output/rmd_runs/`: timestamped complete outputs from the R Markdown pages.
- `outputs/`: paper supplementary workbooks and validation products.
- `legacy_website_assets/`: older data-heavy website assets and their
  standalone HTML report, which are not used by the active site.

`dryad_file_manifest.tsv` records the size and SHA-256 checksum of each file.

## Reconnect the data to the analysis repository

Clone the GitHub repository, download this package, and run:

```bash
export GREGARIA_DIET_DRYAD_DIR=/path/to/gregaria-diet-infection-interaction-dryad
bash scripts/setup_dryad_links.sh
```

The setup script creates links at the project-relative paths already documented
in the R Markdown parameters. It does not run analyses or modify deposited
files.
