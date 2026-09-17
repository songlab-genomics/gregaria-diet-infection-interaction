# Output

This directory is reserved for newly generated local results. Each R Markdown
render creates a fresh timestamped folder under `output/rmd_runs/`; Snakemake
and HPC workflows use `output/runs/`.

These generated trees are ignored by Git. Analysis-ready files required by the
published reports are frozen under `data/derived/`, while the rendered website
and publication-facing figures are committed under `docs/`.
