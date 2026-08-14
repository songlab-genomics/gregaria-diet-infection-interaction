# Data

Small design files required to understand the experiment remain in this GitHub
repository:

- `metadata/`: corrected sample assignments and body-mass metadata.
- `excluded_loci/`: the *S. gregaria* rRNA exclusion list.

Large count tables, functional annotations, reference annotations, external
comparison datasets, and scaffold-audit inputs are stored in the associated
Dryad data package. Run `bash scripts/setup_dryad_links.sh` after downloading
that package. The script recreates the expected `data/...` links without
changing paths in the R Markdown analyses.
