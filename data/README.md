# Data

This directory contains the analysis-ready inputs used by the workflowR
reports:

- `metadata/`: corrected sample assignments and body-mass metadata.
- `excluded_loci/`: *S. gregaria* rRNA and lncRNA identifier lists.
- `raw_read_counts/`: original and merged featureCounts tables and historical
  comparison tables.
- `derived/`: frozen count matrices, DESeq2 results, mapping summaries,
  taxonomy reports, and scaffold-audit evidence promoted from authoritative
  analysis runs.
- `GO_Annotations/`: eggNOG and supporting functional annotations.
- `external/`: phase-reference, time-course, and manually curated inputs.
- `scaffold_origin/`: taxonomy lookup material used by the origin audit.
- `reference/`: assembly reports and the compressed NCBI GFF.

Run `bash scripts/setup_reference_data.sh` once after cloning to expand the GFF
used by the R Markdown pages. `SHA256SUMS` records the checksum of every bundled
file except the generated uncompressed GFF and the manifest itself.
