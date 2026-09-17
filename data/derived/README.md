# Frozen analysis inputs

These files were promoted from completed, authoritative analysis runs so the
workflowR reports do not depend on timestamped local output directories.

- `counts/`: rRNA-filtered competitive and host-only matrices for the primary
  transcript-plus-exon analysis and the exon-only sensitivity analysis.
- `deseq2/`: complete contrast tables and compact model-audit summaries.
- `mapping/`: gene placement, competitive mapping summaries, and STAR final
  logs used by quality-control pages.
- `taxonomy/`: Kraken2 reports, Bracken family tables, and the 45-library
  manifest used for the microbial-composition report.
- `scaffold_audit/`: transferred cross-species, Kraken/FCS, BLASTn, and DIAMOND
  evidence used to assess unresolved scaffolds.
- `analysis_inputs/`: frozen cluster assignments, homology evidence, curated
  gene audits, and compressed placed/unplaced result tables.

The biological analysis cohort contains 44 retained libraries; sample 1044 is
included only in quality-control and taxonomy source material. All files are
covered by `data/SHA256SUMS`.
