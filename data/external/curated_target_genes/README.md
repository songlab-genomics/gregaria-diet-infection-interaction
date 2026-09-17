# Curated target-gene Figure 4 inputs

These files are frozen copies of the coauthor materials supplied on 2026-08-05.
They are retained unchanged so that the original Figure 4 can be compared with
the main-chromosome reanalysis.

| File | Role | SHA-256 |
| --- | --- | --- |
| `2026-07-26_curated-hypoth_DEGs_v14.xlsx` | Curated immune and protein-anabolism target list plus the original contrast values | `d5d912da409d95691cb951e1bd3cb7f3240c462894fbc34f8ed28c7580da7995` |
| `CuratedDEGs_CompositeFig_v03.pdf` | Original composite figure used as the layout reference | `632a280b5ca7d60ab6ffe72bc90184481395ba5837f2157c52e340a51c752c37` |
| `hypoth-curated_DEGs_figures.docx` | Coauthor draft headings for the curated DEG figures | `823683759389950e92375e3edf13b8a32c02e2fc291359b9e06e98af5c638a73` |

The original R plotting script is stored under
`scripts/reference/curated_target_figure_v14_original.R` with SHA-256
`02b1834646caf1a19734b5d918a1123b46b073097ffa70662e5def3503404d20`.

The active analysis is `analysis/30-curated-target-gene-figure.Rmd`. It does not
reuse the pre-summarized workbook counts. Instead, it keeps the target-gene
categories, joins each gene to the official assembly-derived placement table,
and recalculates every plotted status from the current competitive-host DESeq2
results. Genes on unplaced scaffolds are exported for follow-up but are not
included in the primary Figure 4 percentages.
