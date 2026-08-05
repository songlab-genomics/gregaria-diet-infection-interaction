# Scaffold-origin reference data

`ncbi_taxonomy_lookup_20260804.tsv` is a compact lookup for the taxids present
in the transferred BLASTn result. Scientific names, ranks, and lineages were
derived from the NCBI Taxonomy dump downloaded on 2026-08-04 from:

`https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`

The lookup is used by `analysis/27-scaffold-origin-audit.Rmd` to classify
external nucleotide hits after removing *Schistocerca gregaria* taxid 7010.
Its SHA-256 checksum is:

`060785edf48455de0d6d010f607b66979dddd678127f97da17ff21986239d73e`
