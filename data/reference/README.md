# Reference annotation

The reports use the NCBI RefSeq annotation for assembly
`GCF_023897955.1_iqSchGreg1.2`.

Git tracks `GCF_023897955.1_iqSchGreg1.2_genomic.gff.gz`; the expanded `.gff`
is intentionally ignored because it is approximately 304 MB. Prepare the
working copy after cloning with:

```bash
bash scripts/setup_reference_data.sh
```

The setup script checks the gzip stream and refuses to replace a non-matching
existing GFF.
