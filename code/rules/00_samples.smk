### =================================================================
### SAMPLE DISCOVERY
### =================================================================

# Detect samples per tissue based on your naming convention.
# Expected raw-read input:
#   {DATADIR}/00-reads-{tissue}/{locust}_MERGE_1.fastq.gz
#   {DATADIR}/00-reads-{tissue}/{locust}_MERGE_2.fastq.gz
# Expected processed-BAM input:
#   {WORKDIR}/02-{tissue}-star/{locust}_Aligned.sortedByCoord.out.bam

def discover_locusts_for_tissue(tissue):
    """Find samples from raw reads or existing STAR BAMs.

    Use SAMPLE_DISCOVERY = "raw_reads" to force FASTQ discovery,
    "star_bam" to force completed STAR BAM discovery, or "auto" to use raw
    reads first and fall back to BAMs. This lets variant calling reuse the
    completed /data/songlab/maeva/{species}-timecourse/data alignments.
    """
    raw_samples = glob_wildcards(
        DATADIR + f"/00-reads-{tissue}/{{locust}}_MERGE_1.fastq.gz"
    ).locust
    bam_samples = glob_wildcards(
        WORKDIR + f"/02-{tissue}-star/{{locust}}_Aligned.sortedByCoord.out.bam"
    ).locust

    if SAMPLE_DISCOVERY == "raw_reads":
        return raw_samples
    if SAMPLE_DISCOVERY == "star_bam":
        return bam_samples
    if len(raw_samples) > 0:
        return raw_samples
    return bam_samples


LOCUSTS = {tissue: discover_locusts_for_tissue(tissue) for tissue in TISSUES}

print("Detected LOCUSTS per tissue:")
for t, ls in LOCUSTS.items():
    print(f"  {t}: {len(ls)} samples")
