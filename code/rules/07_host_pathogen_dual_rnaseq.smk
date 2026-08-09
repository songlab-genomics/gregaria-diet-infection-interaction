### =================================================================
### FAT-BODY HOST-PATHOGEN DUAL RNA-SEQ RULES
### =================================================================

import csv
import hashlib
import os
import tarfile
from pathlib import Path


ALLOWED_DUAL_INPUT_MODES = {"trimmed", "raw"}
ALLOWED_DUAL_KRAKEN_SOURCES = {
    "all_trimmed",
    "host_unmapped",
    "competitive_unmapped",
}

if DUAL_INPUT_MODE not in ALLOWED_DUAL_INPUT_MODES:
    raise ValueError(
        f"DUAL_INPUT_MODE must be one of {sorted(ALLOWED_DUAL_INPUT_MODES)}; "
        f"received {DUAL_INPUT_MODE}"
    )
if not set(DUAL_KRAKEN_SOURCES).issubset(ALLOWED_DUAL_KRAKEN_SOURCES):
    raise ValueError(
        "DUAL_KRAKEN_SOURCES contains unsupported values: "
        + ", ".join(sorted(set(DUAL_KRAKEN_SOURCES) - ALLOWED_DUAL_KRAKEN_SOURCES))
    )


with open(DUAL_SAMPLE_TABLE, newline="") as handle:
    DUAL_RECORDS = list(csv.DictReader(handle, delimiter="\t"))

required_sample_fields = {"sample_id", "raw_prefix", "tissue", "treatment", "diet"}
missing_sample_fields = required_sample_fields - set(DUAL_RECORDS[0] if DUAL_RECORDS else {})
if missing_sample_fields:
    raise ValueError(
        f"DUAL_SAMPLE_TABLE is missing fields: {', '.join(sorted(missing_sample_fields))}"
    )

DUAL_SAMPLES = [record["sample_id"] for record in DUAL_RECORDS]
DUAL_RECORD_BY_SAMPLE = {record["sample_id"]: record for record in DUAL_RECORDS}

if len(DUAL_RECORD_BY_SAMPLE) != len(DUAL_SAMPLES):
    duplicate_samples = sorted(
        sample for sample in set(DUAL_SAMPLES) if DUAL_SAMPLES.count(sample) > 1
    )
    raise ValueError(
        "DUAL_SAMPLE_TABLE contains duplicate sample IDs: "
        + ", ".join(duplicate_samples)
    )

unexpected_treatments = sorted(
    {record["treatment"] for record in DUAL_RECORDS} - {"Control", "Infected"}
)
unexpected_diets = sorted(
    {record["diet"] for record in DUAL_RECORDS} - {"33", "50", "83"}
)
if unexpected_treatments or unexpected_diets:
    raise ValueError(
        "Unexpected metadata values. Treatments: "
        + ", ".join(unexpected_treatments or ["none"])
        + "; diets: "
        + ", ".join(unexpected_diets or ["none"])
    )

DUAL_COMBINED_FASTA = DUAL_REFERENCE_DIR + "/combined_host_metarhizium.fna"
DUAL_COMBINED_GTF = DUAL_REFERENCE_DIR + "/combined_host_metarhizium.gtf"
DUAL_HOST_PREFIXED_GTF = DUAL_REFERENCE_DIR + "/host.prefixed.gtf"
DUAL_FUNGUS_PREFIXED_GTF = DUAL_REFERENCE_DIR + "/metarhizium.prefixed.gtf"
DUAL_SEQUENCE_MAP = DUAL_REFERENCE_DIR + "/combined_sequence_id_map.tsv"
DUAL_GENE_MAP = DUAL_REFERENCE_DIR + "/combined_gene_id_map.tsv"
DUAL_PROVENANCE = DUAL_DIR + "/run_provenance.tsv"

DUAL_HOST_ONLY_PRIMARY_MATRIX = (
    DUAL_DIR
    + "/04-count-matrices/host-only/host_transcript_exon_counts.tsv"
)
DUAL_COMPETITIVE_HOST_PRIMARY_MATRIX = (
    DUAL_DIR
    + "/04-count-matrices/competitive-host/host_transcript_exon_counts.tsv"
)
DUAL_FUNGUS_PRIMARY_MATRIX = (
    DUAL_DIR
    + "/04-count-matrices/competitive-fungus/metarhizium_transcript_exon_counts.tsv"
)
# Secondary host matrices support the explicit exon-only sensitivity pages.
# They never replace the transcript+exon primary analysis.
DUAL_HOST_ONLY_EXON_MATRIX = (
    DUAL_DIR + "/04-count-matrices/host-only/host_exon_counts.tsv"
)
DUAL_COMPETITIVE_HOST_EXON_MATRIX = (
    DUAL_DIR + "/04-count-matrices/competitive-host/host_exon_counts.tsv"
)
DUAL_HOST_ONLY_MAPPING_QC = (
    DUAL_DIR + "/05-mapping-comparison/host_only_mapping_summary.tsv"
)
DUAL_COMPETITIVE_MAPPING_QC = (
    DUAL_DIR + "/05-mapping-comparison/competitive_mapping_summary.tsv"
)
DUAL_TAXONOMY_MANIFEST = DUAL_DIR + "/08-taxonomy/family_abundance_manifest.tsv"

DUAL_HOST_DE_FILENAMES = [
    "host_counts_no_rrna.tsv",
    "host_deseq2_filter_audit.tsv",
    "host_deseq2_all_contrasts.csv",
    "host_deseq2_contrast_summary.tsv",
    "host_transcriptome_pca.png",
    "host_deg_counts.png",
]
DUAL_HOST_ONLY_DE_DIR = DUAL_DIR + "/06-host-deseq2/host-only"
DUAL_COMPETITIVE_HOST_DE_DIR = DUAL_DIR + "/06-host-deseq2/competitive-host"
DUAL_HOST_ONLY_DE_OUTPUTS = [
    DUAL_HOST_ONLY_DE_DIR + "/" + filename for filename in DUAL_HOST_DE_FILENAMES
]
DUAL_COMPETITIVE_HOST_DE_OUTPUTS = [
    DUAL_COMPETITIVE_HOST_DE_DIR + "/" + filename
    for filename in DUAL_HOST_DE_FILENAMES
]

DUAL_FUNGAL_DE_DIR = DUAL_DIR + "/07-metarhizium-deseq2"
DUAL_FUNGAL_DE_OUTPUTS = [
    DUAL_FUNGAL_DE_DIR + "/fungal_deseq2_sample_eligibility.tsv",
    DUAL_FUNGAL_DE_DIR + "/fungal_deseq2_all_contrasts.csv",
    DUAL_FUNGAL_DE_DIR + "/fungal_deseq2_contrast_summary.tsv",
    DUAL_FUNGAL_DE_DIR + "/fungal_transcriptome_pca.png",
    DUAL_FUNGAL_DE_DIR + "/fungal_deg_counts_by_diet.png",
]

DUAL_COUNT_RECONCILIATION_SAMPLE = (
    DUAL_DIR + "/05-mapping-comparison/host_only_vs_competitive_sample_counts.tsv"
)
DUAL_COUNT_RECONCILIATION_GENE = (
    DUAL_DIR + "/05-mapping-comparison/host_only_vs_competitive_gene_counts.tsv"
)
DUAL_DE_RECONCILIATION = (
    DUAL_DIR + "/10-final-evidence-catalogue/host_mapping_deg_reconciliation.csv"
)
DUAL_DE_RECONCILIATION_SUMMARY = (
    DUAL_DIR + "/10-final-evidence-catalogue/host_mapping_deg_reconciliation_summary.tsv"
)
DUAL_FINAL_ORIGIN_EVIDENCE = (
    DUAL_DIR + "/10-final-evidence-catalogue/host_deg_origin_evidence_catalogue.tsv"
)
DUAL_TRANSFER_MANIFEST = DUAL_DIR + "/11-local-transfer/transfer_manifest.tsv"
DUAL_TRANSFER_BUNDLE = (
    DUAL_DIR + "/11-local-transfer/fatbody_hpc_results_for_local.tar.gz"
)

DUAL_UNPLACED_DEG_GENES = (
    DUAL_DIR + "/09-unplaced-origin/host_only_unplaced_deg_genes.tsv"
)
DUAL_UNPLACED_SCAFFOLD_LIST = (
    DUAL_DIR + "/09-unplaced-origin/host_only_unplaced_deg_scaffolds.txt"
)
DUAL_UNPLACED_FASTA = (
    DUAL_DIR + "/09-unplaced-origin/candidate_unplaced_scaffolds.fna"
)
DUAL_UNPLACED_METRICS = DUAL_DIR + "/09-unplaced-origin/candidate_unplaced_scaffolds.tsv"
DUAL_UNPLACED_KRAKEN_REPORT = (
    DUAL_DIR + "/09-unplaced-origin/candidate_unplaced_scaffolds.kraken2.report"
)
DUAL_UNPLACED_KRAKEN_OUTPUT = (
    DUAL_DIR + "/09-unplaced-origin/candidate_unplaced_scaffolds.kraken2.output"
)
DUAL_UNPLACED_TAXONOMY = (
    DUAL_DIR + "/09-unplaced-origin/candidate_unplaced_scaffold_taxonomy.tsv"
)

DUAL_FAMILY_PLOT_OUTPUTS = [
    DUAL_DIR + f"/08-taxonomy/figures/{source}/{filename}"
    for source in DUAL_KRAKEN_SOURCES
    for filename in [
        "Family_Composition_Reads.png",
        "Family_Composition_Reads.svg",
        "Family_Composition_Percent.png",
        "Family_Composition_Percent.svg",
        "family_composition_plot_data.tsv",
    ]
]


def dual_raw_read(wildcards, mate):
    record = DUAL_RECORD_BY_SAMPLE[wildcards.sample]
    return f"{DUAL_RAW_DIR}/{record['raw_prefix']}_MERGE_{mate}.fq.gz"


def dual_trimmed_read(wildcards, mate):
    if DUAL_INPUT_MODE == "raw":
        return (
            DUAL_DIR
            + f"/01-trimmed-fastp/{wildcards.sample}_{mate}.trimmed.fastq.gz"
        )
    return f"{DUAL_TRIMMED_DIR}/{wildcards.sample}_{mate}.trimmed.fastq.gz"


def dual_kraken_read(wildcards, mate):
    if wildcards.source == "all_trimmed":
        return dual_trimmed_read(wildcards, mate)
    if wildcards.source == "host_unmapped":
        return (
            DUAL_DIR
            + f"/02-host-only-star/unmapped/{wildcards.sample}_host_unmapped_R{mate}.fastq.gz"
        )
    return (
        DUAL_DIR
        + f"/02-competitive-star/unmapped/{wildcards.sample}_competitive_unmapped_R{mate}.fastq.gz"
    )


def dual_host_only_star_bams():
    return expand(
        DUAL_DIR + "/02-host-only-star/{sample}_Aligned.sortedByCoord.out.bam",
        sample=DUAL_SAMPLES,
    )


def dual_competitive_star_bams():
    return expand(
        DUAL_DIR + "/02-competitive-star/{sample}_Aligned.sortedByCoord.out.bam",
        sample=DUAL_SAMPLES,
    )


def dual_featurecounts(pattern):
    return expand(pattern, sample=DUAL_SAMPLES)


def dual_kraken_reports():
    return expand(
        DUAL_DIR + "/08-taxonomy/{source}/{sample}.kraken2.report",
        source=DUAL_KRAKEN_SOURCES,
        sample=DUAL_SAMPLES,
    )


def dual_bracken_outputs():
    if not DUAL_RUN_BRACKEN:
        return []
    return expand(
        DUAL_DIR + "/08-taxonomy/{source}/{sample}.bracken.family.tsv",
        source=DUAL_KRAKEN_SOURCES,
        sample=DUAL_SAMPLES,
    )


def dual_local_transfer_files():
    """Small HPC results needed for local QC, plotting, and DESeq2."""
    featurecounts_summaries = [
        DUAL_DIR + f"/03-featurecounts/{count_type}/{sample}.featureCounts.txt.summary"
        for count_type in [
            "host-only-transcript-exon",
            "competitive-host-transcript-exon",
            "competitive-fungus-transcript-exon",
            "host-only-exon",
            "competitive-host-exon",
        ]
        for sample in DUAL_SAMPLES
    ]
    return (
        [
            DUAL_PROVENANCE,
            DUAL_SEQUENCE_MAP,
            DUAL_GENE_MAP,
            DUAL_HOST_ONLY_PRIMARY_MATRIX,
            DUAL_COMPETITIVE_HOST_PRIMARY_MATRIX,
            DUAL_FUNGUS_PRIMARY_MATRIX,
            DUAL_HOST_ONLY_EXON_MATRIX,
            DUAL_COMPETITIVE_HOST_EXON_MATRIX,
            DUAL_HOST_ONLY_MAPPING_QC,
            DUAL_COMPETITIVE_MAPPING_QC,
            DUAL_COUNT_RECONCILIATION_SAMPLE,
            DUAL_COUNT_RECONCILIATION_GENE,
            DUAL_TAXONOMY_MANIFEST,
        ]
        + expand(
            DUAL_DIR + "/01-trimmed-fastp/qc/{sample}.fastp.json",
            sample=DUAL_SAMPLES,
        )
        + expand(
            DUAL_DIR + "/01-trimmed-fastp/qc/{sample}.fastp.html",
            sample=DUAL_SAMPLES,
        )
        + expand(
            DUAL_DIR + "/02-host-only-star/{sample}_Log.final.out",
            sample=DUAL_SAMPLES,
        )
        + expand(
            DUAL_DIR + "/02-competitive-star/{sample}_Log.final.out",
            sample=DUAL_SAMPLES,
        )
        + featurecounts_summaries
        + dual_kraken_reports()
        + dual_bracken_outputs()
    )


def dual_rnaseq_final_outputs():
    return (
        [
            DUAL_PROVENANCE,
            DUAL_COMBINED_FASTA,
            DUAL_SEQUENCE_MAP,
            DUAL_GENE_MAP,
            DUAL_HOST_ONLY_PRIMARY_MATRIX,
            DUAL_COMPETITIVE_HOST_PRIMARY_MATRIX,
            DUAL_FUNGUS_PRIMARY_MATRIX,
            DUAL_HOST_ONLY_EXON_MATRIX,
            DUAL_COMPETITIVE_HOST_EXON_MATRIX,
            DUAL_HOST_ONLY_MAPPING_QC,
            DUAL_COMPETITIVE_MAPPING_QC,
            DUAL_COUNT_RECONCILIATION_SAMPLE,
            DUAL_COUNT_RECONCILIATION_GENE,
            DUAL_TAXONOMY_MANIFEST,
            DUAL_TRANSFER_MANIFEST,
            DUAL_TRANSFER_BUNDLE,
            DUAL_UNPLACED_DEG_GENES,
            DUAL_UNPLACED_TAXONOMY,
            DUAL_DE_RECONCILIATION,
            DUAL_DE_RECONCILIATION_SUMMARY,
            DUAL_FINAL_ORIGIN_EVIDENCE,
        ]
        + dual_host_only_star_bams()
        + dual_competitive_star_bams()
        + dual_kraken_reports()
        + dual_bracken_outputs()
        + DUAL_FAMILY_PLOT_OUTPUTS
        + DUAL_HOST_ONLY_DE_OUTPUTS
        + DUAL_COMPETITIVE_HOST_DE_OUTPUTS
        + DUAL_FUNGAL_DE_OUTPUTS
    )


def dual_hpc_preprocessing_outputs():
    """Compute-heavy outputs produced on Sol before local R analysis."""
    return (
        [
            DUAL_PROVENANCE,
            DUAL_COMBINED_FASTA,
            DUAL_SEQUENCE_MAP,
            DUAL_GENE_MAP,
            DUAL_HOST_ONLY_PRIMARY_MATRIX,
            DUAL_COMPETITIVE_HOST_PRIMARY_MATRIX,
            DUAL_FUNGUS_PRIMARY_MATRIX,
            DUAL_HOST_ONLY_MAPPING_QC,
            DUAL_COMPETITIVE_MAPPING_QC,
            DUAL_COUNT_RECONCILIATION_SAMPLE,
            DUAL_COUNT_RECONCILIATION_GENE,
            DUAL_TAXONOMY_MANIFEST,
            DUAL_TRANSFER_MANIFEST,
            DUAL_TRANSFER_BUNDLE,
        ]
        + dual_host_only_star_bams()
        + dual_competitive_star_bams()
        + dual_kraken_reports()
        + dual_bracken_outputs()
    )


rule host_pathogen_dual_hpc_preprocessing:
    input:
        dual_hpc_preprocessing_outputs()


rule host_pathogen_dual_all:
    input:
        dual_rnaseq_final_outputs()


rule dual_run_provenance:
    input:
        sample_table=DUAL_SAMPLE_TABLE,
        host_fasta=DUAL_HOST_FASTA,
        host_gtf=DUAL_HOST_GTF,
        fungus_fasta=DUAL_FUNGUS_FASTA,
        fungus_gtf=DUAL_FUNGUS_GTF,
        rrna=DUAL_RRNA_LIST
    output:
        DUAL_PROVENANCE
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        rows = [
            ("run_id", DUAL_RUN_ID),
            ("run_directory", DUAL_DIR),
            ("input_mode", DUAL_INPUT_MODE),
            ("raw_read_directory", DUAL_RAW_DIR),
            ("trimmed_read_directory", DUAL_TRIMMED_DIR),
            ("sample_table", DUAL_SAMPLE_TABLE),
            ("host_fasta", DUAL_HOST_FASTA),
            ("host_gtf", DUAL_HOST_GTF),
            ("experimental_fungal_isolate", DUAL_EXPERIMENTAL_ISOLATE),
            ("fungal_mapping_proxy", DUAL_FUNGUS_ISOLATE),
            ("fungal_reference_accession", DUAL_FUNGUS_ACCESSION),
            ("fungal_fasta", DUAL_FUNGUS_FASTA),
            ("fungal_gtf", DUAL_FUNGUS_GTF),
            ("reference_cache_dir", DUAL_REFERENCE_CACHE_DIR or "not used"),
            ("mapping_branches", "host-only,competitive-host-pathogen"),
            ("host_only_star_index", DUAL_HOST_ONLY_STAR_INDEX_DIR),
            ("competitive_star_index", DUAL_COMPETITIVE_STAR_INDEX_DIR),
            ("featurecounts_strand", DUAL_FEATURECOUNTS_STRAND),
            ("host_rrna_exclusion_list", input.rrna),
            ("outlier_policy", "all samples retained; no automatic replacement"),
            ("kraken_database", DUAL_KRAKEN_DB),
            ("kraken_sources", ",".join(DUAL_KRAKEN_SOURCES)),
            ("kraken_focus_families", ",".join(DUAL_KRAKEN_FOCUS_FAMILIES)),
            ("bracken_enabled", str(DUAL_RUN_BRACKEN)),
        ]
        with open(output[0], "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["parameter", "value"])
            writer.writerows(rows)


rule dual_fastp:
    input:
        r1=lambda wildcards: dual_raw_read(wildcards, 1),
        r2=lambda wildcards: dual_raw_read(wildcards, 2)
    output:
        r1=DUAL_DIR + "/01-trimmed-fastp/{sample}_1.trimmed.fastq.gz",
        r2=DUAL_DIR + "/01-trimmed-fastp/{sample}_2.trimmed.fastq.gz",
        json=DUAL_DIR + "/01-trimmed-fastp/qc/{sample}.fastp.json",
        html=DUAL_DIR + "/01-trimmed-fastp/qc/{sample}.fastp.html"
    params:
        fastp_module=FASTP_MODULE,
        conda_module=FASTP_CONDA_MODULE,
        conda_env=FASTP_CONDA_ENV,
        cmd=FASTP_CMD
    threads: 16
    shell:
        r"""
        set -euo pipefail
        if command -v module >/dev/null 2>&1; then
          module purge
          if [ -n "{params.fastp_module}" ]; then module load {params.fastp_module}; fi
          if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
        fi
        if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
        mkdir -p $(dirname {output.r1}) $(dirname {output.json})
        {params.cmd} \
          --thread {threads} \
          --in1 {input.r1} \
          --in2 {input.r2} \
          --out1 {output.r1} \
          --out2 {output.r2} \
          --trim_front1 2 \
          --trim_front2 2 \
          --detect_adapter_for_pe \
          --length_required 50 \
          --json {output.json} \
          --html {output.html}
        """


rule dual_build_combined_reference:
    input:
        host_fasta=DUAL_HOST_FASTA,
        host_gtf=DUAL_HOST_GTF,
        fungus_fasta=DUAL_FUNGUS_FASTA,
        fungus_gtf=DUAL_FUNGUS_GTF
    output:
        combined_fasta=DUAL_COMBINED_FASTA,
        combined_gtf=DUAL_COMBINED_GTF,
        host_gtf=DUAL_HOST_PREFIXED_GTF,
        fungus_gtf=DUAL_FUNGUS_PREFIXED_GTF,
        sequence_map=DUAL_SEQUENCE_MAP,
        gene_map=DUAL_GENE_MAP
    shell:
        r"""
        python3 scripts/build_combined_reference.py \
          --host-fasta {input.host_fasta} \
          --host-gtf {input.host_gtf} \
          --fungus-fasta {input.fungus_fasta} \
          --fungus-gtf {input.fungus_gtf} \
          --host-prefix '{DUAL_HOST_PREFIX}' \
          --fungus-prefix '{DUAL_FUNGUS_PREFIX}' \
          --combined-fasta {output.combined_fasta} \
          --combined-gtf {output.combined_gtf} \
          --host-prefixed-gtf {output.host_gtf} \
          --fungus-prefixed-gtf {output.fungus_gtf} \
          --sequence-map {output.sequence_map} \
          --gene-map {output.gene_map}
        """


rule dual_host_only_star_index:
    input:
        fasta=DUAL_HOST_FASTA,
        gtf=DUAL_HOST_GTF
    output:
        index=directory(DUAL_HOST_ONLY_STAR_INDEX_DIR)
    threads: 20
    shell:
        r"""
        set -euo pipefail
        module purge
        module load star-2.7.10b-gcc-12.1.0
        mkdir -p {output.index}
        STAR \
          --runMode genomeGenerate \
          --runThreadN {threads} \
          --genomeDir {output.index} \
          --genomeFastaFiles {input.fasta} \
          --sjdbGTFfile {input.gtf} \
          --alignIntronMax 2500000 \
          --sjdbOverhang 149
        """


rule dual_competitive_star_index:
    input:
        fasta=DUAL_COMBINED_FASTA,
        gtf=DUAL_COMBINED_GTF
    output:
        index=directory(DUAL_COMPETITIVE_STAR_INDEX_DIR)
    threads: 20
    shell:
        r"""
        set -euo pipefail
        module purge
        module load star-2.7.10b-gcc-12.1.0
        mkdir -p {output.index}
        STAR \
          --runMode genomeGenerate \
          --runThreadN {threads} \
          --genomeDir {output.index} \
          --genomeFastaFiles {input.fasta} \
          --sjdbGTFfile {input.gtf} \
          --alignIntronMax 2500000 \
          --sjdbOverhang 149
        """


rule dual_host_only_star_align:
    input:
        index=DUAL_HOST_ONLY_STAR_INDEX_DIR,
        r1=lambda wildcards: dual_trimmed_read(wildcards, 1),
        r2=lambda wildcards: dual_trimmed_read(wildcards, 2)
    output:
        bam=DUAL_DIR + "/02-host-only-star/{sample}_Aligned.sortedByCoord.out.bam",
        csi=DUAL_DIR + "/02-host-only-star/{sample}_Aligned.sortedByCoord.out.bam.csi",
        log=DUAL_DIR + "/02-host-only-star/{sample}_Log.final.out",
        sj=DUAL_DIR + "/02-host-only-star/{sample}_SJ.out.tab"
    params:
        prefix=DUAL_DIR + "/02-host-only-star/{sample}_"
    threads: 16
    shell:
        r"""
        set -euo pipefail
        module purge
        module load star-2.7.10b-gcc-12.1.0
        mkdir -p $(dirname {output.bam})
        STAR \
          --runThreadN {threads} \
          --genomeDir {input.index} \
          --genomeLoad NoSharedMemory \
          --limitBAMsortRAM 64000000000 \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} LB:Stranded_Total_RNA_RiboZero PL:Illumina PU:NovaSeqXPlus \
          --outSAMattributes NH HI AS NM MD \
          --outSAMmapqUnique {DUAL_STAR_UNIQUE_MAPQ} \
          --outSAMunmapped Within KeepPairs \
          --outFilterMultimapNmax 20 \
          --twopassMode Basic \
          --alignIntronMax 2500000 \
          --readFilesCommand zcat \
          --readFilesIn {input.r1} {input.r2} \
          --outFileNamePrefix {params.prefix}
        module purge
        module load {SAMTOOLS_MODULE}
        samtools index -c -@ {threads} {output.bam}
        """


rule dual_competitive_star_align:
    input:
        index=DUAL_COMPETITIVE_STAR_INDEX_DIR,
        r1=lambda wildcards: dual_trimmed_read(wildcards, 1),
        r2=lambda wildcards: dual_trimmed_read(wildcards, 2)
    output:
        bam=DUAL_DIR + "/02-competitive-star/{sample}_Aligned.sortedByCoord.out.bam",
        csi=DUAL_DIR + "/02-competitive-star/{sample}_Aligned.sortedByCoord.out.bam.csi",
        log=DUAL_DIR + "/02-competitive-star/{sample}_Log.final.out",
        sj=DUAL_DIR + "/02-competitive-star/{sample}_SJ.out.tab"
    params:
        prefix=DUAL_DIR + "/02-competitive-star/{sample}_"
    threads: 16
    shell:
        r"""
        set -euo pipefail
        module purge
        module load star-2.7.10b-gcc-12.1.0
        mkdir -p $(dirname {output.bam})
        STAR \
          --runThreadN {threads} \
          --genomeDir {input.index} \
          --genomeLoad NoSharedMemory \
          --limitBAMsortRAM 64000000000 \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} LB:Stranded_Total_RNA_RiboZero PL:Illumina PU:NovaSeqXPlus \
          --outSAMattributes NH HI AS NM MD \
          --outSAMmapqUnique {DUAL_STAR_UNIQUE_MAPQ} \
          --outSAMunmapped Within KeepPairs \
          --outFilterMultimapNmax 20 \
          --twopassMode Basic \
          --alignIntronMax 2500000 \
          --readFilesCommand zcat \
          --readFilesIn {input.r1} {input.r2} \
          --outFileNamePrefix {params.prefix}
        module purge
        module load {SAMTOOLS_MODULE}
        samtools index -c -@ {threads} {output.bam}
        """


rule dual_featurecounts_host_only_transcript_exon:
    input:
        bam=DUAL_DIR + "/02-host-only-star/{sample}_Aligned.sortedByCoord.out.bam",
        gtf=DUAL_HOST_GTF
    output:
        counts=DUAL_DIR + "/03-featurecounts/host-only-transcript-exon/{sample}.featureCounts.txt",
        summary=DUAL_DIR + "/03-featurecounts/host-only-transcript-exon/{sample}.featureCounts.txt.summary"
    threads: 12
    shell:
        r"""
        set -euo pipefail
        module purge
        module load mamba
        source activate subread
        mkdir -p $(dirname {output.counts})
        featureCounts \
          -p --countReadPairs -B -C \
          -s {DUAL_FEATURECOUNTS_STRAND} \
          -t transcript,exon \
          -g gene_id \
          --extraAttributes gene_name \
          --primary \
          -Q 10 \
          -T {threads} \
          -a {input.gtf} \
          -o {output.counts} \
          {input.bam}
        """


rule dual_featurecounts_competitive_host_transcript_exon:
    input:
        bam=DUAL_DIR + "/02-competitive-star/{sample}_Aligned.sortedByCoord.out.bam",
        gtf=DUAL_HOST_PREFIXED_GTF
    output:
        counts=DUAL_DIR + "/03-featurecounts/competitive-host-transcript-exon/{sample}.featureCounts.txt",
        summary=DUAL_DIR + "/03-featurecounts/competitive-host-transcript-exon/{sample}.featureCounts.txt.summary"
    threads: 12
    shell:
        r"""
        set -euo pipefail
        module purge
        module load mamba
        source activate subread
        mkdir -p $(dirname {output.counts})
        featureCounts \
          -p --countReadPairs -B -C \
          -s {DUAL_FEATURECOUNTS_STRAND} \
          -t transcript,exon \
          -g gene_id \
          --extraAttributes gene_name \
          --primary \
          -Q 10 \
          -T {threads} \
          -a {input.gtf} \
          -o {output.counts} \
          {input.bam}
        """


rule dual_featurecounts_fungus_transcript_exon:
    input:
        bam=DUAL_DIR + "/02-competitive-star/{sample}_Aligned.sortedByCoord.out.bam",
        gtf=DUAL_FUNGUS_PREFIXED_GTF
    output:
        counts=DUAL_DIR + "/03-featurecounts/competitive-fungus-transcript-exon/{sample}.featureCounts.txt",
        summary=DUAL_DIR + "/03-featurecounts/competitive-fungus-transcript-exon/{sample}.featureCounts.txt.summary"
    threads: 12
    shell:
        r"""
        set -euo pipefail
        module purge
        module load mamba
        source activate subread
        mkdir -p $(dirname {output.counts})
        featureCounts \
          -p --countReadPairs -B -C \
          -s {DUAL_FEATURECOUNTS_STRAND} \
          -t transcript,exon \
          -g gene_id \
          --extraAttributes gene_name \
          --primary \
          -Q 10 \
          -T {threads} \
          -a {input.gtf} \
          -o {output.counts} \
          {input.bam}
        """


rule dual_featurecounts_host_only_exon:
    input:
        bam=DUAL_DIR + "/02-host-only-star/{sample}_Aligned.sortedByCoord.out.bam",
        gtf=DUAL_HOST_GTF
    output:
        counts=DUAL_DIR + "/03-featurecounts/host-only-exon/{sample}.featureCounts.txt",
        summary=DUAL_DIR + "/03-featurecounts/host-only-exon/{sample}.featureCounts.txt.summary"
    threads: 12
    shell:
        r"""
        set -euo pipefail
        module purge
        module load mamba
        source activate subread
        mkdir -p $(dirname {output.counts})
        featureCounts \
          -p --countReadPairs -B -C \
          -s {DUAL_FEATURECOUNTS_STRAND} \
          -t exon \
          -g gene_id \
          --extraAttributes gene_name \
          --primary \
          -Q 10 \
          -T {threads} \
          -a {input.gtf} \
          -o {output.counts} \
          {input.bam}
        """


rule dual_featurecounts_competitive_host_exon:
    input:
        bam=DUAL_DIR + "/02-competitive-star/{sample}_Aligned.sortedByCoord.out.bam",
        gtf=DUAL_HOST_PREFIXED_GTF
    output:
        counts=DUAL_DIR + "/03-featurecounts/competitive-host-exon/{sample}.featureCounts.txt",
        summary=DUAL_DIR + "/03-featurecounts/competitive-host-exon/{sample}.featureCounts.txt.summary"
    threads: 12
    shell:
        r"""
        set -euo pipefail
        module purge
        module load mamba
        source activate subread
        mkdir -p $(dirname {output.counts})
        featureCounts \
          -p --countReadPairs -B -C \
          -s {DUAL_FEATURECOUNTS_STRAND} \
          -t exon \
          -g gene_id \
          --extraAttributes gene_name \
          --primary \
          -Q 10 \
          -T {threads} \
          -a {input.gtf} \
          -o {output.counts} \
          {input.bam}
        """


rule dual_merge_host_only_transcript_exon_counts:
    input:
        dual_featurecounts(
            DUAL_DIR
            + "/03-featurecounts/host-only-transcript-exon/{sample}.featureCounts.txt"
        )
    output:
        DUAL_HOST_ONLY_PRIMARY_MATRIX
    params:
        sample_args=" ".join(
            "--sample-file "
            + sample
            + "="
            + DUAL_DIR
            + f"/03-featurecounts/host-only-transcript-exon/{sample}.featureCounts.txt"
            for sample in DUAL_SAMPLES
        )
    shell:
        r"""
        python3 scripts/merge_featurecounts.py \
          {params.sample_args} \
          --output {output}
        """


rule dual_merge_competitive_host_transcript_exon_counts:
    input:
        dual_featurecounts(
            DUAL_DIR
            + "/03-featurecounts/competitive-host-transcript-exon/{sample}.featureCounts.txt"
        )
    output:
        DUAL_COMPETITIVE_HOST_PRIMARY_MATRIX
    params:
        sample_args=" ".join(
            "--sample-file "
            + sample
            + "="
            + DUAL_DIR
            + f"/03-featurecounts/competitive-host-transcript-exon/{sample}.featureCounts.txt"
            for sample in DUAL_SAMPLES
        )
    shell:
        r"""
        python3 scripts/merge_featurecounts.py \
          {params.sample_args} \
          --strip-prefix '{DUAL_HOST_PREFIX}' \
          --output {output}
        """


rule dual_merge_fungus_transcript_exon_counts:
    input:
        dual_featurecounts(
            DUAL_DIR
            + "/03-featurecounts/competitive-fungus-transcript-exon/{sample}.featureCounts.txt"
        )
    output:
        DUAL_FUNGUS_PRIMARY_MATRIX
    params:
        sample_args=" ".join(
            "--sample-file "
            + sample
            + "="
            + DUAL_DIR
            + f"/03-featurecounts/competitive-fungus-transcript-exon/{sample}.featureCounts.txt"
            for sample in DUAL_SAMPLES
        )
    shell:
        r"""
        python3 scripts/merge_featurecounts.py \
          {params.sample_args} \
          --strip-prefix '{DUAL_FUNGUS_PREFIX}' \
          --output {output}
        """


rule dual_merge_host_only_exon_counts:
    input:
        dual_featurecounts(
            DUAL_DIR + "/03-featurecounts/host-only-exon/{sample}.featureCounts.txt"
        )
    output:
        DUAL_HOST_ONLY_EXON_MATRIX
    params:
        sample_args=" ".join(
            "--sample-file "
            + sample
            + "="
            + DUAL_DIR
            + f"/03-featurecounts/host-only-exon/{sample}.featureCounts.txt"
            for sample in DUAL_SAMPLES
        )
    shell:
        r"""
        python3 scripts/merge_featurecounts.py \
          {params.sample_args} \
          --output {output}
        """


rule dual_merge_competitive_host_exon_counts:
    input:
        dual_featurecounts(
            DUAL_DIR + "/03-featurecounts/competitive-host-exon/{sample}.featureCounts.txt"
        )
    output:
        DUAL_COMPETITIVE_HOST_EXON_MATRIX
    params:
        sample_args=" ".join(
            "--sample-file "
            + sample
            + "="
            + DUAL_DIR
            + f"/03-featurecounts/competitive-host-exon/{sample}.featureCounts.txt"
            for sample in DUAL_SAMPLES
        )
    shell:
        r"""
        python3 scripts/merge_featurecounts.py \
          {params.sample_args} \
          --strip-prefix '{DUAL_HOST_PREFIX}' \
          --output {output}
        """


rule dual_export_host_unmapped:
    input:
        bam=DUAL_DIR + "/02-host-only-star/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        r1=DUAL_DIR + "/02-host-only-star/unmapped/{sample}_host_unmapped_R1.fastq.gz",
        r2=DUAL_DIR + "/02-host-only-star/unmapped/{sample}_host_unmapped_R2.fastq.gz",
        stats=DUAL_DIR + "/02-host-only-star/unmapped/{sample}.unmapped_pairs.tsv"
    params:
        tmp=DUAL_DIR + "/02-host-only-star/unmapped/{sample}.collate"
    threads: 4
    shell:
        r"""
        set -euo pipefail
        module purge
        module load {SAMTOOLS_MODULE}
        mkdir -p $(dirname {output.r1})
        samtools collate -@ {threads} -u -O {input.bam} {params.tmp} \
          | samtools fastq -@ {threads} -f 12 -F 256 -n \
              -1 >(gzip -c > {output.r1}) \
              -2 >(gzip -c > {output.r2}) \
              -0 /dev/null -s /dev/null -
        pairs=$(gzip -cd {output.r1} | awk 'END {{print int(NR/4)}}')
        printf "sample_id\thost_unmapped_pairs\n{wildcards.sample}\t%s\n" "$pairs" \
          > {output.stats}
        """


rule dual_export_competitive_unmapped:
    input:
        bam=DUAL_DIR + "/02-competitive-star/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        r1=DUAL_DIR + "/02-competitive-star/unmapped/{sample}_competitive_unmapped_R1.fastq.gz",
        r2=DUAL_DIR + "/02-competitive-star/unmapped/{sample}_competitive_unmapped_R2.fastq.gz",
        stats=DUAL_DIR + "/02-competitive-star/unmapped/{sample}.unmapped_pairs.tsv"
    params:
        tmp=DUAL_DIR + "/02-competitive-star/unmapped/{sample}.collate"
    threads: 4
    shell:
        r"""
        set -euo pipefail
        module purge
        module load {SAMTOOLS_MODULE}
        mkdir -p $(dirname {output.r1})
        samtools collate -@ {threads} -u -O {input.bam} {params.tmp} \
          | samtools fastq -@ {threads} -f 12 -F 256 -n \
              -1 >(gzip -c > {output.r1}) \
              -2 >(gzip -c > {output.r2}) \
              -0 /dev/null -s /dev/null -
        pairs=$(gzip -cd {output.r1} | awk 'END {{print int(NR/4)}}')
        printf "sample_id\tcompetitive_unmapped_pairs\n{wildcards.sample}\t%s\n" "$pairs" \
          > {output.stats}
        """


rule dual_host_only_mapping_summary_sample:
    input:
        log=DUAL_DIR + "/02-host-only-star/{sample}_Log.final.out"
    output:
        DUAL_DIR + "/05-mapping-comparison/host-only/{sample}.star_mapping.tsv"
    shell:
        r"""
        python3 scripts/summarize_star_log.py \
          --sample {wildcards.sample} \
          --branch host-only \
          --input {input.log} \
          --output {output}
        """


rule dual_merge_host_only_mapping_summary:
    input:
        expand(
            DUAL_DIR + "/05-mapping-comparison/host-only/{sample}.star_mapping.tsv",
            sample=DUAL_SAMPLES,
        )
    output:
        DUAL_HOST_ONLY_MAPPING_QC
    params:
        input_args=" ".join(
            "--input "
            + DUAL_DIR
            + f"/05-mapping-comparison/host-only/{sample}.star_mapping.tsv"
            for sample in DUAL_SAMPLES
        )
    shell:
        r"""
        python3 scripts/merge_tsv_rows.py \
          {params.input_args} \
          --output {output}
        """


rule dual_competitive_mapping_summary_sample:
    input:
        bam=DUAL_DIR + "/02-competitive-star/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        DUAL_DIR + "/05-mapping-comparison/competitive/{sample}.mapping.tsv"
    threads: 4
    shell:
        r"""
        set -euo pipefail
        module purge
        module load {SAMTOOLS_MODULE}
        mkdir -p $(dirname {output})
        samtools view -@ {threads} -F 2308 {input.bam} \
          | python3 scripts/summarize_competitive_mapping.py \
              --sample {wildcards.sample} \
              --host-prefix '{DUAL_HOST_PREFIX}' \
              --fungus-prefix '{DUAL_FUNGUS_PREFIX}' \
              --output {output}
        """


rule dual_merge_competitive_mapping_summary:
    input:
        expand(
            DUAL_DIR + "/05-mapping-comparison/competitive/{sample}.mapping.tsv",
            sample=DUAL_SAMPLES,
        )
    output:
        DUAL_COMPETITIVE_MAPPING_QC
    params:
        input_args=" ".join(
            "--input "
            + DUAL_DIR
            + f"/05-mapping-comparison/competitive/{sample}.mapping.tsv"
            for sample in DUAL_SAMPLES
        )
    shell:
        r"""
        python3 scripts/merge_tsv_rows.py \
          {params.input_args} \
          --output {output}
        """


rule dual_kraken2_classify:
    input:
        r1=lambda wildcards: dual_kraken_read(wildcards, 1),
        r2=lambda wildcards: dual_kraken_read(wildcards, 2)
    output:
        report=DUAL_DIR + "/08-taxonomy/{source}/{sample}.kraken2.report",
        classified=DUAL_DIR + "/08-taxonomy/{source}/{sample}.kraken2.output.gz"
    params:
        db=DUAL_KRAKEN_DB,
        confidence=DUAL_KRAKEN_CONFIDENCE,
        min_hit_groups=DUAL_KRAKEN_MIN_HIT_GROUPS,
        conda_module=METATX_CONDA_MODULE,
        conda_env=METATX_CONDA_ENV,
        cmd=KRAKEN2_CMD
    threads: 16
    wildcard_constraints:
        source="all_trimmed|host_unmapped|competitive_unmapped"
    shell:
        r"""
        set -euo pipefail
        if command -v module >/dev/null 2>&1; then
          module purge
          if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
        fi
        if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
        mkdir -p $(dirname {output.report})
        pairs=$(gzip -cd {input.r1} | awk 'END {{print int(NR/4)}}')
        if [ "${{pairs:-0}}" -eq 0 ]; then
          printf "100.00\t0\t0\tU\t0\tunclassified\n" > {output.report}
          gzip -c </dev/null > {output.classified}
        else
          {params.cmd} \
            --db {params.db} \
            --threads {threads} \
            --paired \
            --gzip-compressed \
            --confidence {params.confidence} \
            --minimum-hit-groups {params.min_hit_groups} \
            --report {output.report} \
            {input.r1} {input.r2} \
            | gzip -c > {output.classified}
        fi
        """


rule dual_bracken_family:
    input:
        report=DUAL_DIR + "/08-taxonomy/{source}/{sample}.kraken2.report"
    output:
        abundance=DUAL_DIR + "/08-taxonomy/{source}/{sample}.bracken.family.tsv",
        report=DUAL_DIR + "/08-taxonomy/{source}/{sample}.bracken.family.report"
    params:
        db=DUAL_KRAKEN_DB,
        cmd=DUAL_BRACKEN_CMD,
        read_length=DUAL_BRACKEN_READ_LENGTH,
        threshold=DUAL_BRACKEN_THRESHOLD,
        conda_module=METATX_CONDA_MODULE,
        conda_env=METATX_CONDA_ENV
    threads: 4
    shell:
        r"""
        set -euo pipefail
        if command -v module >/dev/null 2>&1; then
          module purge
          if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
        fi
        if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
        mkdir -p $(dirname {output.abundance})
        total_reads=$(awk '{{ total += $3 }} END {{ print total + 0 }}' {input.report})
        if [ "${{total_reads}}" -eq 0 ]; then
          printf 'name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n' \
            > {output.abundance}
          cp {input.report} {output.report}
        else
          {params.cmd} \
            -d {params.db} \
            -i {input.report} \
            -o {output.abundance} \
            -w {output.report} \
            -r {params.read_length} \
            -l F \
            -t {params.threshold}
        fi
        """


rule dual_taxonomy_manifest:
    input:
        kraken=dual_kraken_reports(),
        bracken=dual_bracken_outputs()
    output:
        DUAL_TAXONOMY_MANIFEST
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w", newline="") as handle:
            fieldnames = [
                "sample_id",
                "source",
                "treatment",
                "diet",
                "abundance_format",
                "abundance_file",
                "kraken_report",
            ]
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for source in DUAL_KRAKEN_SOURCES:
                for record in DUAL_RECORDS:
                    sample = record["sample_id"]
                    kraken_report = (
                        DUAL_DIR
                        + f"/08-taxonomy/{source}/{sample}.kraken2.report"
                    )
                    if DUAL_RUN_BRACKEN:
                        abundance_format = "bracken"
                        abundance_file = (
                            DUAL_DIR
                            + f"/08-taxonomy/{source}/{sample}.bracken.family.tsv"
                        )
                    else:
                        abundance_format = "kraken"
                        abundance_file = kraken_report
                    writer.writerow(
                        {
                            "sample_id": sample,
                            "source": source,
                            "treatment": record["treatment"],
                            "diet": record["diet"],
                            "abundance_format": abundance_format,
                            "abundance_file": abundance_file,
                            "kraken_report": kraken_report,
                        }
                    )


rule dual_prepare_local_transfer_bundle:
    input:
        results=dual_local_transfer_files(),
        metadata=DUAL_SAMPLE_TABLE,
        rrna=DUAL_RRNA_LIST
    output:
        manifest=DUAL_TRANSFER_MANIFEST,
        bundle=DUAL_TRANSFER_BUNDLE
    run:
        output_dir = Path(output.manifest).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        run_root = Path(DUAL_DIR).resolve()
        transfer_entries = []

        def add_entry(source, archive_path):
            source = Path(str(source)).resolve()
            digest = hashlib.sha256()
            with source.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
            transfer_entries.append(
                {
                    "source_file": str(source),
                    "archive_path": archive_path,
                    "size_bytes": source.stat().st_size,
                    "sha256": digest.hexdigest(),
                }
            )

        for source in input.results:
            source_path = Path(str(source)).resolve()
            try:
                archive_path = str(source_path.relative_to(run_root))
            except ValueError:
                # Add-on runs may read immutable maps from the prior reference
                # cache. Preserve them in the bundle without pretending they
                # were generated in the new sample-output directory.
                reference_root = Path(DUAL_REFERENCE_DIR).resolve()
                archive_path = "00-reference-cache/" + str(
                    source_path.relative_to(reference_root)
                )
            add_entry(source_path, archive_path)
        add_entry(input.metadata, "inputs/fatbody_dual_rnaseq_samples.tsv")
        add_entry(input.rrna, "inputs/gregaria_rrna_list.txt")

        with open(output.manifest, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["source_file", "archive_path", "size_bytes", "sha256"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(transfer_entries)

        with tarfile.open(output.bundle, "w:gz") as archive:
            for entry in transfer_entries:
                archive.add(entry["source_file"], arcname=entry["archive_path"])
            archive.add(output.manifest, arcname="transfer_manifest.tsv")


rule dual_plot_family_composition:
    input:
        manifest=DUAL_TAXONOMY_MANIFEST
    output:
        DUAL_FAMILY_PLOT_OUTPUTS
    params:
        outdir=DUAL_DIR + "/08-taxonomy/figures",
        top=DUAL_KRAKEN_TOP_FAMILIES,
        focus=",".join(DUAL_KRAKEN_FOCUS_FAMILIES),
        r=DUAL_R_CMD
    shell:
        r"""
        {params.r} scripts/plot_kraken_family_composition.R \
          --manifest {input.manifest} \
          --output-dir {params.outdir} \
          --top-families {params.top} \
          --focus-families '{params.focus}'
        """


rule dual_host_only_deseq2:
    input:
        counts=DUAL_HOST_ONLY_PRIMARY_MATRIX,
        metadata=DUAL_SAMPLE_TABLE,
        rrna=DUAL_RRNA_LIST
    output:
        DUAL_HOST_ONLY_DE_OUTPUTS
    params:
        outdir=DUAL_HOST_ONLY_DE_DIR,
        branch="host-only",
        alpha=DUAL_HOST_DE_ALPHA,
        lfc=DUAL_HOST_DE_LFC,
        min_count=DUAL_HOST_DE_MIN_COUNT,
        min_samples=DUAL_HOST_DE_MIN_SAMPLES,
        r=DUAL_R_CMD
    shell:
        r"""
        {params.r} scripts/run_host_deseq2.R \
          --counts {input.counts} \
          --metadata {input.metadata} \
          --rrna-list {input.rrna} \
          --mapping-branch {params.branch} \
          --output-dir {params.outdir} \
          --alpha {params.alpha} \
          --lfc-threshold {params.lfc} \
          --min-count {params.min_count} \
          --min-samples {params.min_samples}
        """


rule dual_competitive_host_deseq2:
    input:
        counts=DUAL_COMPETITIVE_HOST_PRIMARY_MATRIX,
        metadata=DUAL_SAMPLE_TABLE,
        rrna=DUAL_RRNA_LIST
    output:
        DUAL_COMPETITIVE_HOST_DE_OUTPUTS
    params:
        outdir=DUAL_COMPETITIVE_HOST_DE_DIR,
        branch="competitive-host",
        alpha=DUAL_HOST_DE_ALPHA,
        lfc=DUAL_HOST_DE_LFC,
        min_count=DUAL_HOST_DE_MIN_COUNT,
        min_samples=DUAL_HOST_DE_MIN_SAMPLES,
        r=DUAL_R_CMD
    shell:
        r"""
        {params.r} scripts/run_host_deseq2.R \
          --counts {input.counts} \
          --metadata {input.metadata} \
          --rrna-list {input.rrna} \
          --mapping-branch {params.branch} \
          --output-dir {params.outdir} \
          --alpha {params.alpha} \
          --lfc-threshold {params.lfc} \
          --min-count {params.min_count} \
          --min-samples {params.min_samples}
        """


rule dual_compare_host_count_matrices:
    input:
        host_only=DUAL_HOST_ONLY_PRIMARY_MATRIX,
        competitive_host=DUAL_COMPETITIVE_HOST_PRIMARY_MATRIX,
        competitive_fungus=DUAL_FUNGUS_PRIMARY_MATRIX
    output:
        sample=DUAL_COUNT_RECONCILIATION_SAMPLE,
        gene=DUAL_COUNT_RECONCILIATION_GENE
    shell:
        r"""
        python3 scripts/compare_host_count_matrices.py \
          --host-only {input.host_only} \
          --competitive-host {input.competitive_host} \
          --competitive-fungus {input.competitive_fungus} \
          --sample-output {output.sample} \
          --gene-output {output.gene}
        """


rule dual_compare_host_deseq2:
    input:
        host_only=DUAL_HOST_ONLY_DE_DIR + "/host_deseq2_all_contrasts.csv",
        competitive=DUAL_COMPETITIVE_HOST_DE_DIR + "/host_deseq2_all_contrasts.csv"
    output:
        rows=DUAL_DE_RECONCILIATION,
        summary=DUAL_DE_RECONCILIATION_SUMMARY
    shell:
        r"""
        python3 scripts/compare_host_deseq2_results.py \
          --host-only {input.host_only} \
          --competitive {input.competitive} \
          --output {output.rows} \
          --summary {output.summary}
        """


rule dual_identify_unplaced_host_degs:
    input:
        results=DUAL_HOST_ONLY_DE_DIR + "/host_deseq2_all_contrasts.csv",
        reconciliation=DUAL_DE_RECONCILIATION,
        gtf=DUAL_HOST_GTF
    output:
        genes=DUAL_UNPLACED_DEG_GENES,
        scaffolds=DUAL_UNPLACED_SCAFFOLD_LIST
    shell:
        r"""
        python3 scripts/identify_unplaced_deg_scaffolds.py \
          --deseq-results {input.results} \
          --mapping-reconciliation {input.reconciliation} \
          --gtf {input.gtf} \
          --output-genes {output.genes} \
          --output-scaffolds {output.scaffolds}
        """


rule dual_extract_candidate_unplaced_scaffolds:
    input:
        fasta=DUAL_HOST_FASTA,
        sequence_list=DUAL_UNPLACED_SCAFFOLD_LIST
    output:
        fasta=DUAL_UNPLACED_FASTA,
        metrics=DUAL_UNPLACED_METRICS
    shell:
        r"""
        python3 scripts/extract_unplaced_scaffolds.py \
          --fasta {input.fasta} \
          --sequence-list {input.sequence_list} \
          --output-fasta {output.fasta} \
          --output-table {output.metrics}
        """


rule dual_kraken_candidate_unplaced_scaffolds:
    input:
        fasta=DUAL_UNPLACED_FASTA
    output:
        report=DUAL_UNPLACED_KRAKEN_REPORT,
        classified=DUAL_UNPLACED_KRAKEN_OUTPUT
    params:
        db=DUAL_KRAKEN_DB,
        confidence=DUAL_KRAKEN_CONFIDENCE,
        min_hit_groups=DUAL_KRAKEN_MIN_HIT_GROUPS,
        conda_module=METATX_CONDA_MODULE,
        conda_env=METATX_CONDA_ENV,
        cmd=KRAKEN2_CMD
    threads: 16
    shell:
        r"""
        set -euo pipefail
        if command -v module >/dev/null 2>&1; then
          module purge
          if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
        fi
        if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
        mkdir -p $(dirname {output.report})
        if grep -q '^>' {input.fasta}; then
          {params.cmd} \
            --db {params.db} \
            --threads {threads} \
            --confidence {params.confidence} \
            --minimum-hit-groups {params.min_hit_groups} \
            --report {output.report} \
            --output {output.classified} \
            {input.fasta}
        else
          printf "100.00\t0\t0\tU\t0\tunclassified\n" > {output.report}
          : > {output.classified}
        fi
        """


rule dual_summarize_candidate_unplaced_scaffolds:
    input:
        metrics=DUAL_UNPLACED_METRICS,
        classified=DUAL_UNPLACED_KRAKEN_OUTPUT,
        report=DUAL_UNPLACED_KRAKEN_REPORT
    output:
        DUAL_UNPLACED_TAXONOMY
    shell:
        r"""
        python3 scripts/summarize_scaffold_kraken.py \
          --metrics {input.metrics} \
          --kraken-output {input.classified} \
          --kraken-report {input.report} \
          --output {output}
        """


rule dual_build_host_origin_evidence_catalogue:
    input:
        reconciliation=DUAL_DE_RECONCILIATION,
        gene_counts=DUAL_COUNT_RECONCILIATION_GENE,
        unplaced_genes=DUAL_UNPLACED_DEG_GENES,
        scaffold_taxonomy=DUAL_UNPLACED_TAXONOMY
    output:
        DUAL_FINAL_ORIGIN_EVIDENCE
    shell:
        r"""
        python3 scripts/build_host_origin_evidence_catalogue.py \
          --deg-reconciliation {input.reconciliation} \
          --gene-count-comparison {input.gene_counts} \
          --unplaced-genes {input.unplaced_genes} \
          --scaffold-taxonomy {input.scaffold_taxonomy} \
          --output {output}
        """


rule dual_fungal_deseq2:
    input:
        counts=DUAL_FUNGUS_PRIMARY_MATRIX,
        metadata=DUAL_SAMPLE_TABLE
    output:
        DUAL_FUNGAL_DE_OUTPUTS
    params:
        outdir=DUAL_FUNGAL_DE_DIR,
        min_assigned=DUAL_FUNGAL_DE_MIN_ASSIGNED,
        min_genes=DUAL_FUNGAL_DE_MIN_GENES,
        min_replicates=DUAL_FUNGAL_DE_MIN_REPLICATES,
        alpha=DUAL_FUNGAL_DE_ALPHA,
        lfc=DUAL_FUNGAL_DE_LFC,
        r=DUAL_R_CMD
    shell:
        r"""
        {params.r} scripts/run_fungal_deseq2.R \
          --counts {input.counts} \
          --metadata {input.metadata} \
          --output-dir {params.outdir} \
          --min-assigned {params.min_assigned} \
          --min-genes {params.min_genes} \
          --min-replicates {params.min_replicates} \
          --alpha {params.alpha} \
          --lfc-threshold {params.lfc}
        """
