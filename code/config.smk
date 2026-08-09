### =================================================================
### CONFIG
### =================================================================

import os
from datetime import datetime


def env_str(name, default):
    """Let HPC/local launch scripts override paths without editing this file."""
    return os.environ.get(name, default)


def env_list(name, default):
    value = os.environ.get(name)
    if value is None or value.strip() == "":
        return default
    return [item.strip() for item in value.split(",") if item.strip()]


def env_bool(name, default):
    value = os.environ.get(name)
    if value is None or value.strip() == "":
        return default
    return value.strip().lower() in {"1", "true", "yes", "y", "on"}


# --------- Species / tissues ----------
SPECIES_TO_PROCESS = env_list("LOCUST_SPECIES", ["gregaria"])
TISSUES = env_list("LOCUST_TISSUES", ["AGY"])

# --------- Paths ----------
DATADIR = env_str("DATADIR", f"/data/songlab/sequencing_data/RNAseq/timecourse/{SPECIES_TO_PROCESS[0]}")
PROJECTDIR = env_str("PROJECTDIR", f"/data/songlab/maeva/{SPECIES_TO_PROCESS[0]}-timecourse")
REFDIR  = env_str("REFDIR", f"{PROJECTDIR}/reference")

# Existing processed RNA-seq folders on the HPC live in /data/songlab/maeva.
# STAR BAMs and count files are read from here when using completed alignments.
WORKDIR = env_str("WORKDIR", f"{PROJECTDIR}/data")

# Sample discovery can use raw reads, completed STAR BAMs, or auto fallback.
# Use "star_bam" when calling variants from existing alignments only.
SAMPLE_DISCOVERY = env_str("SAMPLE_DISCOVERY", "star_bam")

# Scratch is used for new, heavy intermediates so existing processed folders are
# not overwritten by exploratory variant-calling runs.
SCRATCHDIR = env_str("SCRATCHDIR", f"/scratch/mtecher/{SPECIES_TO_PROCESS[0]}-timecourse")

KAIJUDIR  = "/data/songlab/maeva/kaijudb"
KRAKENDIR = "/data/songlab/maeva/kraken2"

# --------- Reference mapping ----------
species_to_genome = {
    "gregaria":    "GCF_023897955.1_iqSchGreg1.2",
    "cancellata":  "GCF_023864275.1_iqSchCanc2.1",
    "piceifrons":  "GCF_021461385.2_iqSchPice1.1",
    "americana":   "GCF_021461395.2_iqSchAmer2.1",
    "cubense":     "GCF_023864345.2_iqSchSeri2.2",
    "nitens":      "GCF_023898315.1_iqSchNite1.1",
}

species_to_label = {
    "gregaria":    "Schistocerca gregaria",
    "cancellata":  "Schistocerca cancellata",
    "piceifrons":  "Schistocerca piceifrons",
    "americana":   "Schistocerca americana",
    "cubense":     "Schistocerca serialis cubense",
    "nitens":      "Schistocerca nitens",
}

# Convenience (single species mode like your current pipeline)
SPECIES = SPECIES_TO_PROCESS[0]
GENOME_ID = species_to_genome[SPECIES]

REF_GENOME = f"{REFDIR}/{GENOME_ID}_genomic.fna"
REF_GTF    = f"{REFDIR}/{GENOME_ID}_genomic.gtf"
REF_TRANSCRIPTOME= f"{REFDIR}/{GENOME_ID}_rna_from_genomic.fna"
REF_PROTEIN = f"{REFDIR}/{GENOME_ID}_protein.faa"

STAR_INDEX_DIR = f"{REFDIR}/index/{SPECIES}/STAR"
SALMON_DIR       = f"{REFDIR}/index/{SPECIES}/SALMON"
SALMON_INDEX_DIR = f"{SALMON_DIR}/{SPECIES}_index"

# --------- RNA-seq variant calling ----------
# Keep this False for ordinary expression-only runs. Set to True, or target
# rule rnaseq_variants_all directly, when you want transcriptome-derived SNPs.
RUN_RNASEQ_VARIANTS = False

# Change this label for every independent run so old callsets stay untouched.
RNASEQ_VARIANT_RUN_ID = env_str("RNASEQ_VARIANT_RUN_ID", "rnaseq_variants_v1")
RNASEQ_VARIANT_DIR = env_str("RNASEQ_VARIANT_DIR", f"{SCRATCHDIR}/rnaseq-variants/{RNASEQ_VARIANT_RUN_ID}")

# GATK needs a sequence dictionary next to the genome FASTA. This name follows
# the reference basename without copying or editing the genome file itself.
REF_DICT = REF_GENOME.rsplit(".", 1)[0] + ".dict"

# Tool modules are kept configurable because HPC module names vary by cluster.
# Leave a value empty ("") if the tool is already available in your environment.
SAMTOOLS_MODULE = env_str("SAMTOOLS_MODULE", "samtools-1.21-gcc-12.1.0")
GATK_MODULE = env_str("GATK_MODULE", "")

# fastp is used in the first trimming rule. Keep these configurable so the
# pilot can use the same controller environment, a separate conda environment,
# a module-only installation, or a full executable path on the cluster.
FASTP_MODULE = env_str("FASTP_MODULE", "")
FASTP_CONDA_MODULE = env_str("FASTP_CONDA_MODULE", "mamba")
FASTP_CONDA_ENV = env_str("FASTP_CONDA_ENV", "myENV")
FASTP_CMD = env_str("FASTP_CMD", "fastp")

# Command used after loading GATK_MODULE. If your cluster has no GATK module,
# set GATK_MODULE = "" and point this to a conda executable or wrapper script.
GATK_CMD = env_str("GATK_CMD", "gatk")

# GATK is available from this mamba/conda environment on the HPC. The workflow
# activates it inside each batch job because SLURM is launched with --export=NONE.
GATK_CONDA_MODULE = env_str("GATK_CONDA_MODULE", "mamba")
GATK_CONDA_ENV = env_str("GATK_CONDA_ENV", "piceifrons-popgenomics")

# --------- eggNOG-mapper functional annotation ----------
# This annotation run is deliberately separate from the RNA-seq mapping rules.
# It reads the NCBI protein FASTA and writes a new eggNOG output folder.
RUN_EGGNOG_ANNOTATION = False

# Change this label for every independent eggNOG run to keep older annotations.
EGGNOG_RUN_ID = env_str("EGGNOG_RUN_ID", f"eggnog_{SPECIES}_desertlocustr_v1")
EGGNOG_DIR = env_str("EGGNOG_DIR", f"{SCRATCHDIR}/eggnog/{EGGNOG_RUN_ID}")
EGGNOG_PREFIX = env_str("EGGNOG_PREFIX", SPECIES)

# Input protein FASTA. If your HPC reference folder uses a different name,
# change only this variable before launching the eggNOG-only Snakefile.
EGGNOG_PROTEIN_FASTA = env_str("EGGNOG_PROTEIN_FASTA", REF_PROTEIN)

# Protein-to-LOCID map used after eggNOG to convert protein IDs to GeneID/LOCID.
EGGNOG_PROTEIN2GENE = env_str("EGGNOG_PROTEIN2GENE", "../data/list/allspecies_protein2geneid.tsv")
EGGNOG_SPECIES_LABEL = env_str("EGGNOG_SPECIES_LABEL", species_to_label.get(SPECIES, SPECIES))

# Optional GO-name source used only to label and split eggNOG GO terms into
# BP/CC/MF tables. If absent, the LOCID universe is still written with GO IDs.
EGGNOG_GO_NAME_TABLE = env_str(
    "EGGNOG_GO_NAME_TABLE",
    f"../data/list/GO_Annotations/blast2go_{SPECIES}_custom.txt",
)

# eggNOG database directory. Run the download target once if this folder has not
# already been prepared on the HPC.
EGGNOG_DATA_DIR = env_str("EGGNOG_DATA_DIR", "/scratch/mtecher/eggnog-mapper-data")

# Tool loading is configurable because HPC module/env names vary. Leave module
# blank when emapper.py is provided by the conda/mamba environment.
EGGNOG_MODULE = env_str("EGGNOG_MODULE", "")
EGGNOG_CONDA_MODULE = env_str("EGGNOG_CONDA_MODULE", "mamba")
EGGNOG_CONDA_ENV = env_str("EGGNOG_CONDA_ENV", "eggnog-mapper")
EGGNOG_CMD = env_str("EGGNOG_CMD", "emapper.py")
EGGNOG_DOWNLOAD_CMD = env_str("EGGNOG_DOWNLOAD_CMD", "download_eggnog_data.py")

# Conservative defaults: keep electronic GO terms out and require closer
# orthology evidence. Set EGGNOG_TAX_SCOPE = "50557" to restrict to Insecta if
# your installed eggNOG-mapper version supports --tax_scope.
EGGNOG_GO_EVIDENCE = env_str("EGGNOG_GO_EVIDENCE", "non-electronic")
EGGNOG_TARGET_ORTHOLOGS = env_str("EGGNOG_TARGET_ORTHOLOGS", "one2one,many2one")
EGGNOG_TAX_SCOPE = env_str("EGGNOG_TAX_SCOPE", "")

# --------- Metatranscriptomics and cross-sample QC ----------
# This branch reuses completed STAR BAMs and trimmed FASTQs. It writes a new
# scratch output folder and does not modify read counts or host transcriptomes.
METATX_RUN_ID = env_str("METATX_RUN_ID", "metatx_qc_v1")
METATX_SPECIES = env_list("METATX_SPECIES", ["gregaria", "piceifrons", "americana", "cubense"])

# Leave empty to discover all tissues with completed STAR BAMs. Use comma lists
# such as METATX_TISSUES=AGY,ALB,MTG to restrict the run.
METATX_TISSUES = env_list("METATX_TISSUES", [])
METATX_PROJECT_BASE = env_str("METATX_PROJECT_BASE", "/data/songlab/maeva")
METATX_SCRATCH_BASE = env_str("METATX_SCRATCH_BASE", "/scratch/mtecher/locust-time-course-RNAseq")
METATX_DIR = env_str("METATX_DIR", f"{METATX_SCRATCH_BASE}/metatranscriptomics/{METATX_RUN_ID}")

# Host-unmapped reads are more specific for microbes/viruses/fungi than all
# trimmed reads. The QC tables still report total and host-mapped read numbers.
METATX_CLASSIFY_SOURCE = env_str("METATX_CLASSIFY_SOURCE", "host_unmapped")
METATX_RUN_KRAKEN2 = env_bool("METATX_RUN_KRAKEN2", True)
METATX_RUN_KAIJU = env_bool("METATX_RUN_KAIJU", False)

METATX_CONDA_MODULE = env_str("METATX_CONDA_MODULE", "mamba")
METATX_CONDA_ENV = env_str(
    "METATX_CONDA_ENV",
    "/scratch/mtecher/conda-envs/metatranscriptomics-qc",
)

KRAKEN2_DB = env_str("KRAKEN2_DB", KRAKENDIR)
KRAKEN2_MODULE = env_str("KRAKEN2_MODULE", "")
KRAKEN2_CMD = env_str("KRAKEN2_CMD", "kraken2")

KAIJU_MODULE = env_str("KAIJU_MODULE", "")
KAIJU_CMD = env_str("KAIJU_CMD", "kaiju")
KAIJU2TABLE_CMD = env_str("KAIJU2TABLE_CMD", "kaiju2table")
KAIJU_DB_FMI = env_str("KAIJU_DB_FMI", f"{KAIJUDIR}/kaiju_db.fmi")
KAIJU_DB_NODES = env_str("KAIJU_DB_NODES", f"{KAIJUDIR}/nodes.dmp")
KAIJU_DB_NAMES = env_str("KAIJU_DB_NAMES", f"{KAIJUDIR}/names.dmp")
KAIJU_TAXON_LEVEL = env_str("KAIJU_TAXON_LEVEL", "species")

# --------- Fat-body host-pathogen dual RNA-seq ----------
# This branch is isolated from the historical host-only workflow. Every launch
# should set a new DUAL_RUN_ID or DUAL_DIR so previous results remain intact.
DUAL_RUN_ID = env_str(
    "DUAL_RUN_ID",
    f"host_pathogen_dual_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
)
DUAL_DIR = env_str(
    "DUAL_DIR",
    f"/scratch/mtecher/gregaria-diet-infection-interaction/output/runs/{DUAL_RUN_ID}",
)
DUAL_SAMPLE_TABLE = env_str(
    "DUAL_SAMPLE_TABLE",
    "config/fatbody_dual_rnaseq_samples.tsv",
)

# Start from the protected Mehreen raw reads by default. Every raw-mode run
# creates its own fastp outputs inside the new timestamped run directory.
DUAL_INPUT_MODE = env_str("DUAL_INPUT_MODE", "raw")
DUAL_RAW_DIR = env_str(
    "DUAL_RAW_DIR",
    "/data/songlab/sequencing_data/RNAseq/mehreen",
)
DUAL_TRIMMED_DIR = env_str(
    "DUAL_TRIMMED_DIR",
    "/scratch/mtecher/gregaria-diet-infection-interaction/input/trimmed-fastp",
)

# Host and fungal references. DWR2009/ARSEF 10343 has no public assembly found
# under that isolate name, so ARSEF 23 is recorded as the mapping proxy.
DUAL_HOST_FASTA = env_str("DUAL_HOST_FASTA", REF_GENOME)
DUAL_HOST_GTF = env_str("DUAL_HOST_GTF", REF_GTF)
DUAL_FUNGUS_ACCESSION = env_str("DUAL_FUNGUS_ACCESSION", "GCF_000187425.2")
DUAL_FUNGUS_ISOLATE = env_str("DUAL_FUNGUS_ISOLATE", "ARSEF 23")
DUAL_EXPERIMENTAL_ISOLATE = env_str("DUAL_EXPERIMENTAL_ISOLATE", "DWR2009 / ARSEF 10343")
DUAL_FUNGUS_REFERENCE_DIR = env_str(
    "DUAL_FUNGUS_REFERENCE_DIR",
    "/scratch/mtecher/gregaria-diet-infection-interaction/reference/metarhizium/GCF_000187425.2",
)
DUAL_FUNGUS_FASTA = env_str(
    "DUAL_FUNGUS_FASTA",
    f"{DUAL_FUNGUS_REFERENCE_DIR}/GCF_000187425.2_genomic.fna",
)
DUAL_FUNGUS_GTF = env_str(
    "DUAL_FUNGUS_GTF",
    f"{DUAL_FUNGUS_REFERENCE_DIR}/GCF_000187425.2_genomic.gtf",
)

# A one-library add-on run can reuse the immutable reference products and STAR
# indexes from the completed 44-library run. Sample outputs still go to the new
# DUAL_DIR, so the historical run is never modified.
DUAL_REFERENCE_CACHE_DIR = env_str("DUAL_REFERENCE_CACHE_DIR", "")
DUAL_REFERENCE_DIR = (
    DUAL_REFERENCE_CACHE_DIR
    if DUAL_REFERENCE_CACHE_DIR
    else f"{DUAL_DIR}/00-reference"
)
DUAL_HOST_ONLY_STAR_INDEX_DIR = f"{DUAL_REFERENCE_DIR}/STAR-host-only"
DUAL_COMPETITIVE_STAR_INDEX_DIR = f"{DUAL_REFERENCE_DIR}/STAR-competitive"
DUAL_HOST_PREFIX = env_str("DUAL_HOST_PREFIX", "HOST__")
DUAL_FUNGUS_PREFIX = env_str("DUAL_FUNGUS_PREFIX", "MR__")

# Illumina Stranded Total RNA libraries are normally reverse stranded. The
# preflight asks users to confirm this with infer_experiment.py; 2 is the
# featureCounts setting for reverse-stranded paired-end data.
DUAL_FEATURECOUNTS_STRAND = int(env_str("DUAL_FEATURECOUNTS_STRAND", "2"))
DUAL_STAR_UNIQUE_MAPQ = int(env_str("DUAL_STAR_UNIQUE_MAPQ", "60"))

# Kraken is run on all trimmed reads, reads unmapped to the host-only reference,
# and reads unmapped to the competitive host-pathogen reference.
DUAL_KRAKEN_SOURCES = env_list(
    "DUAL_KRAKEN_SOURCES",
    ["all_trimmed", "host_unmapped", "competitive_unmapped"],
)
DUAL_KRAKEN_DB = env_str(
    "DUAL_KRAKEN_DB",
    "/scratch/mtecher/kraken2/pluspf_20260626",
)
DUAL_KRAKEN_CONFIDENCE = env_str("DUAL_KRAKEN_CONFIDENCE", "0.05")
DUAL_KRAKEN_MIN_HIT_GROUPS = int(env_str("DUAL_KRAKEN_MIN_HIT_GROUPS", "2"))
DUAL_KRAKEN_TOP_FAMILIES = int(env_str("DUAL_KRAKEN_TOP_FAMILIES", "25"))
DUAL_KRAKEN_FOCUS_FAMILIES = env_list(
    "DUAL_KRAKEN_FOCUS_FAMILIES",
    ["Clavicipitaceae"],
)
DUAL_RUN_BRACKEN = env_bool("DUAL_RUN_BRACKEN", True)
DUAL_BRACKEN_CMD = env_str("DUAL_BRACKEN_CMD", "bracken")
DUAL_BRACKEN_READ_LENGTH = int(env_str("DUAL_BRACKEN_READ_LENGTH", "150"))
DUAL_BRACKEN_THRESHOLD = int(env_str("DUAL_BRACKEN_THRESHOLD", "10"))
DUAL_R_CMD = env_str("DUAL_R_CMD", "Rscript")
DUAL_RRNA_LIST = env_str(
    "DUAL_RRNA_LIST",
    "../data/excluded_loci/gregaria_rrna_list.txt",
)

# Both host count matrices receive the same all-samples DESeq2 analysis. No
# samples are replaced or excluded; the explicit rRNA list is removed first.
DUAL_HOST_DE_ALPHA = float(env_str("DUAL_HOST_DE_ALPHA", "0.05"))
DUAL_HOST_DE_LFC = float(env_str("DUAL_HOST_DE_LFC", "1"))
DUAL_HOST_DE_MIN_COUNT = int(env_str("DUAL_HOST_DE_MIN_COUNT", "10"))
DUAL_HOST_DE_MIN_SAMPLES = int(env_str("DUAL_HOST_DE_MIN_SAMPLES", "3"))

# Fungal DE is attempted only for infected libraries with adequate fungal
# counts. Thresholds remain configurable because the pilot determines whether
# fat body contains enough fungal RNA for formal inference.
DUAL_FUNGAL_DE_MIN_ASSIGNED = int(env_str("DUAL_FUNGAL_DE_MIN_ASSIGNED", "50000"))
DUAL_FUNGAL_DE_MIN_GENES = int(env_str("DUAL_FUNGAL_DE_MIN_GENES", "1000"))
DUAL_FUNGAL_DE_MIN_REPLICATES = int(env_str("DUAL_FUNGAL_DE_MIN_REPLICATES", "3"))
DUAL_FUNGAL_DE_ALPHA = float(env_str("DUAL_FUNGAL_DE_ALPHA", "0.05"))
DUAL_FUNGAL_DE_LFC = float(env_str("DUAL_FUNGAL_DE_LFC", "1"))

# --------- Independent unplaced-scaffold origin audit ----------
# This target starts from a frozen local DEG scaffold set. It is intentionally
# independent of fastp, STAR, featureCounts, and R so those completed stages
# cannot be resubmitted while investigating scaffold provenance.
ORIGIN_RUN_ID = env_str(
    "ORIGIN_RUN_ID",
    f"scaffold_origin_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
)
ORIGIN_DIR = env_str(
    "ORIGIN_DIR",
    f"/scratch/mtecher/gregaria-diet-infection-interaction/output/runs/{ORIGIN_RUN_ID}",
)
ORIGIN_CANDIDATE_LIST = env_str(
    "ORIGIN_CANDIDATE_LIST",
    "config/scaffold_origin/unified_unplaced_deg_scaffolds_20260729.txt",
)
ORIGIN_CANDIDATE_GENES = env_str(
    "ORIGIN_CANDIDATE_GENES",
    "config/scaffold_origin/unified_unplaced_deg_genes_20260729.tsv",
)
ORIGIN_CANDIDATE_PROVENANCE = env_str(
    "ORIGIN_CANDIDATE_PROVENANCE",
    "config/scaffold_origin/unified_unplaced_deg_scaffold_provenance_20260729.tsv",
)
ORIGIN_NCBI_FCS_REPORT = env_str(
    "ORIGIN_NCBI_FCS_REPORT",
    "config/scaffold_origin/GCF_023897955.1_iqSchGreg1.2_fcs_report_20230323.txt",
)
ORIGIN_HOST_FASTA = env_str("ORIGIN_HOST_FASTA", DUAL_HOST_FASTA)
ORIGIN_HOST_ANNOTATION = env_str("ORIGIN_HOST_ANNOTATION", DUAL_HOST_GTF)
ORIGIN_EXPECTED_SCAFFOLDS = int(env_str("ORIGIN_EXPECTED_SCAFFOLDS", "77"))
ORIGIN_KRAKEN_DB = env_str("ORIGIN_KRAKEN_DB", DUAL_KRAKEN_DB)
ORIGIN_KRAKEN_CONFIDENCE = env_str(
    "ORIGIN_KRAKEN_CONFIDENCE", DUAL_KRAKEN_CONFIDENCE
)
ORIGIN_KRAKEN_MIN_HIT_GROUPS = int(
    env_str("ORIGIN_KRAKEN_MIN_HIT_GROUPS", str(DUAL_KRAKEN_MIN_HIT_GROUPS))
)
ORIGIN_CONDA_MODULE = env_str("ORIGIN_CONDA_MODULE", METATX_CONDA_MODULE)
ORIGIN_CONDA_ENV = env_str("ORIGIN_CONDA_ENV", METATX_CONDA_ENV)
ORIGIN_KRAKEN_CMD = env_str("ORIGIN_KRAKEN_CMD", KRAKEN2_CMD)
ORIGIN_RUN_CROSS_SPECIES = env_bool("ORIGIN_RUN_CROSS_SPECIES", True)
ORIGIN_CROSS_SPECIES_REFERENCE_ITEMS = env_list(
    "ORIGIN_CROSS_SPECIES_REFERENCES",
    [
        "piceifrons=/data/songlab/maeva/piceifrons-timecourse/reference/GCF_021461385.2_iqSchPice1.1_genomic.fna",
        "americana=/data/songlab/maeva/americana-timecourse/reference/GCF_021461395.2_iqSchAmer2.1_genomic.fna",
        "serialis_cubense=/data/songlab/maeva/cubense-timecourse/reference/GCF_023864345.2_iqSchSeri2.2_genomic.fna",
    ],
)
ORIGIN_CROSS_SPECIES_REFERENCES = {}
for origin_reference_item in ORIGIN_CROSS_SPECIES_REFERENCE_ITEMS:
    if "=" not in origin_reference_item:
        raise ValueError(
            "ORIGIN_CROSS_SPECIES_REFERENCES entries must use species=/path/to/genome.fna"
        )
    origin_species, origin_reference = origin_reference_item.split("=", 1)
    ORIGIN_CROSS_SPECIES_REFERENCES[origin_species] = origin_reference
ORIGIN_MINIMAP2_CMD = env_str("ORIGIN_MINIMAP2_CMD", "minimap2")

# --------- Focused BLAST/DIAMOND origin audit ----------
# This is a separate run starting from the frozen unresolved DEG set for the
# infected diet 33 versus 50 contrast. It never invokes RNA-seq preprocessing.
HOMOLOGY_RUN_ID = env_str(
    "HOMOLOGY_RUN_ID",
    f"scaffold_homology_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
)
HOMOLOGY_DIR = env_str(
    "HOMOLOGY_DIR",
    f"/scratch/mtecher/gregaria-diet-infection-interaction/output/runs/{HOMOLOGY_RUN_ID}",
)
HOMOLOGY_CANDIDATE_SCAFFOLDS = env_str(
    "HOMOLOGY_CANDIDATE_SCAFFOLDS",
    "config/scaffold_homology/infected_diet33_vs50_unresolved_scaffolds_20260731.txt",
)
HOMOLOGY_CANDIDATE_GENES = env_str(
    "HOMOLOGY_CANDIDATE_GENES",
    "config/scaffold_homology/infected_diet33_vs50_unresolved_genes_20260731.tsv",
)
HOMOLOGY_CANDIDATE_PROVENANCE = env_str(
    "HOMOLOGY_CANDIDATE_PROVENANCE",
    "config/scaffold_homology/infected_diet33_vs50_unresolved_provenance_20260731.tsv",
)
HOMOLOGY_HOST_FASTA = env_str("HOMOLOGY_HOST_FASTA", ORIGIN_HOST_FASTA)
HOMOLOGY_HOST_ANNOTATION = env_str(
    "HOMOLOGY_HOST_ANNOTATION", ORIGIN_HOST_ANNOTATION
)
HOMOLOGY_HOST_PROTEIN_FASTA = env_str(
    "HOMOLOGY_HOST_PROTEIN_FASTA",
    "/data/songlab/maeva/gregaria-timecourse/reference/GCF_023897955.1_iqSchGreg1.2_protein.faa",
)
HOMOLOGY_BLAST_DB_ROOT = env_str(
    "HOMOLOGY_BLAST_DB_ROOT",
    "/scratch/mtecher/scaffold_homology_databases/ncbi_blastdb",
)
HOMOLOGY_CORE_NT_DB = env_str(
    "HOMOLOGY_CORE_NT_DB", f"{HOMOLOGY_BLAST_DB_ROOT}/core_nt"
)
HOMOLOGY_NR_CLUSTER_SEQ_DB = env_str(
    "HOMOLOGY_NR_CLUSTER_SEQ_DB", f"{HOMOLOGY_BLAST_DB_ROOT}/nr_cluster_seq"
)
HOMOLOGY_DATABASE_READY = env_str(
    "HOMOLOGY_DATABASE_READY",
    f"{HOMOLOGY_BLAST_DB_ROOT}/download_complete.tsv",
)
HOMOLOGY_TAXONOMY_NODES = env_str(
    "HOMOLOGY_TAXONOMY_NODES", f"{ORIGIN_KRAKEN_DB}/nodes.dmp"
)
HOMOLOGY_TAXONOMY_NAMES = env_str(
    "HOMOLOGY_TAXONOMY_NAMES", f"{ORIGIN_KRAKEN_DB}/names.dmp"
)
HOMOLOGY_CONDA_MODULE = env_str("HOMOLOGY_CONDA_MODULE", "mamba")
HOMOLOGY_CONDA_ENV = env_str(
    "HOMOLOGY_CONDA_ENV",
    "/scratch/mtecher/conda-envs/scaffold-homology",
)
HOMOLOGY_BLASTN_CMD = env_str("HOMOLOGY_BLASTN_CMD", "blastn")
HOMOLOGY_DIAMOND_CMD = env_str("HOMOLOGY_DIAMOND_CMD", "diamond")
HOMOLOGY_SELF_TAXID = int(env_str("HOMOLOGY_SELF_TAXID", "7010"))
HOMOLOGY_EXPECTED_SCAFFOLDS = int(
    env_str("HOMOLOGY_EXPECTED_SCAFFOLDS", "59")
)
HOMOLOGY_EXPECTED_GENES = int(env_str("HOMOLOGY_EXPECTED_GENES", "948"))
HOMOLOGY_MAX_TARGET_SEQS = int(env_str("HOMOLOGY_MAX_TARGET_SEQS", "25"))
HOMOLOGY_TOP_FRACTION = float(env_str("HOMOLOGY_TOP_FRACTION", "0.05"))
