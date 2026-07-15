### =================================================================
### CONFIG
### =================================================================

import os


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
METATX_CONDA_ENV = env_str("METATX_CONDA_ENV", "metatranscriptomics-qc")

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
