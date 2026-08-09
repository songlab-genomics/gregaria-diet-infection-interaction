#!/usr/bin/env bash
set -euo pipefail

# Validate references, all sample FASTQs, and Kraken/Bracken resources before
# submitting the R-free HPC preprocessing workflow.

PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
DUAL_SAMPLE_TABLE="${DUAL_SAMPLE_TABLE:-config/fatbody_dual_rnaseq_samples.tsv}"
DUAL_REFERENCE_CACHE_DIR="${DUAL_REFERENCE_CACHE_DIR:-}"
DUAL_INPUT_MODE="${DUAL_INPUT_MODE:-raw}"
DUAL_RAW_DIR="${DUAL_RAW_DIR:-/data/songlab/sequencing_data/RNAseq/mehreen}"
DUAL_TRIMMED_DIR="${DUAL_TRIMMED_DIR:-${PROJECT_ROOT}/input/trimmed-fastp}"
REFDIR="${REFDIR:-/data/songlab/maeva/gregaria-timecourse/reference}"
DUAL_HOST_FASTA="${DUAL_HOST_FASTA:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.fna}"
DUAL_HOST_GTF="${DUAL_HOST_GTF:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.gtf}"
DUAL_FUNGUS_REFERENCE_DIR="${DUAL_FUNGUS_REFERENCE_DIR:-${PROJECT_ROOT}/reference/metarhizium/GCF_000187425.2}"
DUAL_FUNGUS_FASTA="${DUAL_FUNGUS_FASTA:-${DUAL_FUNGUS_REFERENCE_DIR}/GCF_000187425.2_genomic.fna}"
DUAL_FUNGUS_GTF="${DUAL_FUNGUS_GTF:-${DUAL_FUNGUS_REFERENCE_DIR}/GCF_000187425.2_genomic.gtf}"
DUAL_KRAKEN_DB="${DUAL_KRAKEN_DB:-/scratch/mtecher/kraken2/pluspf_20260626}"
DUAL_RUN_BRACKEN="${DUAL_RUN_BRACKEN:-true}"
DUAL_BRACKEN_READ_LENGTH="${DUAL_BRACKEN_READ_LENGTH:-150}"
DUAL_RRNA_LIST="${DUAL_RRNA_LIST:-../data/excluded_loci/gregaria_rrna_list.txt}"
DUAL_REQUIRE_R="${DUAL_REQUIRE_R:-false}"
FASTP_CMD="${FASTP_CMD:-fastp}"
METATX_CONDA_MODULE="${METATX_CONDA_MODULE:-mamba}"
METATX_CONDA_ENV="${METATX_CONDA_ENV:-/scratch/mtecher/conda-envs/metatranscriptomics-qc}"
KRAKEN2_CMD="${KRAKEN2_CMD:-kraken2}"
DUAL_BRACKEN_CMD="${DUAL_BRACKEN_CMD:-bracken}"

if [ ! -s "${DUAL_SAMPLE_TABLE}" ]; then
  echo "ERROR: missing sample table: ${DUAL_SAMPLE_TABLE}" >&2
  exit 2
fi

if [ -n "${DUAL_REFERENCE_CACHE_DIR}" ]; then
  for cached_reference in \
    "${DUAL_REFERENCE_CACHE_DIR}/combined_host_metarhizium.fna" \
    "${DUAL_REFERENCE_CACHE_DIR}/combined_host_metarhizium.gtf" \
    "${DUAL_REFERENCE_CACHE_DIR}/host.prefixed.gtf" \
    "${DUAL_REFERENCE_CACHE_DIR}/metarhizium.prefixed.gtf" \
    "${DUAL_REFERENCE_CACHE_DIR}/STAR-host-only/Genome" \
    "${DUAL_REFERENCE_CACHE_DIR}/STAR-competitive/Genome"; do
    if [ ! -s "${cached_reference}" ]; then
      echo "ERROR: cached reference product missing: ${cached_reference}" >&2
      exit 15
    fi
  done
fi

for reference in \
  "${DUAL_HOST_FASTA}" \
  "${DUAL_HOST_GTF}" \
  "${DUAL_FUNGUS_FASTA}" \
  "${DUAL_FUNGUS_GTF}"; do
  if [ ! -s "${reference}" ]; then
    echo "ERROR: missing reference file: ${reference}" >&2
    exit 3
  fi
done

if [ ! -s "${DUAL_RRNA_LIST}" ]; then
  echo "ERROR: missing host rRNA exclusion list: ${DUAL_RRNA_LIST}" >&2
  exit 11
fi

missing_reads=0
sample_count=0
while IFS=$'\t' read -r sample_id raw_prefix tissue treatment diet; do
  if [ "${sample_id}" = "sample_id" ]; then
    continue
  fi
  sample_count=$((sample_count + 1))
  for mate in 1 2; do
    if [ "${DUAL_INPUT_MODE}" = "trimmed" ]; then
      read_path="${DUAL_TRIMMED_DIR}/${sample_id}_${mate}.trimmed.fastq.gz"
    elif [ "${DUAL_INPUT_MODE}" = "raw" ]; then
      read_path="${DUAL_RAW_DIR}/${raw_prefix}_MERGE_${mate}.fq.gz"
    else
      echo "ERROR: DUAL_INPUT_MODE must be trimmed or raw." >&2
      exit 4
    fi
    if [ ! -s "${read_path}" ]; then
      echo "MISSING READ: ${read_path}" >&2
      missing_reads=$((missing_reads + 1))
    fi
  done
done < "${DUAL_SAMPLE_TABLE}"

if [ "${missing_reads}" -ne 0 ]; then
  echo "ERROR: ${missing_reads} FASTQ files are missing." >&2
  exit 5
fi

for database_file in hash.k2d opts.k2d taxo.k2d; do
  if [ ! -s "${DUAL_KRAKEN_DB}/${database_file}" ]; then
    echo "ERROR: Kraken2 database file missing: ${DUAL_KRAKEN_DB}/${database_file}" >&2
    exit 6
  fi
done

if [[ "${DUAL_RUN_BRACKEN,,}" =~ ^(1|true|yes|y|on)$ ]]; then
  bracken_distribution="${DUAL_KRAKEN_DB}/database${DUAL_BRACKEN_READ_LENGTH}mers.kmer_distrib"
  if [ ! -s "${bracken_distribution}" ]; then
    echo "ERROR: Bracken database file missing: ${bracken_distribution}" >&2
    echo "Build Bracken for read length ${DUAL_BRACKEN_READ_LENGTH}, or set DUAL_RUN_BRACKEN=false." >&2
    exit 7
  fi
fi

# Test the same environment activation performed inside Kraken/Bracken workers.
# Database files alone are insufficient if the classification commands are
# absent from the worker environment.
(
  if command -v module >/dev/null 2>&1; then
    module purge
    if [ -n "${METATX_CONDA_MODULE}" ]; then
      module load "${METATX_CONDA_MODULE}"
    fi
  fi

  if [ -n "${METATX_CONDA_ENV}" ]; then
    if ! source activate "${METATX_CONDA_ENV}"; then
      echo "ERROR: cannot activate Kraken/Bracken environment: ${METATX_CONDA_ENV}" >&2
      echo "Install it with: sbatch scripts/install_metatranscriptomics_qc_env.slurm" >&2
      exit 12
    fi
  fi

  if ! command -v "${KRAKEN2_CMD}" >/dev/null 2>&1; then
    echo "ERROR: Kraken2 command unavailable in ${METATX_CONDA_ENV}: ${KRAKEN2_CMD}" >&2
    exit 13
  fi

  if [[ "${DUAL_RUN_BRACKEN,,}" =~ ^(1|true|yes|y|on)$ ]] \
    && ! command -v "${DUAL_BRACKEN_CMD}" >/dev/null 2>&1; then
    echo "ERROR: Bracken command unavailable in ${METATX_CONDA_ENV}: ${DUAL_BRACKEN_CMD}" >&2
    exit 14
  fi

  echo "  Kraken2 executable: $(command -v "${KRAKEN2_CMD}")"
  "${KRAKEN2_CMD}" --version | sed -n '1p'
  if [[ "${DUAL_RUN_BRACKEN,,}" =~ ^(1|true|yes|y|on)$ ]]; then
    echo "  Bracken executable: $(command -v "${DUAL_BRACKEN_CMD}")"
  fi
)

if ! command -v snakemake >/dev/null 2>&1; then
  echo "ERROR: snakemake is unavailable in the active controller environment." >&2
  exit 8
fi

if [ "${DUAL_INPUT_MODE}" = "raw" ] && ! command -v "${FASTP_CMD}" >/dev/null 2>&1; then
  echo "ERROR: ${FASTP_CMD} is unavailable in the active controller environment." >&2
  exit 10
fi

if [[ "${DUAL_REQUIRE_R,,}" =~ ^(1|true|yes|y|on)$ ]]; then
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript is unavailable but DUAL_REQUIRE_R=true." >&2
    exit 9
  fi
  Rscript -e 'required <- c("DESeq2","dplyr","ggplot2","readr","tidyr","tibble","scales"); missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop("Missing R packages: ", paste(missing, collapse=", "))'
fi

echo "Dual RNA-seq preflight passed."
echo "  Samples: ${sample_count}"
echo "  Input mode: ${DUAL_INPUT_MODE}"
echo "  Raw reads: ${DUAL_RAW_DIR}"
echo "  External trimmed reads (trimmed mode only): ${DUAL_TRIMMED_DIR}"
echo "  Host reference: ${DUAL_HOST_FASTA}"
echo "  Reference cache: ${DUAL_REFERENCE_CACHE_DIR:-not used}"
echo "  Fungal reference: ${DUAL_FUNGUS_FASTA}"
echo "  Kraken database: ${DUAL_KRAKEN_DB}"
echo "  Kraken/Bracken environment: ${METATX_CONDA_ENV}"
echo "  Bracken enabled: ${DUAL_RUN_BRACKEN}"
echo "  Host rRNA exclusion list: ${DUAL_RRNA_LIST}"
echo "  Mapping branches: host-only and competitive host-pathogen"
echo "  R analysis required on HPC: ${DUAL_REQUIRE_R}"
