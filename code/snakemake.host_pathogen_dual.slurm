#!/usr/bin/env bash
#SBATCH --job-name=fatbody_dual_smk
#SBATCH --output=slurm-fatbody-dual-%j.out
#SBATCH --error=slurm-fatbody-dual-%j.err
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=6-23:59:00
#SBATCH --export=ALL

set -euo pipefail

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -s "${SLURM_SUBMIT_DIR}/cluster.json" ]; then
  cd "${SLURM_SUBMIT_DIR}"
else
  cd "$(dirname "${BASH_SOURCE[0]}")"
fi

TARGET="${1:-host_pathogen_dual_hpc_preprocessing}"
JOBS="${2:-200}"
SNAKEMAKE_ENV="${3:-myENV}"
SNAKEFILE="${4:-Snakefile.host_pathogen_dual_rnaseq}"
LATENCY_WAIT="${LATENCY_WAIT:-120}"
CLUSTER_CONFIG="${CLUSTER_CONFIG:-$(pwd)/cluster.json}"

PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
if [ -z "${DUAL_RUN_ID:-}" ]; then
  echo "ERROR: DUAL_RUN_ID must be exported before submission." >&2
  echo "For a new run, use: export DUAL_RUN_ID=host_pathogen_dual_\$(date +%Y%m%d_%H%M%S)" >&2
  echo "For a resume, export the exact existing run ID and DUAL_RESUME_EXISTING=true." >&2
  exit 2
fi
DUAL_DIR="${DUAL_DIR:-${PROJECT_ROOT}/output/runs/${DUAL_RUN_ID}}"
DUAL_SAMPLE_TABLE="${DUAL_SAMPLE_TABLE:-config/fatbody_dual_rnaseq_samples.tsv}"
DUAL_INPUT_MODE="${DUAL_INPUT_MODE:-raw}"
DUAL_RAW_DIR="${DUAL_RAW_DIR:-/data/songlab/sequencing_data/RNAseq/mehreen}"
DUAL_TRIMMED_DIR="${DUAL_TRIMMED_DIR:-${PROJECT_ROOT}/input/trimmed-fastp}"
REFDIR="${REFDIR:-/data/songlab/maeva/gregaria-timecourse/reference}"
PROJECTDIR="${PROJECTDIR:-$(dirname "${REFDIR}")}"
SCRATCHDIR="${SCRATCHDIR:-${PROJECT_ROOT}/scratch/${DUAL_RUN_ID}}"

DUAL_HOST_FASTA="${DUAL_HOST_FASTA:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.fna}"
DUAL_HOST_GTF="${DUAL_HOST_GTF:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.gtf}"
DUAL_FUNGUS_REFERENCE_DIR="${DUAL_FUNGUS_REFERENCE_DIR:-${PROJECT_ROOT}/reference/metarhizium/GCF_000187425.2}"
DUAL_FUNGUS_FASTA="${DUAL_FUNGUS_FASTA:-${DUAL_FUNGUS_REFERENCE_DIR}/GCF_000187425.2_genomic.fna}"
DUAL_FUNGUS_GTF="${DUAL_FUNGUS_GTF:-${DUAL_FUNGUS_REFERENCE_DIR}/GCF_000187425.2_genomic.gtf}"
DUAL_KRAKEN_DB="${DUAL_KRAKEN_DB:-/scratch/mtecher/kraken2/pluspf_20260626}"
DUAL_RUN_BRACKEN="${DUAL_RUN_BRACKEN:-true}"
DUAL_KRAKEN_SOURCES="${DUAL_KRAKEN_SOURCES:-all_trimmed,host_unmapped,competitive_unmapped}"
DUAL_FEATURECOUNTS_STRAND="${DUAL_FEATURECOUNTS_STRAND:-2}"
DUAL_RRNA_LIST="${DUAL_RRNA_LIST:-../data/excluded_loci/gregaria_rrna_list.txt}"
DUAL_RESUME_EXISTING="${DUAL_RESUME_EXISTING:-false}"

CLUSTER_LOG_DIR="${CLUSTER_LOG_DIR:-${PROJECT_ROOT}/logs/slurm/${DUAL_RUN_ID}}"

if [ "${TARGET}" = "host_pathogen_dual_hpc_preprocessing" ]; then
  DUAL_REQUIRE_R=false
else
  DUAL_REQUIRE_R="${DUAL_REQUIRE_R:-true}"
fi

if [ -e "${DUAL_DIR}" ]; then
  if [[ "${DUAL_RESUME_EXISTING,,}" =~ ^(1|true|yes|y|on)$ ]]; then
    echo "Resuming existing run directory: ${DUAL_DIR}"
    echo "Completed Snakemake outputs will be retained."
  else
    echo "ERROR: run directory already exists: ${DUAL_DIR}" >&2
    echo "Set DUAL_RESUME_EXISTING=true only after confirming no controller is active." >&2
    exit 2
  fi
fi

if [ ! -s "${SNAKEFILE}" ] || [ ! -s "${CLUSTER_CONFIG}" ]; then
  echo "ERROR: Snakefile or cluster configuration is missing." >&2
  exit 3
fi

export PROJECT_ROOT
export DUAL_RUN_ID
export DUAL_DIR
export DUAL_SAMPLE_TABLE
export DUAL_INPUT_MODE
export DUAL_RAW_DIR
export DUAL_TRIMMED_DIR
export DUAL_HOST_FASTA
export DUAL_HOST_GTF
export DUAL_FUNGUS_REFERENCE_DIR
export DUAL_FUNGUS_FASTA
export DUAL_FUNGUS_GTF
export DUAL_KRAKEN_DB
export DUAL_RUN_BRACKEN
export DUAL_KRAKEN_SOURCES
export DUAL_FEATURECOUNTS_STRAND
export DUAL_RRNA_LIST
export DUAL_RESUME_EXISTING
export DUAL_REQUIRE_R
export LOCUST_SPECIES="gregaria"
export LOCUST_TISSUES="fatbody"
export PROJECTDIR
export REFDIR
export SCRATCHDIR
export METATX_CONDA_MODULE="${METATX_CONDA_MODULE:-mamba}"
export METATX_CONDA_ENV="${METATX_CONDA_ENV:-/scratch/mtecher/conda-envs/metatranscriptomics-qc}"
export FASTP_MODULE="${FASTP_MODULE:-}"
export FASTP_CONDA_MODULE="${FASTP_CONDA_MODULE:-mamba}"
export FASTP_CONDA_ENV="${FASTP_CONDA_ENV:-${SNAKEMAKE_ENV}}"
export FASTP_CMD="${FASTP_CMD:-fastp}"

mkdir -p "${CLUSTER_LOG_DIR}"

if command -v module >/dev/null 2>&1; then
  module purge
  module load mamba
fi
source activate "${SNAKEMAKE_ENV}"

echo "Run ID: ${DUAL_RUN_ID}"
echo "Target: ${TARGET}"
echo "Output: ${DUAL_DIR}"
echo "Resume existing run: ${DUAL_RESUME_EXISTING}"
echo "Raw reads: ${DUAL_RAW_DIR}"
echo "Input mode: ${DUAL_INPUT_MODE}"
echo "External trimmed reads (trimmed mode only): ${DUAL_TRIMMED_DIR}"
echo "Host reference: ${DUAL_HOST_FASTA}"
echo "Fungal reference: ${DUAL_FUNGUS_FASTA}"
echo "Kraken database: ${DUAL_KRAKEN_DB}"
echo "Host rRNA exclusion list: ${DUAL_RRNA_LIST}"
echo "Mapping branches: fresh host-only and competitive host-pathogen"
echo "R analysis required on HPC: ${DUAL_REQUIRE_R}"

bash pilot/check_host_pathogen_dual_hpc.sh

snakemake -s "${SNAKEFILE}" "${TARGET}" \
  --dry-run \
  --rerun-incomplete \
  --rerun-triggers mtime

snakemake -s "${SNAKEFILE}" "${TARGET}" \
  --jobs "${JOBS}" \
  --latency-wait "${LATENCY_WAIT}" \
  --rerun-incomplete \
  --keep-going \
  --rerun-triggers mtime \
  --cluster-cancel scancel \
  --cluster-config "${CLUSTER_CONFIG}" \
  --cluster "sbatch --output=${CLUSTER_LOG_DIR}/%x-%j.out --error=${CLUSTER_LOG_DIR}/%x-%j.err -N {cluster.nodes} -n {cluster.ntasks} -c {cluster.cpus-per-task} --mem={cluster.mem} -t {cluster.time} -p {cluster.partition} --qos={cluster.qos} --export=ALL"

echo "Finished: $(date)"
