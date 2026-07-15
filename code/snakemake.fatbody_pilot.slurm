#!/usr/bin/env bash
#SBATCH --job-name=fatbody_pilot_smk
#SBATCH --output=slurm-fatbody-pilot-%j.out
#SBATCH --error=slurm-fatbody-pilot-%j.err
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --export=NONE

set -euo pipefail

# Submit this script from the repository code/ folder on the cluster.
if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -s "${SLURM_SUBMIT_DIR}/cluster.json" ]; then
  cd "${SLURM_SUBMIT_DIR}"
else
  cd "$(dirname "${BASH_SOURCE[0]}")"
fi

TARGET="${1:-expression_all}"
JOBS="${2:-12}"
SNAKEMAKE_ENV="${3:-myENV}"
SNAKEFILE="${4:-Snakefile}"
LATENCY_WAIT="${LATENCY_WAIT:-120}"
CLUSTER_CONFIG="${CLUSTER_CONFIG:-$(pwd)/cluster.json}"

RUN_ID="${RUN_ID:-fatbody_pilot_$(date +%Y%m%d_%H%M%S)}"
RAW_DIR="${RAW_DIR:-/data/songlab/sequencing_data/RNAseq/mehreen}"
PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
LINK_DATADIR="${LINK_DATADIR:-${PROJECT_ROOT}/raw_links/${RUN_ID}}"
WORKDIR="${WORKDIR:-${PROJECT_ROOT}/output/runs/${RUN_ID}}"
SCRATCHDIR="${SCRATCHDIR:-${PROJECT_ROOT}/scratch/${RUN_ID}}"
CLUSTER_LOG_DIR="${CLUSTER_LOG_DIR:-${PROJECT_ROOT}/logs/slurm/${RUN_ID}}"
SLURM_EXPORT_VARS="${SLURM_EXPORT_VARS:-ALL,LOCUST_SPECIES,LOCUST_TISSUES,SAMPLE_DISCOVERY,DATADIR,WORKDIR,SCRATCHDIR,PROJECTDIR,REFDIR,FASTP_MODULE,FASTP_CONDA_MODULE,FASTP_CONDA_ENV,FASTP_CMD}"

GENOME_ID="GCF_023897955.1_iqSchGreg1.2"

resolve_refdir() {
  # Respect an explicit REFDIR first. The AGY/time-course pipeline keeps
  # references under /data/songlab/maeva/{species}-timecourse/reference.
  if [ -n "${REFDIR:-}" ]; then
    echo "${REFDIR}"
    return
  fi

  for candidate in \
    "/data/songlab/maeva/gregaria-timecourse/reference" \
    "/scratch/mtecher/gregaria-timecourse/reference" \
    "/scratch/mtecher/locust-time-course-RNAseq/reference" \
    "/scratch/mtecher/locust-time-course-RNAseq/data/reference" \
    "/scratch/mtecher/locust-time-course-RNAseq/references"; do
    if [ -s "${candidate}/${GENOME_ID}_genomic.fna" ] && [ -s "${candidate}/${GENOME_ID}_genomic.gtf" ]; then
      echo "${candidate}"
      return
    fi
  done

  # If no complete reference is found, return the most likely location so the
  # downstream error message prints a useful path.
  echo "/data/songlab/maeva/gregaria-timecourse/reference"
}

REFDIR="$(resolve_refdir)"
PROJECTDIR="${PROJECTDIR:-$(dirname "${REFDIR}")}"

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Run ID: ${RUN_ID}"
echo "Snakefile: ${SNAKEFILE}"
echo "Target: ${TARGET}"
echo "Concurrent jobs: ${JOBS}"
echo "Raw reads: ${RAW_DIR}"
echo "Symlink DATADIR: ${LINK_DATADIR}"
echo "Output WORKDIR: ${WORKDIR}"
echo "Reference directory: ${REFDIR}"
echo "Cluster job logs: ${CLUSTER_LOG_DIR}"
echo "Worker export variables: ${SLURM_EXPORT_VARS}"

if [ ! -s "${CLUSTER_CONFIG}" ]; then
  echo "ERROR: cluster config was not found: ${CLUSTER_CONFIG}" >&2
  exit 2
fi

if [ ! -s "${SNAKEFILE}" ]; then
  echo "ERROR: Snakefile was not found: ${SNAKEFILE}" >&2
  exit 2
fi

bash pilot/prepare_fatbody_pilot_reads.sh "${RAW_DIR}" "${LINK_DATADIR}" "pilot/fatbody_pilot_samples.tsv"
mkdir -p "${CLUSTER_LOG_DIR}"

export LOCUST_SPECIES="gregaria"
export LOCUST_TISSUES="fatbody"
export SAMPLE_DISCOVERY="raw_reads"
export DATADIR="${LINK_DATADIR}"
export WORKDIR="${WORKDIR}"
export SCRATCHDIR="${SCRATCHDIR}"
export PROJECTDIR="${PROJECTDIR}"
export REFDIR="${REFDIR}"
export FASTP_MODULE="${FASTP_MODULE:-}"
export FASTP_CONDA_MODULE="${FASTP_CONDA_MODULE:-mamba}"
export FASTP_CONDA_ENV="${FASTP_CONDA_ENV:-${SNAKEMAKE_ENV}}"
export FASTP_CMD="${FASTP_CMD:-fastp}"

if command -v module >/dev/null 2>&1; then
  module purge
  module load mamba
fi
source activate "${SNAKEMAKE_ENV}"

snakemake -s "${SNAKEFILE}" "${TARGET}" \
  --jobs "${JOBS}" \
  --latency-wait "${LATENCY_WAIT}" \
  --rerun-incomplete \
  --rerun-triggers mtime \
  --cluster-cancel scancel \
  --cluster-config "${CLUSTER_CONFIG}" \
  --cluster "sbatch --output=${CLUSTER_LOG_DIR}/%x-%j.out --error=${CLUSTER_LOG_DIR}/%x-%j.err -N {cluster.nodes} -n {cluster.ntasks} -c {cluster.cpus-per-task} --mem={cluster.mem} -t {cluster.time} -p {cluster.partition} --qos={cluster.qos} --export=${SLURM_EXPORT_VARS}"

echo "Finished: $(date)"
