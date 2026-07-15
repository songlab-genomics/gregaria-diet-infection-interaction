#!/usr/bin/env bash
#SBATCH --job-name=eggnog_smk
#SBATCH --output=slurm-eggnog-%j.out
#SBATCH --error=slurm-eggnog-%j.err
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --export=NONE

set -euo pipefail

# SLURM copies the submitted script into /var/spool before running it, so $0 is
# not a reliable way to find the repository. Use the submit directory first.
if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -s "${SLURM_SUBMIT_DIR}/cluster.json" ]; then
  cd "${SLURM_SUBMIT_DIR}"
else
  cd "$(dirname "${BASH_SOURCE[0]}")"
fi

TARGET="${1:-eggnog_desertlocustr_all}"
JOBS="${2:-8}"
SNAKEMAKE_ENV="${3:-myENV}"
SNAKEFILE="${4:-Snakefile.eggnog_annotations}"
EGGNOG_RUN_ID_ARG="${5:-}"
EGGNOG_DATA_DIR_ARG="${6:-}"
LATENCY_WAIT="${LATENCY_WAIT:-120}"
CLUSTER_CONFIG="${CLUSTER_CONFIG:-$(pwd)/cluster.json}"

if [ -n "${EGGNOG_RUN_ID_ARG}" ]; then
  export EGGNOG_RUN_ID="${EGGNOG_RUN_ID_ARG}"
fi
if [ -n "${EGGNOG_DATA_DIR_ARG}" ]; then
  export EGGNOG_DATA_DIR="${EGGNOG_DATA_DIR_ARG}"
fi

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Snakefile: ${SNAKEFILE}"
echo "Target: ${TARGET}"
echo "Concurrent jobs: ${JOBS}"
echo "Snakemake controller environment: ${SNAKEMAKE_ENV}"
echo "EGGNOG_RUN_ID: ${EGGNOG_RUN_ID:-config default}"
echo "EGGNOG_DATA_DIR: ${EGGNOG_DATA_DIR:-config default}"
echo "Cluster config: ${CLUSTER_CONFIG}"

if [ ! -s "${CLUSTER_CONFIG}" ]; then
  echo "ERROR: cluster config was not found: ${CLUSTER_CONFIG}" >&2
  echo "Copy code/cluster.json to this HPC code directory before launching." >&2
  exit 2
fi

if [ ! -s "${SNAKEFILE}" ]; then
  echo "ERROR: Snakefile was not found: ${SNAKEFILE}" >&2
  exit 2
fi

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
  --cluster "sbatch -N {cluster.nodes} -n {cluster.ntasks} -c {cluster.cpus-per-task} --mem={cluster.mem} -t {cluster.time} -p {cluster.partition} --qos={cluster.qos} --export={cluster.export}"

echo "Finished: $(date)"
