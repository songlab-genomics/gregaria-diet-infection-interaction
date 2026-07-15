#!/usr/bin/env bash
#SBATCH --job-name=metatx_qc_smk
#SBATCH --output=slurm-metatx-qc-%j.out
#SBATCH --error=slurm-metatx-qc-%j.err
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --export=NONE

set -euo pipefail

# SLURM runs the copied script from /var/spool, so use the submit directory to
# find Snakefiles and cluster.json.
if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -s "${SLURM_SUBMIT_DIR}/cluster.json" ]; then
  cd "${SLURM_SUBMIT_DIR}"
else
  cd "$(dirname "${BASH_SOURCE[0]}")"
fi

TARGET="${1:-metatranscriptomics_qc_all}"
JOBS="${2:-40}"
SNAKEMAKE_ENV="${3:-myENV}"
SNAKEFILE="${4:-Snakefile.metatranscriptomics_qc}"
LATENCY_WAIT="${LATENCY_WAIT:-120}"
CLUSTER_CONFIG="${CLUSTER_CONFIG:-$(pwd)/cluster.json}"

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Snakefile: ${SNAKEFILE}"
echo "Target: ${TARGET}"
echo "Concurrent jobs: ${JOBS}"
echo "Snakemake controller environment: ${SNAKEMAKE_ENV}"
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
