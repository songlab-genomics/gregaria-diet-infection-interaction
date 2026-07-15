#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_eggnog_desertlocustr.sh local [target] [cores] [run_id] [eggnog_data_dir]
#   ./run_eggnog_desertlocustr.sh slurm [target] [jobs] [controller_env] [run_id] [eggnog_data_dir]

cd "$(dirname "$0")"

MODE="${1:-local}"
TARGET="${2:-eggnog_desertlocustr_all}"
WORKERS="${3:-8}"
SNAKEFILE="${SNAKEFILE:-Snakefile.eggnog_annotations}"

case "${MODE}" in
  local)
    RUN_ID="${4:-${EGGNOG_RUN_ID:-}}"
    DATA_DIR="${5:-${EGGNOG_DATA_DIR:-}}"
    if [ -n "${RUN_ID}" ]; then export EGGNOG_RUN_ID="${RUN_ID}"; fi
    if [ -n "${DATA_DIR}" ]; then export EGGNOG_DATA_DIR="${DATA_DIR}"; fi
    echo "Running locally with ${WORKERS} cores"
    echo "Target: ${TARGET}"
    echo "EGGNOG_RUN_ID: ${EGGNOG_RUN_ID:-config default}"
    echo "EGGNOG_DATA_DIR: ${EGGNOG_DATA_DIR:-config default}"
    snakemake -s "${SNAKEFILE}" "${TARGET}" \
      --cores "${WORKERS}" \
      --latency-wait 120 \
      --rerun-incomplete \
      --rerun-triggers mtime \
      --printshellcmds \
      --reason
    ;;
  slurm)
    CONTROLLER_ENV="${4:-myENV}"
    RUN_ID="${5:-${EGGNOG_RUN_ID:-}}"
    DATA_DIR="${6:-${EGGNOG_DATA_DIR:-}}"
    echo "Submitting SLURM controller for target ${TARGET}"
    sbatch snakemake.eggnog.slurm "${TARGET}" "${WORKERS}" "${CONTROLLER_ENV}" "${SNAKEFILE}" "${RUN_ID}" "${DATA_DIR}"
    ;;
  dryrun)
    RUN_ID="${3:-${EGGNOG_RUN_ID:-}}"
    DATA_DIR="${4:-${EGGNOG_DATA_DIR:-}}"
    if [ -n "${RUN_ID}" ]; then export EGGNOG_RUN_ID="${RUN_ID}"; fi
    if [ -n "${DATA_DIR}" ]; then export EGGNOG_DATA_DIR="${DATA_DIR}"; fi
    echo "Dry-running locally"
    snakemake -s "${SNAKEFILE}" -n "${TARGET}" \
      --cores 1 \
      --rerun-incomplete \
      --rerun-triggers mtime \
      --printshellcmds \
      --reason
    ;;
  *)
    echo "Unknown mode: ${MODE}" >&2
    echo "Use local, dryrun, or slurm." >&2
    exit 2
    ;;
esac
