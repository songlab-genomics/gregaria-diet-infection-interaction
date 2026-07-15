#!/usr/bin/env bash
set -euo pipefail

# Preflight check for the fat body pilot run on the cluster. Run this from the
# repository code/ folder before submitting snakemake.fatbody_pilot.slurm.

RAW_DIR="${RAW_DIR:-/data/songlab/sequencing_data/RNAseq/mehreen}"
PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
RUN_ID="${RUN_ID:-fatbody_pilot_preflight}"
LINK_DATADIR="${LINK_DATADIR:-${PROJECT_ROOT}/raw_links/${RUN_ID}}"
WORKDIR="${WORKDIR:-${PROJECT_ROOT}/output/runs/${RUN_ID}}"
SCRATCHDIR="${SCRATCHDIR:-${PROJECT_ROOT}/scratch/${RUN_ID}}"
CLUSTER_LOG_DIR="${CLUSTER_LOG_DIR:-${PROJECT_ROOT}/logs/slurm/${RUN_ID}}"
SAMPLE_TABLE="${SAMPLE_TABLE:-pilot/fatbody_pilot_samples.tsv}"
SNAKEMAKE_ENV="${SNAKEMAKE_ENV:-myENV}"
TARGET="${TARGET:-expression_all}"

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
REF_GENOME="${REFDIR}/${GENOME_ID}_genomic.fna"
REF_GTF="${REFDIR}/${GENOME_ID}_genomic.gtf"
STAR_INDEX_DIR="${REFDIR}/index/gregaria/STAR"

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

echo "Fat body pilot preflight"
echo "Working directory: $(pwd)"
echo "Raw reads: ${RAW_DIR}"
echo "Project root: ${PROJECT_ROOT}"
echo "Reference directory: ${REFDIR}"
echo "Pilot sample table: ${SAMPLE_TABLE}"
echo "Preflight link DATADIR: ${LINK_DATADIR}"
echo "Preflight WORKDIR: ${WORKDIR}"
echo "Cluster log directory: ${CLUSTER_LOG_DIR}"
echo

[ -d "${RAW_DIR}" ] || fail "raw read directory not found: ${RAW_DIR}"
[ -d "${PROJECT_ROOT}" ] || fail "project root not found: ${PROJECT_ROOT}"
[ -s "${SAMPLE_TABLE}" ] || fail "pilot sample table not found: ${SAMPLE_TABLE}"
[ -s "Snakefile" ] || fail "Snakefile not found; run this from the code/ folder"
[ -s "cluster.json" ] || fail "cluster.json not found; run this from the code/ folder"
[ -s "${REF_GENOME}" ] || fail "reference genome not found: ${REF_GENOME}"
[ -s "${REF_GTF}" ] || fail "reference GTF not found: ${REF_GTF}"

echo "Checking pilot raw reads..."
sample_count=0
while IFS=$'\t' read -r sample_id raw_prefix tissue treatment diet note; do
  if [ "${sample_id}" = "sample_id" ]; then
    continue
  fi
  sample_count=$((sample_count + 1))
  for read_num in 1 2; do
    raw_file="${RAW_DIR}/${raw_prefix}_MERGE_${read_num}.fq.gz"
    [ -s "${raw_file}" ] || fail "missing or empty raw read file: ${raw_file}"
  done
done < "${SAMPLE_TABLE}"
echo "  Found ${sample_count} pilot samples with both read pairs."

echo "Preparing run-specific symlink layout..."
bash pilot/prepare_fatbody_pilot_reads.sh "${RAW_DIR}" "${LINK_DATADIR}" "${SAMPLE_TABLE}"

echo "Checking symlink layout expected by Snakemake..."
link_count=$(find "${LINK_DATADIR}/00-reads-fatbody" -type l -name '*_MERGE_*.fastq.gz' | wc -l | tr -d ' ')
expected_links=$((sample_count * 2))
[ "${link_count}" -eq "${expected_links}" ] || fail "expected ${expected_links} FASTQ symlinks, found ${link_count}"
echo "  Found ${link_count} FASTQ symlinks."

if [ -d "${STAR_INDEX_DIR}" ]; then
  echo "STAR index exists: ${STAR_INDEX_DIR}"
else
  echo "STAR index not found at ${STAR_INDEX_DIR}; Snakemake will build it before alignment."
fi

if command -v sbatch >/dev/null 2>&1; then
  echo "SLURM sbatch is available."
else
  echo "WARNING: sbatch not found in this shell; run on a login node before submitting."
fi

if command -v module >/dev/null 2>&1; then
  module purge
  module load mamba
else
  echo "WARNING: module command is not available; cannot load mamba automatically."
fi

echo "Activating Snakemake controller environment: ${SNAKEMAKE_ENV}"
source activate "${SNAKEMAKE_ENV}"

export LOCUST_SPECIES="gregaria"
export LOCUST_TISSUES="fatbody"
export SAMPLE_DISCOVERY="raw_reads"
export DATADIR="${LINK_DATADIR}"
export WORKDIR="${WORKDIR}"
export SCRATCHDIR="${SCRATCHDIR}"
export PROJECTDIR="${PROJECTDIR}"
export REFDIR="${REFDIR}"
export FASTP_CONDA_ENV="${FASTP_CONDA_ENV:-${SNAKEMAKE_ENV}}"

echo "Checking fastp in environment: ${FASTP_CONDA_ENV}"
if ! command -v fastp >/dev/null 2>&1; then
  fail "fastp was not found after activating ${FASTP_CONDA_ENV}; set FASTP_CONDA_ENV, FASTP_MODULE, or FASTP_CMD before launching"
fi
fastp --version || true

echo "Running Snakemake dry-run for target: ${TARGET}"
snakemake -s Snakefile "${TARGET}" \
  --dry-run \
  --printshellcmds \
  --reason \
  --cores 1

echo
echo "Preflight passed. You can launch the pilot with:"
echo "  sbatch snakemake.fatbody_pilot.slurm"
