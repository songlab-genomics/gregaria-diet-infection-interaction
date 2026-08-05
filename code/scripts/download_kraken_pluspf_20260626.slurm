#!/usr/bin/env bash
#SBATCH --job-name=kraken_pluspf
#SBATCH --output=slurm-kraken-pluspf-%j.out
#SBATCH --error=slurm-kraken-pluspf-%j.err
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL

set -euo pipefail

# The 24-hour general job keeps the slow transfer off the login node. Marker
# files preserve completed stages if the job ever needs to be resubmitted.
KDB="${DUAL_KRAKEN_DB:-/scratch/mtecher/kraken2/pluspf_20260626}"
ARCHIVE="${KDB}/k2_pluspf_20260626.tar.gz"
ARCHIVE_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260626.tar.gz"
METADATA_URL="https://genome-idx.s3.amazonaws.com/kraken/pluspf_20260626"

mkdir -p "${KDB}"
cd "${KDB}"

for command_name in wget tar pigz md5sum flock; do
  if ! command -v "${command_name}" >/dev/null 2>&1; then
    echo "ERROR: required command is unavailable: ${command_name}" >&2
    exit 2
  fi
done

# Avoid two download jobs modifying the same database simultaneously.
exec 9>"${KDB}/.download.lock"
if ! flock -n 9; then
  echo "ERROR: another PlusPF download/extraction job is active." >&2
  exit 3
fi

echo "Kraken PlusPF release: 2026-06-26"
echo "Database directory: ${KDB}"
echo "Started: $(date)"

if [ ! -e "${KDB}/.archive_download_complete" ]; then
  echo "Downloading or resuming the 84.8 GB PlusPF archive..."
  wget \
    --continue \
    --tries=0 \
    --timeout=60 \
    --output-document="${ARCHIVE}" \
    "${ARCHIVE_URL}"
  test -s "${ARCHIVE}"
  touch "${KDB}/.archive_download_complete"
else
  echo "Archive download stage already complete."
fi

for metadata_file in inspect.txt library_report.tsv pluspf.md5; do
  if [ ! -s "${KDB}/${metadata_file}" ]; then
    wget \
      --continue \
      --tries=0 \
      --timeout=60 \
      --output-document="${KDB}/${metadata_file}" \
      "${METADATA_URL}/${metadata_file}"
  fi
done

if [ ! -e "${KDB}/.archive_extraction_complete" ]; then
  echo "Extracting the approximately 110 GB Kraken/Bracken database..."
  echo "Using pigz with ${SLURM_CPUS_PER_TASK:-8} threads."
  tar \
    --use-compress-program="pigz -d -p ${SLURM_CPUS_PER_TASK:-8}" \
    -xf "${ARCHIVE}" \
    -C "${KDB}"
  touch "${KDB}/.archive_extraction_complete"
else
  echo "Archive extraction stage already complete."
fi

if [ ! -e "${KDB}/.checksum_complete" ]; then
  echo "Validating extracted database files..."
  (
    cd "${KDB}"
    md5sum -c pluspf.md5
  )
  touch "${KDB}/.checksum_complete"
else
  echo "Checksum stage already complete."
fi

for database_file in \
  hash.k2d \
  opts.k2d \
  taxo.k2d \
  database150mers.kmer_distrib; do
  if [ ! -s "${KDB}/${database_file}" ]; then
    echo "ERROR: expected database file is missing: ${KDB}/${database_file}" >&2
    exit 4
  fi
done

echo "Focal fungal entries in the published inspection table:"
grep -Ei -m 10 'Metarhizium|Clavicipitaceae' "${KDB}/inspect.txt" || true

echo "Kraken PlusPF database is ready."
echo "Finished: $(date)"
