#!/usr/bin/env bash
set -euo pipefail

# Read-only status report for the shared core_nt and ClusteredNR databases.

DB_DIR="${HOMOLOGY_BLAST_DB_ROOT:-/scratch/mtecher/scaffold_homology_databases/ncbi_blastdb}"
ARCHIVE_DIR="${HOMOLOGY_CLUSTER_ARCHIVE_DIR:-${DB_DIR}/archives/nr_cluster_seq}"
EXTRACT_MARKER_DIR="${DB_DIR}/.nr_cluster_seq_extracted"
ENV_DIR="${HOMOLOGY_CONDA_ENV:-/scratch/mtecher/conda-envs/scaffold-homology}"

count_files() {
  local directory="$1"
  local pattern="$2"
  find "${directory}" -maxdepth 1 -type f -name "${pattern}" -size +0c \
    2>/dev/null | wc -l | tr -d ' '
}

latest_volume_list=$(find "${DB_DIR}/metadata" -maxdepth 1 -type f \
  -name 'nr_cluster_seq-volumes_*.txt' -size +0c 2>/dev/null \
  | sort | tail -n 1)

expected="unknown"
if [ -n "${latest_volume_list}" ]; then
  expected=$(awk 'NF {n++} END {print n + 0}' "${latest_volume_list}")
fi

archives=$(count_files "${ARCHIVE_DIR}" 'nr_cluster_seq.*.tar.gz')
checksums=$(count_files "${ARCHIVE_DIR}" 'nr_cluster_seq.*.tar.gz.md5')
extracted=$(count_files "${EXTRACT_MARKER_DIR}" 'nr_cluster_seq.*.tar.gz.ok')

echo "Scaffold homology database status"
echo "Database root: ${DB_DIR}"
echo "ClusteredNR expected volumes: ${expected}"
echo "ClusteredNR archives present: ${archives}"
echo "ClusteredNR checksums present: ${checksums}"
echo "ClusteredNR volumes extracted: ${extracted}"
echo "Completion marker: ${DB_DIR}/download_complete.tsv"
if [ -s "${DB_DIR}/download_complete.tsv" ]; then
  sed 's/^/  /' "${DB_DIR}/download_complete.tsv"
else
  echo "  NOT COMPLETE"
fi

echo
df -h "${DB_DIR}"

echo
if command -v module >/dev/null 2>&1; then
  module purge
  module load mamba
fi
source activate "${ENV_DIR}"

for database in core_nt nr_cluster_seq; do
  echo
  echo "blastdbcmd validation: ${database}"
  if blastdbcmd -db "${DB_DIR}/${database}" -info >/tmp/scaffold_homology_dbinfo.$$ 2>&1; then
    sed -n '1,4p' /tmp/scaffold_homology_dbinfo.$$
  else
    echo "  NOT READY"
    sed -n '1,3p' /tmp/scaffold_homology_dbinfo.$$
  fi
done
rm -f /tmp/scaffold_homology_dbinfo.$$

echo
echo "This check is read-only; it does not submit, modify, or remove database files."
