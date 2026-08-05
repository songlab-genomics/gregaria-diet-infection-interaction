#!/usr/bin/env bash
#SBATCH --job-name=homology_db
#SBATCH --output=slurm-homology-db-%j.out
#SBATCH --error=slurm-homology-db-%j.err
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=6-23:59:00
#SBATCH --export=ALL

set -euo pipefail

ENV_DIR="${HOMOLOGY_CONDA_ENV:-/scratch/mtecher/conda-envs/scaffold-homology}"
DB_DIR="${HOMOLOGY_BLAST_DB_ROOT:-/scratch/mtecher/scaffold_homology_databases/ncbi_blastdb}"
ARCHIVE_DIR="${HOMOLOGY_CLUSTER_ARCHIVE_DIR:-${DB_DIR}/archives/nr_cluster_seq}"
METADATA_DIR="${DB_DIR}/metadata"
EXTRACT_MARKER_DIR="${DB_DIR}/.nr_cluster_seq_extracted"
NR_BASE_URL="https://ftp.ncbi.nlm.nih.gov/blast/db/experimental"
NR_METADATA_URL="${NR_BASE_URL}/nr_cluster_seq-prot-metadata.json"

module purge
module load mamba
source activate "${ENV_DIR}"

mkdir -p "${DB_DIR}" "${ARCHIVE_DIR}" "${METADATA_DIR}" "${EXTRACT_MARKER_DIR}"
cd "${DB_DIR}"

echo "Downloading or resuming NCBI core_nt..."
update_blastdb.pl --decompress --num_threads "${SLURM_CPUS_PER_TASK:-8}" core_nt

timestamp=$(date +%Y%m%d_%H%M%S)
metadata_file="${METADATA_DIR}/nr_cluster_seq-prot-metadata_${timestamp}.json"
volume_list="${METADATA_DIR}/nr_cluster_seq-volumes_${timestamp}.txt"

echo "Reading the current NCBI nr_cluster_seq volume manifest..."
wget --no-verbose -O "${metadata_file}" "${NR_METADATA_URL}"
python3 -c '
import json
import pathlib
import sys

metadata = json.load(open(sys.argv[1], encoding="utf-8"))
if metadata.get("dbname") != "nr_cluster_seq":
    raise SystemExit("Unexpected NCBI database name in metadata")
files = sorted(metadata.get("files", []))
if len(files) != int(metadata.get("number-of-volumes", -1)):
    raise SystemExit("NCBI volume count does not match the file list")
for url in files:
    print(pathlib.PurePosixPath(url).name)
' "${metadata_file}" > "${volume_list}"

volume_count=$(wc -l < "${volume_list}" | tr -d ' ')
echo "Downloading or resuming ${volume_count} NCBI nr_cluster_seq volumes..."
echo "Archive directory: ${ARCHIVE_DIR}"

while IFS= read -r archive_name; do
  [ -n "${archive_name}" ] || continue
  archive_path="${ARCHIVE_DIR}/${archive_name}"
  md5_path="${archive_path}.md5"

  if [ ! -s "${md5_path}" ]; then
    wget --no-verbose -O "${md5_path}.partial" \
      "${NR_BASE_URL}/${archive_name}.md5"
    mv "${md5_path}.partial" "${md5_path}"
  fi

  (
    cd "${ARCHIVE_DIR}"
    wget --continue --no-verbose "${NR_BASE_URL}/${archive_name}"
  )

  (
    cd "${ARCHIVE_DIR}"
    md5sum --check "${archive_name}.md5"
  )

  extract_marker="${EXTRACT_MARKER_DIR}/${archive_name}.ok"
  if [ ! -s "${extract_marker}" ]; then
    echo "Extracting ${archive_name}..."
    tar -xzf "${archive_path}" -C "${DB_DIR}"
    printf 'extracted_at\t%s\narchive\t%s\n' \
      "$(date --iso-8601=seconds)" "${archive_path}" > "${extract_marker}"
  fi
done < "${volume_list}"

echo "Validating both BLAST databases..."
blastdbcmd -db core_nt -info
blastdbcmd -db nr_cluster_seq -info

completion_tmp="${DB_DIR}/download_complete.tsv.partial"
printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)" > "${completion_tmp}"
printf 'core_nt_db\t%s\n' "${DB_DIR}/core_nt" >> "${completion_tmp}"
printf 'clustered_nr_db\t%s\n' "${DB_DIR}/nr_cluster_seq" >> "${completion_tmp}"
printf 'nr_cluster_seq_metadata\t%s\n' "${metadata_file}" >> "${completion_tmp}"
printf 'nr_cluster_seq_volumes\t%s\n' "${volume_count}" >> "${completion_tmp}"
mv "${completion_tmp}" "${DB_DIR}/download_complete.tsv"

echo "Homology databases are ready: ${DB_DIR}"
echo "ClusteredNR local database basename: ${DB_DIR}/nr_cluster_seq"
