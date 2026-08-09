#!/usr/bin/env bash
#SBATCH --job-name=merge_1044
#SBATCH --output=slurm-merge-1044-%j.out
#SBATCH --error=slurm-merge-1044-%j.err
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00

set -euo pipefail

# Set SOURCE_DIR to the uploaded folder that contains the four L1/L2 FASTQs.
SOURCE_DIR="${SOURCE_DIR:?Export SOURCE_DIR for the uploaded mehreen_1044 folder}"
DEST_DIR="${DEST_DIR:-/data/songlab/sequencing_data/RNAseq/mehreen}"
SAMPLE_PREFIX="${SAMPLE_PREFIX:-mehreen_1044}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mapfile -t R1_LANES < <(find "${SOURCE_DIR}" -maxdepth 1 -type f \
  -name "${SAMPLE_PREFIX}_*_L*_1.fq.gz" -print | sort)
mapfile -t R2_LANES < <(find "${SOURCE_DIR}" -maxdepth 1 -type f \
  -name "${SAMPLE_PREFIX}_*_L*_2.fq.gz" -print | sort)

if [ "${#R1_LANES[@]}" -ne 2 ] || [ "${#R2_LANES[@]}" -ne 2 ]; then
  echo "ERROR: expected two R1 lanes and two R2 lanes in ${SOURCE_DIR}." >&2
  printf 'R1: %s\n' "${R1_LANES[@]:-none}" >&2
  printf 'R2: %s\n' "${R2_LANES[@]:-none}" >&2
  exit 2
fi

for input_fastq in "${R1_LANES[@]}" "${R2_LANES[@]}"; do
  pigz -t "${input_fastq}"
done

mkdir -p "${DEST_DIR}"
OUT_R1="${DEST_DIR}/${SAMPLE_PREFIX}_MERGE_1.fq.gz"
OUT_R2="${DEST_DIR}/${SAMPLE_PREFIX}_MERGE_2.fq.gz"

# Existing merged reads are protected. Rename them yourself only after an
# independent audit; this job never overwrites inputs.
if [ -e "${OUT_R1}" ] || [ -e "${OUT_R2}" ]; then
  echo "ERROR: a destination file already exists; nothing was overwritten." >&2
  echo "  ${OUT_R1}" >&2
  echo "  ${OUT_R2}" >&2
  exit 3
fi

TMP_R1="${OUT_R1}.partial.${SLURM_JOB_ID:-manual}"
TMP_R2="${OUT_R2}.partial.${SLURM_JOB_ID:-manual}"
trap 'rm -f "${TMP_R1}" "${TMP_R2}"' EXIT

# Decompress both lanes in lane order and recompress one FASTQ per mate.
pigz -dc "${R1_LANES[@]}" | pigz -p "${THREADS}" > "${TMP_R1}"
pigz -dc "${R2_LANES[@]}" | pigz -p "${THREADS}" > "${TMP_R2}"
pigz -t "${TMP_R1}"
pigz -t "${TMP_R2}"

mv "${TMP_R1}" "${OUT_R1}"
mv "${TMP_R2}" "${OUT_R2}"
trap - EXIT

echo "Created protected merged inputs:"
ls -lh "${OUT_R1}" "${OUT_R2}"
echo "Detailed lane-content reconciliation is deferred and can be run later."
