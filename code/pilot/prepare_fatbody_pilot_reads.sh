#!/usr/bin/env bash
set -euo pipefail

# Create run-specific symlinks from the Mehreen raw-read folder into the naming
# layout used by the copied time-course Snakemake rules.

RAW_DIR="${1:-/data/songlab/sequencing_data/RNAseq/mehreen}"
LINK_DATADIR="${2:-/scratch/mtecher/gregaria-diet-infection-interaction/raw_links/fatbody_pilot}"
SAMPLE_TABLE="${3:-pilot/fatbody_pilot_samples.tsv}"

if [ ! -s "${SAMPLE_TABLE}" ]; then
  echo "ERROR: pilot sample table not found: ${SAMPLE_TABLE}" >&2
  exit 2
fi

while IFS=$'\t' read -r sample_id raw_prefix tissue treatment diet note; do
  if [ "${sample_id}" = "sample_id" ]; then
    continue
  fi

  dest_dir="${LINK_DATADIR}/00-reads-${tissue}"
  mkdir -p "${dest_dir}"

  for read_num in 1 2; do
    src="${RAW_DIR}/${raw_prefix}_MERGE_${read_num}.fq.gz"
    dst="${dest_dir}/${sample_id}_MERGE_${read_num}.fastq.gz"

    if [ ! -s "${src}" ]; then
      echo "ERROR: missing raw read file: ${src}" >&2
      exit 3
    fi

    ln -sfn "${src}" "${dst}"
  done
done < "${SAMPLE_TABLE}"

echo "Pilot raw-read links are ready under: ${LINK_DATADIR}"
