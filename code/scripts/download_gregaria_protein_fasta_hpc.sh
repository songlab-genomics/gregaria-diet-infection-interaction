#!/usr/bin/env bash

set -euo pipefail

REF_DIR="${REFDIR:-/data/songlab/maeva/gregaria-timecourse/reference}"
ASSEMBLY="GCF_023897955.1_iqSchGreg1.2"
FTP_ROOT="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/897/955/GCF_023897955.1_iqSchGreg1.2"
ARCHIVE="${REF_DIR}/${ASSEMBLY}_protein.faa.gz"
OUTPUT="${REF_DIR}/${ASSEMBLY}_protein.faa"

mkdir -p "${REF_DIR}"

if [ -s "${OUTPUT}" ] && grep -q '^>' "${OUTPUT}"; then
  echo "Protein FASTA already present; leaving it unchanged: ${OUTPUT}"
  grep -c '^>' "${OUTPUT}"
  exit 0
fi

# -c resumes an interrupted download without replacing a completed file.
wget -c "${FTP_ROOT}/${ASSEMBLY}_protein.faa.gz" -O "${ARCHIVE}"
gzip -t "${ARCHIVE}"
pigz -dc "${ARCHIVE}" > "${OUTPUT}.tmp"
grep -q '^>' "${OUTPUT}.tmp"
mv "${OUTPUT}.tmp" "${OUTPUT}"

echo "Protein FASTA ready: ${OUTPUT}"
grep -c '^>' "${OUTPUT}"
