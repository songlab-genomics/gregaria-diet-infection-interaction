#!/usr/bin/env bash
set -euo pipefail

# Download the annotated M. robertsii ARSEF 23 RefSeq package used as the
# mapping proxy for experimental isolate DWR2009 / ARSEF 10343.

ACCESSION="${ACCESSION:-GCF_000187425.2}"
PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
REFERENCE_ROOT="${REFERENCE_ROOT:-${PROJECT_ROOT}/reference/metarhizium/${ACCESSION}}"
PACKAGE_ZIP="${REFERENCE_ROOT}/${ACCESSION}_ncbi_dataset.zip"
UNPACK_DIR="${REFERENCE_ROOT}/ncbi_package"

if ! command -v datasets >/dev/null 2>&1; then
  echo "ERROR: NCBI datasets is not available." >&2
  echo "Install with: mamba install -c conda-forge ncbi-datasets-cli" >&2
  exit 2
fi

if [ -e "${REFERENCE_ROOT}/${ACCESSION}_genomic.fna" ] || \
   [ -e "${REFERENCE_ROOT}/${ACCESSION}_genomic.gtf" ]; then
  echo "ERROR: standardized reference files already exist in ${REFERENCE_ROOT}" >&2
  echo "Choose a different REFERENCE_ROOT rather than overwriting them." >&2
  exit 3
fi

mkdir -p "${REFERENCE_ROOT}"

datasets download genome accession "${ACCESSION}" \
  --include genome,gtf,protein,rna,seq-report \
  --filename "${PACKAGE_ZIP}" \
  --no-progressbar

unzip -q "${PACKAGE_ZIP}" -d "${UNPACK_DIR}"

ASSEMBLY_DIR="${UNPACK_DIR}/ncbi_dataset/data/${ACCESSION}"
GENOME_FILE="$(find "${ASSEMBLY_DIR}" -maxdepth 1 -type f -name '*_genomic.fna' -print -quit)"
GTF_FILE="$(find "${ASSEMBLY_DIR}" -maxdepth 1 -type f -name 'genomic.gtf' -print -quit)"
PROTEIN_FILE="$(find "${ASSEMBLY_DIR}" -maxdepth 1 -type f -name 'protein.faa' -print -quit)"

if [ -z "${GENOME_FILE}" ] || [ -z "${GTF_FILE}" ]; then
  echo "ERROR: the NCBI package did not contain both genomic FASTA and GTF." >&2
  exit 4
fi

ln -s "${GENOME_FILE}" "${REFERENCE_ROOT}/${ACCESSION}_genomic.fna"
ln -s "${GTF_FILE}" "${REFERENCE_ROOT}/${ACCESSION}_genomic.gtf"
if [ -n "${PROTEIN_FILE}" ]; then
  ln -s "${PROTEIN_FILE}" "${REFERENCE_ROOT}/${ACCESSION}_protein.faa"
fi

echo "Metarhizium reference ready:"
echo "  FASTA: ${REFERENCE_ROOT}/${ACCESSION}_genomic.fna"
echo "  GTF:   ${REFERENCE_ROOT}/${ACCESSION}_genomic.gtf"
echo "Experimental isolate: DWR2009 / ARSEF 10343"
echo "Mapping proxy: ARSEF 23 (${ACCESSION})"
