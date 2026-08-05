#!/usr/bin/env bash
set -euo pipefail

# Validate only the resources needed by the independent scaffold-origin audit.
# Raw reads, STAR indexes, BAMs, featureCounts, fungal references, and R are not
# inputs to this target.

PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
REFDIR="${REFDIR:-/data/songlab/maeva/gregaria-timecourse/reference}"
ORIGIN_CANDIDATE_LIST="${ORIGIN_CANDIDATE_LIST:-config/scaffold_origin/unified_unplaced_deg_scaffolds_20260729.txt}"
ORIGIN_CANDIDATE_GENES="${ORIGIN_CANDIDATE_GENES:-config/scaffold_origin/unified_unplaced_deg_genes_20260729.tsv}"
ORIGIN_CANDIDATE_PROVENANCE="${ORIGIN_CANDIDATE_PROVENANCE:-config/scaffold_origin/unified_unplaced_deg_scaffold_provenance_20260729.tsv}"
ORIGIN_NCBI_FCS_REPORT="${ORIGIN_NCBI_FCS_REPORT:-config/scaffold_origin/GCF_023897955.1_iqSchGreg1.2_fcs_report_20230323.txt}"
ORIGIN_HOST_FASTA="${ORIGIN_HOST_FASTA:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.fna}"
ORIGIN_HOST_ANNOTATION="${ORIGIN_HOST_ANNOTATION:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.gtf}"
ORIGIN_EXPECTED_SCAFFOLDS="${ORIGIN_EXPECTED_SCAFFOLDS:-77}"
ORIGIN_KRAKEN_DB="${ORIGIN_KRAKEN_DB:-/scratch/mtecher/kraken2/pluspf_20260626}"
ORIGIN_CONDA_MODULE="${ORIGIN_CONDA_MODULE:-mamba}"
ORIGIN_CONDA_ENV="${ORIGIN_CONDA_ENV:-/scratch/mtecher/conda-envs/metatranscriptomics-qc}"
ORIGIN_KRAKEN_CMD="${ORIGIN_KRAKEN_CMD:-kraken2}"
ORIGIN_MINIMAP2_CMD="${ORIGIN_MINIMAP2_CMD:-minimap2}"
ORIGIN_RUN_CROSS_SPECIES="${ORIGIN_RUN_CROSS_SPECIES:-true}"
ORIGIN_CROSS_SPECIES_REFERENCES="${ORIGIN_CROSS_SPECIES_REFERENCES:-piceifrons=/data/songlab/maeva/piceifrons-timecourse/reference/GCF_021461385.2_iqSchPice1.1_genomic.fna,americana=/data/songlab/maeva/americana-timecourse/reference/GCF_021461395.2_iqSchAmer2.1_genomic.fna,serialis_cubense=/data/songlab/maeva/cubense-timecourse/reference/GCF_023864345.2_iqSchSeri2.2_genomic.fna}"

for path in \
  "${ORIGIN_CANDIDATE_LIST}" \
  "${ORIGIN_CANDIDATE_GENES}" \
  "${ORIGIN_CANDIDATE_PROVENANCE}" \
  "${ORIGIN_NCBI_FCS_REPORT}" \
  "${ORIGIN_HOST_FASTA}" \
  "${ORIGIN_HOST_ANNOTATION}"; do
  if [ ! -s "${path}" ]; then
    echo "ERROR: missing or empty scaffold-audit input: ${path}" >&2
    exit 2
  fi
done

candidate_count=$(
  awk 'NF && $1 !~ /^#/ {print $1}' "${ORIGIN_CANDIDATE_LIST}" \
    | sort -u \
    | wc -l \
    | awk '{print $1}'
)
listed_count=$(
  awk 'NF && $1 !~ /^#/ {n++} END {print n + 0}' "${ORIGIN_CANDIDATE_LIST}"
)
if [ "${candidate_count}" -ne "${listed_count}" ]; then
  echo "ERROR: candidate scaffold list contains duplicates." >&2
  exit 3
fi
if [ "${candidate_count}" -ne "${ORIGIN_EXPECTED_SCAFFOLDS}" ]; then
  echo "ERROR: expected ${ORIGIN_EXPECTED_SCAFFOLDS} candidate scaffolds, found ${candidate_count}." >&2
  exit 4
fi
if awk 'NF && $1 !~ /^#/ && $1 !~ /^NW_/ {exit 1}' "${ORIGIN_CANDIDATE_LIST}"; then
  :
else
  echo "ERROR: candidate list contains an accession that is not an NW_ scaffold." >&2
  exit 5
fi

gene_count=$(awk -F '\t' 'NR > 1 && $1 != "" {n++} END {print n + 0}' "${ORIGIN_CANDIDATE_GENES}")
provenance_count=$(awk -F '\t' 'NR > 1 && $1 != "" {n++} END {print n + 0}' "${ORIGIN_CANDIDATE_PROVENANCE}")
if [ "${provenance_count}" -ne "${candidate_count}" ]; then
  echo "ERROR: provenance has ${provenance_count} rows for ${candidate_count} candidates." >&2
  exit 6
fi

for database_file in hash.k2d opts.k2d taxo.k2d nodes.dmp names.dmp; do
  if [ ! -s "${ORIGIN_KRAKEN_DB}/${database_file}" ]; then
    echo "ERROR: Kraken2 database file missing: ${ORIGIN_KRAKEN_DB}/${database_file}" >&2
    exit 7
  fi
done

cross_species_count=0
case "${ORIGIN_RUN_CROSS_SPECIES}" in
1|true|TRUE|yes|YES|y|Y|on|ON)
  IFS=',' read -r -a cross_species_items <<< "${ORIGIN_CROSS_SPECIES_REFERENCES}"
  for item in "${cross_species_items[@]}"; do
    if [[ "${item}" != *=* ]]; then
      echo "ERROR: cross-species reference must use species=/path/to/genome.fna: ${item}" >&2
      exit 11
    fi
    species="${item%%=*}"
    reference="${item#*=}"
    if [ -z "${species}" ] || [ ! -s "${reference}" ]; then
      echo "ERROR: missing cross-species reference for ${species}: ${reference}" >&2
      exit 11
    fi
    cross_species_count=$((cross_species_count + 1))
  done
  ;;
esac

for script in \
  scripts/stage_scaffold_origin_inputs.py \
  scripts/extract_unplaced_scaffolds.py \
  scripts/summarize_candidate_scaffold_annotation.py \
  scripts/summarize_scaffold_kraken.py \
  scripts/summarize_cross_species_scaffold_homology.py \
  scripts/build_scaffold_origin_audit.py \
  scripts/make_scaffold_origin_transfer_bundle.py; do
  if [ ! -s "${script}" ]; then
    echo "ERROR: workflow helper is missing: ${script}" >&2
    exit 8
  fi
done

(
  if command -v module >/dev/null 2>&1; then
    module purge
    if [ -n "${ORIGIN_CONDA_MODULE}" ]; then
      module load "${ORIGIN_CONDA_MODULE}"
    fi
  fi
  if [ -n "${ORIGIN_CONDA_ENV}" ]; then
    source activate "${ORIGIN_CONDA_ENV}"
  fi
  if ! command -v "${ORIGIN_KRAKEN_CMD}" >/dev/null 2>&1; then
    echo "ERROR: Kraken2 command unavailable: ${ORIGIN_KRAKEN_CMD}" >&2
    exit 9
  fi
  echo "  Kraken2 executable: $(command -v "${ORIGIN_KRAKEN_CMD}")"
  "${ORIGIN_KRAKEN_CMD}" --version | sed -n '1p'
  case "${ORIGIN_RUN_CROSS_SPECIES}" in
  1|true|TRUE|yes|YES|y|Y|on|ON)
    if ! command -v "${ORIGIN_MINIMAP2_CMD}" >/dev/null 2>&1; then
      echo "ERROR: minimap2 is unavailable in ${ORIGIN_CONDA_ENV}." >&2
      echo "Update the environment with: sbatch scripts/install_metatranscriptomics_qc_env.slurm" >&2
      exit 12
    fi
    echo "  minimap2 executable: $(command -v "${ORIGIN_MINIMAP2_CMD}")"
    "${ORIGIN_MINIMAP2_CMD}" --version | sed -n '1p'
    ;;
  esac
)

if ! command -v snakemake >/dev/null 2>&1; then
  echo "ERROR: snakemake is unavailable in the controller environment." >&2
  exit 10
fi

echo "Scaffold-origin audit preflight passed."
echo "  Candidate scaffolds: ${candidate_count}"
echo "  Candidate DEG genes: ${gene_count}"
echo "  Host FASTA: ${ORIGIN_HOST_FASTA}"
echo "  Host annotation: ${ORIGIN_HOST_ANNOTATION}"
echo "  Kraken database: ${ORIGIN_KRAKEN_DB}"
echo "  Cross-Schistocerca comparisons: ${cross_species_count}"
echo "  NCBI FCS report: ${ORIGIN_NCBI_FCS_REPORT}"
echo "  RNA-seq preprocessing dependencies: none"
echo "  R dependencies: none"
