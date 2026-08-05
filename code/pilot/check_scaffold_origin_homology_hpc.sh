#!/usr/bin/env bash
set -euo pipefail

# Validate only the focused BLASTn/DIAMOND audit. No RNA-seq files are needed.

PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
REFDIR="${REFDIR:-/data/songlab/maeva/gregaria-timecourse/reference}"
HOMOLOGY_CANDIDATE_SCAFFOLDS="${HOMOLOGY_CANDIDATE_SCAFFOLDS:-config/scaffold_homology/infected_diet33_vs50_unresolved_scaffolds_20260731.txt}"
HOMOLOGY_CANDIDATE_GENES="${HOMOLOGY_CANDIDATE_GENES:-config/scaffold_homology/infected_diet33_vs50_unresolved_genes_20260731.tsv}"
HOMOLOGY_CANDIDATE_PROVENANCE="${HOMOLOGY_CANDIDATE_PROVENANCE:-config/scaffold_homology/infected_diet33_vs50_unresolved_provenance_20260731.tsv}"
HOMOLOGY_HOST_FASTA="${HOMOLOGY_HOST_FASTA:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.fna}"
HOMOLOGY_HOST_ANNOTATION="${HOMOLOGY_HOST_ANNOTATION:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.gtf}"
HOMOLOGY_HOST_PROTEIN_FASTA="${HOMOLOGY_HOST_PROTEIN_FASTA:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_protein.faa}"
HOMOLOGY_BLAST_DB_ROOT="${HOMOLOGY_BLAST_DB_ROOT:-/scratch/mtecher/scaffold_homology_databases/ncbi_blastdb}"
HOMOLOGY_CORE_NT_DB="${HOMOLOGY_CORE_NT_DB:-${HOMOLOGY_BLAST_DB_ROOT}/core_nt}"
HOMOLOGY_NR_CLUSTER_SEQ_DB="${HOMOLOGY_NR_CLUSTER_SEQ_DB:-${HOMOLOGY_BLAST_DB_ROOT}/nr_cluster_seq}"
HOMOLOGY_DATABASE_READY="${HOMOLOGY_DATABASE_READY:-${HOMOLOGY_BLAST_DB_ROOT}/download_complete.tsv}"
HOMOLOGY_TAXONOMY_NODES="${HOMOLOGY_TAXONOMY_NODES:-/scratch/mtecher/kraken2/pluspf_20260626/nodes.dmp}"
HOMOLOGY_TAXONOMY_NAMES="${HOMOLOGY_TAXONOMY_NAMES:-/scratch/mtecher/kraken2/pluspf_20260626/names.dmp}"
HOMOLOGY_CONDA_MODULE="${HOMOLOGY_CONDA_MODULE:-mamba}"
HOMOLOGY_CONDA_ENV="${HOMOLOGY_CONDA_ENV:-/scratch/mtecher/conda-envs/scaffold-homology}"
HOMOLOGY_EXPECTED_SCAFFOLDS="${HOMOLOGY_EXPECTED_SCAFFOLDS:-59}"
HOMOLOGY_EXPECTED_GENES="${HOMOLOGY_EXPECTED_GENES:-948}"

for path in \
  "${HOMOLOGY_CANDIDATE_SCAFFOLDS}" \
  "${HOMOLOGY_CANDIDATE_GENES}" \
  "${HOMOLOGY_CANDIDATE_PROVENANCE}" \
  "${HOMOLOGY_HOST_FASTA}" \
  "${HOMOLOGY_HOST_ANNOTATION}" \
  "${HOMOLOGY_HOST_PROTEIN_FASTA}" \
  "${HOMOLOGY_TAXONOMY_NODES}" \
  "${HOMOLOGY_TAXONOMY_NAMES}"; do
  if [ ! -s "${path}" ]; then
    echo "ERROR: missing or empty homology-audit input: ${path}" >&2
    exit 2
  fi
done

if [ ! -s "${HOMOLOGY_DATABASE_READY}" ]; then
  echo "ERROR: homology database setup is incomplete: ${HOMOLOGY_DATABASE_READY}" >&2
  echo "The local ClusteredNR database basename must be nr_cluster_seq." >&2
  echo "Submit: sbatch scripts/download_scaffold_homology_databases.slurm" >&2
  exit 2
fi

scaffold_count=$(awk 'NF && $1 !~ /^#/ {n++} END {print n + 0}' "${HOMOLOGY_CANDIDATE_SCAFFOLDS}")
gene_count=$(awk -F '\t' 'NR > 1 && $1 != "" {n++} END {print n + 0}' "${HOMOLOGY_CANDIDATE_GENES}")
if [ "${scaffold_count}" -ne "${HOMOLOGY_EXPECTED_SCAFFOLDS}" ]; then
  echo "ERROR: expected ${HOMOLOGY_EXPECTED_SCAFFOLDS} scaffolds, found ${scaffold_count}." >&2
  exit 3
fi
if [ "${gene_count}" -ne "${HOMOLOGY_EXPECTED_GENES}" ]; then
  echo "ERROR: expected ${HOMOLOGY_EXPECTED_GENES} genes, found ${gene_count}." >&2
  exit 4
fi

for script in \
  scripts/stage_scaffold_homology_inputs.py \
  scripts/extract_unplaced_scaffolds.py \
  scripts/prepare_scaffold_homology_queries.py \
  scripts/summarize_scaffold_homology_hits.py \
  scripts/make_scaffold_origin_transfer_bundle.py; do
  if [ ! -s "${script}" ]; then
    echo "ERROR: homology workflow helper is missing: ${script}" >&2
    exit 5
  fi
done

(
  if command -v module >/dev/null 2>&1; then
    module purge
    module load "${HOMOLOGY_CONDA_MODULE}"
  fi
  source activate "${HOMOLOGY_CONDA_ENV}"

  for command_name in blastn blastdbcmd diamond; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
      echo "ERROR: ${command_name} is unavailable in ${HOMOLOGY_CONDA_ENV}." >&2
      echo "Install with: sbatch scripts/install_scaffold_homology_env.slurm" >&2
      exit 6
    fi
  done

  echo "  BLASTn: $(command -v blastn)"
  blastn -version | sed -n '1p'
  echo "  DIAMOND: $(command -v diamond)"
  diamond version
  blastdbcmd -db "${HOMOLOGY_CORE_NT_DB}" -info | sed -n '1,3p'
  blastdbcmd -db "${HOMOLOGY_NR_CLUSTER_SEQ_DB}" -info | sed -n '1,3p'
)

if ! command -v snakemake >/dev/null 2>&1; then
  echo "ERROR: snakemake is unavailable in the controller environment." >&2
  exit 7
fi

echo "Focused scaffold homology preflight passed."
echo "  Candidate contrast: infected_diet33_vs50"
echo "  Unresolved scaffolds: ${scaffold_count}"
echo "  Candidate DEG genes: ${gene_count}"
echo "  BLASTn database: ${HOMOLOGY_CORE_NT_DB}"
echo "  DIAMOND database: ${HOMOLOGY_NR_CLUSTER_SEQ_DB}"
echo "  S. gregaria self-taxid excluded: 7010"
echo "  RNA-seq preprocessing dependencies: none"
echo "  R dependencies: none"
