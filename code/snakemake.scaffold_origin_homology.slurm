#!/usr/bin/env bash
#SBATCH --job-name=scaffold_homology_smk
#SBATCH --output=slurm-scaffold-homology-%j.out
#SBATCH --error=slurm-scaffold-homology-%j.err
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=6-23:59:00
#SBATCH --export=ALL

set -euo pipefail

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -s "${SLURM_SUBMIT_DIR}/cluster.json" ]; then
  cd "${SLURM_SUBMIT_DIR}"
else
  cd "$(dirname "${BASH_SOURCE[0]}")"
fi

TARGET="${1:-scaffold_origin_homology_hpc}"
JOBS="${2:-8}"
SNAKEMAKE_ENV="${3:-myENV}"
SNAKEFILE="${4:-Snakefile.scaffold_origin_homology}"
LATENCY_WAIT="${LATENCY_WAIT:-120}"
CLUSTER_CONFIG="${CLUSTER_CONFIG:-$(pwd)/cluster.json}"
PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"

if [ -z "${HOMOLOGY_RUN_ID:-}" ]; then
  echo "ERROR: HOMOLOGY_RUN_ID must be exported before submission." >&2
  echo "New run: export HOMOLOGY_RUN_ID=scaffold_homology_\$(date +%Y%m%d_%H%M%S)" >&2
  echo "Resume: export the exact ID and HOMOLOGY_RESUME_EXISTING=true." >&2
  exit 2
fi

HOMOLOGY_DIR="${HOMOLOGY_DIR:-${PROJECT_ROOT}/output/runs/${HOMOLOGY_RUN_ID}}"
HOMOLOGY_RESUME_EXISTING="${HOMOLOGY_RESUME_EXISTING:-false}"
REFDIR="${REFDIR:-/data/songlab/maeva/gregaria-timecourse/reference}"
HOMOLOGY_HOST_FASTA="${HOMOLOGY_HOST_FASTA:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.fna}"
HOMOLOGY_HOST_ANNOTATION="${HOMOLOGY_HOST_ANNOTATION:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.gtf}"
HOMOLOGY_HOST_PROTEIN_FASTA="${HOMOLOGY_HOST_PROTEIN_FASTA:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_protein.faa}"
HOMOLOGY_BLAST_DB_ROOT="${HOMOLOGY_BLAST_DB_ROOT:-/scratch/mtecher/scaffold_homology_databases/ncbi_blastdb}"
HOMOLOGY_CORE_NT_DB="${HOMOLOGY_CORE_NT_DB:-${HOMOLOGY_BLAST_DB_ROOT}/core_nt}"
HOMOLOGY_NR_CLUSTER_SEQ_DB="${HOMOLOGY_NR_CLUSTER_SEQ_DB:-${HOMOLOGY_BLAST_DB_ROOT}/nr_cluster_seq}"
HOMOLOGY_DATABASE_READY="${HOMOLOGY_DATABASE_READY:-${HOMOLOGY_BLAST_DB_ROOT}/download_complete.tsv}"
HOMOLOGY_TAXONOMY_NODES="${HOMOLOGY_TAXONOMY_NODES:-/scratch/mtecher/kraken2/pluspf_20260626/nodes.dmp}"
HOMOLOGY_TAXONOMY_NAMES="${HOMOLOGY_TAXONOMY_NAMES:-/scratch/mtecher/kraken2/pluspf_20260626/names.dmp}"
HOMOLOGY_CONDA_ENV="${HOMOLOGY_CONDA_ENV:-/scratch/mtecher/conda-envs/scaffold-homology}"
HOMOLOGY_CONDA_MODULE="${HOMOLOGY_CONDA_MODULE:-mamba}"
HOMOLOGY_CANDIDATE_SCAFFOLDS="${HOMOLOGY_CANDIDATE_SCAFFOLDS:-config/scaffold_homology/infected_diet33_vs50_unresolved_scaffolds_20260731.txt}"
HOMOLOGY_CANDIDATE_GENES="${HOMOLOGY_CANDIDATE_GENES:-config/scaffold_homology/infected_diet33_vs50_unresolved_genes_20260731.tsv}"
HOMOLOGY_CANDIDATE_PROVENANCE="${HOMOLOGY_CANDIDATE_PROVENANCE:-config/scaffold_homology/infected_diet33_vs50_unresolved_provenance_20260731.tsv}"
CLUSTER_LOG_DIR="${CLUSTER_LOG_DIR:-${PROJECT_ROOT}/logs/slurm/${HOMOLOGY_RUN_ID}}"

if [ -e "${HOMOLOGY_DIR}" ]; then
  case "${HOMOLOGY_RESUME_EXISTING}" in
  1|true|TRUE|yes|YES|y|Y|on|ON)
    echo "Resuming homology run: ${HOMOLOGY_DIR}"
    echo "Completed Snakemake outputs will be retained."
    ;;
  *)
    echo "ERROR: run directory already exists: ${HOMOLOGY_DIR}" >&2
    echo "Use a new HOMOLOGY_RUN_ID, or resume only after checking scheduler state." >&2
    exit 2
    ;;
  esac
fi

export PROJECT_ROOT HOMOLOGY_RUN_ID HOMOLOGY_DIR HOMOLOGY_RESUME_EXISTING
export REFDIR HOMOLOGY_HOST_FASTA HOMOLOGY_HOST_ANNOTATION
export HOMOLOGY_HOST_PROTEIN_FASTA HOMOLOGY_BLAST_DB_ROOT
export HOMOLOGY_CORE_NT_DB HOMOLOGY_NR_CLUSTER_SEQ_DB
export HOMOLOGY_DATABASE_READY
export HOMOLOGY_TAXONOMY_NODES HOMOLOGY_TAXONOMY_NAMES
export HOMOLOGY_CONDA_ENV HOMOLOGY_CONDA_MODULE
export HOMOLOGY_CANDIDATE_SCAFFOLDS HOMOLOGY_CANDIDATE_GENES
export HOMOLOGY_CANDIDATE_PROVENANCE
export LOCUST_SPECIES="gregaria"
export LOCUST_TISSUES="fatbody"

mkdir -p "${CLUSTER_LOG_DIR}"

module purge
module load mamba
source activate "${SNAKEMAKE_ENV}"

echo "Run ID: ${HOMOLOGY_RUN_ID}"
echo "Target: ${TARGET}"
echo "Output: ${HOMOLOGY_DIR}"
echo "Contrast: infected_diet33_vs50"
echo "Target set: 948 unresolved DEG genes on 59 scaffolds"
echo "BLASTn database: ${HOMOLOGY_CORE_NT_DB}"
echo "DIAMOND database: ${HOMOLOGY_NR_CLUSTER_SEQ_DB}"
echo "RNA-seq preprocessing rules included: false"

bash pilot/check_scaffold_origin_homology_hpc.sh

snakemake -s "${SNAKEFILE}" "${TARGET}" \
  --dry-run \
  --rerun-incomplete \
  --rerun-triggers mtime

snakemake -s "${SNAKEFILE}" "${TARGET}" \
  --jobs "${JOBS}" \
  --latency-wait "${LATENCY_WAIT}" \
  --rerun-incomplete \
  --keep-going \
  --rerun-triggers mtime \
  --cluster-cancel scancel \
  --cluster-config "${CLUSTER_CONFIG}" \
  --cluster "sbatch --output=${CLUSTER_LOG_DIR}/%x-%j.out --error=${CLUSTER_LOG_DIR}/%x-%j.err -N {cluster.nodes} -n {cluster.ntasks} -c {cluster.cpus-per-task} --mem={cluster.mem} -t {cluster.time} -p {cluster.partition} --qos={cluster.qos} --export=ALL"

echo "Finished: $(date)"
