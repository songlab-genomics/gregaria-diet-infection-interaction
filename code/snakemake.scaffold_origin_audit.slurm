#!/usr/bin/env bash
#SBATCH --job-name=scaffold_origin_smk
#SBATCH --output=slurm-scaffold-origin-%j.out
#SBATCH --error=slurm-scaffold-origin-%j.err
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

TARGET="${1:-scaffold_origin_audit_hpc}"
JOBS="${2:-20}"
SNAKEMAKE_ENV="${3:-myENV}"
SNAKEFILE="${4:-Snakefile.scaffold_origin_audit}"
LATENCY_WAIT="${LATENCY_WAIT:-120}"
CLUSTER_CONFIG="${CLUSTER_CONFIG:-$(pwd)/cluster.json}"
PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"

if [ -z "${ORIGIN_RUN_ID:-}" ]; then
  echo "ERROR: ORIGIN_RUN_ID must be exported before submission." >&2
  echo "New run: export ORIGIN_RUN_ID=scaffold_origin_\$(date +%Y%m%d_%H%M%S)" >&2
  echo "Resume: export the exact ID and ORIGIN_RESUME_EXISTING=true." >&2
  exit 2
fi

ORIGIN_DIR="${ORIGIN_DIR:-${PROJECT_ROOT}/output/runs/${ORIGIN_RUN_ID}}"
ORIGIN_RESUME_EXISTING="${ORIGIN_RESUME_EXISTING:-false}"
REFDIR="${REFDIR:-/data/songlab/maeva/gregaria-timecourse/reference}"
PROJECTDIR="${PROJECTDIR:-$(dirname "${REFDIR}")}"
ORIGIN_HOST_FASTA="${ORIGIN_HOST_FASTA:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.fna}"
ORIGIN_HOST_ANNOTATION="${ORIGIN_HOST_ANNOTATION:-${REFDIR}/GCF_023897955.1_iqSchGreg1.2_genomic.gtf}"
ORIGIN_KRAKEN_DB="${ORIGIN_KRAKEN_DB:-/scratch/mtecher/kraken2/pluspf_20260626}"
ORIGIN_CONDA_MODULE="${ORIGIN_CONDA_MODULE:-mamba}"
ORIGIN_CONDA_ENV="${ORIGIN_CONDA_ENV:-/scratch/mtecher/conda-envs/metatranscriptomics-qc}"
ORIGIN_RUN_CROSS_SPECIES="${ORIGIN_RUN_CROSS_SPECIES:-true}"
ORIGIN_CROSS_SPECIES_REFERENCES="${ORIGIN_CROSS_SPECIES_REFERENCES:-piceifrons=/data/songlab/maeva/piceifrons-timecourse/reference/GCF_021461385.2_iqSchPice1.1_genomic.fna,americana=/data/songlab/maeva/americana-timecourse/reference/GCF_021461395.2_iqSchAmer2.1_genomic.fna,serialis_cubense=/data/songlab/maeva/cubense-timecourse/reference/GCF_023864345.2_iqSchSeri2.2_genomic.fna}"
ORIGIN_MINIMAP2_CMD="${ORIGIN_MINIMAP2_CMD:-minimap2}"
ORIGIN_CANDIDATE_LIST="${ORIGIN_CANDIDATE_LIST:-config/scaffold_origin/unified_unplaced_deg_scaffolds_20260729.txt}"
ORIGIN_CANDIDATE_GENES="${ORIGIN_CANDIDATE_GENES:-config/scaffold_origin/unified_unplaced_deg_genes_20260729.tsv}"
ORIGIN_CANDIDATE_PROVENANCE="${ORIGIN_CANDIDATE_PROVENANCE:-config/scaffold_origin/unified_unplaced_deg_scaffold_provenance_20260729.tsv}"
ORIGIN_NCBI_FCS_REPORT="${ORIGIN_NCBI_FCS_REPORT:-config/scaffold_origin/GCF_023897955.1_iqSchGreg1.2_fcs_report_20230323.txt}"
CLUSTER_LOG_DIR="${CLUSTER_LOG_DIR:-${PROJECT_ROOT}/logs/slurm/${ORIGIN_RUN_ID}}"

if [ -e "${ORIGIN_DIR}" ]; then
  case "${ORIGIN_RESUME_EXISTING}" in
  1|true|TRUE|yes|YES|y|Y|on|ON)
    echo "Resuming scaffold-origin run: ${ORIGIN_DIR}"
    echo "Completed Snakemake outputs will be retained."
    ;;
  *)
    echo "ERROR: run directory already exists: ${ORIGIN_DIR}" >&2
    echo "Use a new ORIGIN_RUN_ID, or resume only after checking scheduler state." >&2
    exit 2
    ;;
  esac
fi

export PROJECT_ROOT ORIGIN_RUN_ID ORIGIN_DIR ORIGIN_RESUME_EXISTING
export REFDIR PROJECTDIR ORIGIN_HOST_FASTA ORIGIN_HOST_ANNOTATION
export ORIGIN_KRAKEN_DB ORIGIN_CONDA_MODULE ORIGIN_CONDA_ENV
export ORIGIN_RUN_CROSS_SPECIES ORIGIN_CROSS_SPECIES_REFERENCES
export ORIGIN_MINIMAP2_CMD
export ORIGIN_CANDIDATE_LIST ORIGIN_CANDIDATE_GENES
export ORIGIN_CANDIDATE_PROVENANCE ORIGIN_NCBI_FCS_REPORT
export LOCUST_SPECIES="gregaria"
export LOCUST_TISSUES="fatbody"

mkdir -p "${CLUSTER_LOG_DIR}"

if command -v module >/dev/null 2>&1; then
  module purge
  module load mamba
fi
source activate "${SNAKEMAKE_ENV}"

echo "Run ID: ${ORIGIN_RUN_ID}"
echo "Target: ${TARGET}"
echo "Output: ${ORIGIN_DIR}"
echo "Candidate set: ${ORIGIN_CANDIDATE_LIST}"
echo "Host reference: ${ORIGIN_HOST_FASTA}"
echo "Kraken database: ${ORIGIN_KRAKEN_DB}"
echo "Cross-Schistocerca alignment enabled: ${ORIGIN_RUN_CROSS_SPECIES}"
echo "RNA-seq preprocessing rules included: false"
echo "R analysis required on HPC: false"

bash pilot/check_scaffold_origin_audit_hpc.sh

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
