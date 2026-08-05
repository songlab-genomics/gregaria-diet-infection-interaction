#!/usr/bin/env bash
set -euo pipefail

# Read-only status report. It never submits, unlocks, or removes files.

PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
RUN_ID="${1:-${HOMOLOGY_RUN_ID:-}}"
if [ -z "${RUN_ID}" ]; then
  echo "ERROR: provide HOMOLOGY_RUN_ID as the first argument." >&2
  exit 2
fi

RUN_DIR="${HOMOLOGY_DIR:-${PROJECT_ROOT}/output/runs/${RUN_ID}}"
LOG_DIR="${CLUSTER_LOG_DIR:-${PROJECT_ROOT}/logs/slurm/${RUN_ID}}"
if [ ! -d "${RUN_DIR}" ]; then
  echo "ERROR: homology run directory not found: ${RUN_DIR}" >&2
  exit 2
fi

count_rows() {
  local path="$1"
  if [ ! -s "${path}" ]; then
    echo 0
  else
    awk 'NR > 1 {n++} END {print n + 0}' "${path}"
  fi
}

count_fasta() {
  local path="$1"
  if [ ! -s "${path}" ]; then
    echo 0
  else
    grep -c '^>' "${path}" || true
  fi
}

file_done() {
  if [ -e "$1" ]; then echo 1; else echo 0; fi
}

status_row() {
  local label="$1"
  local done="$2"
  local expected="$3"
  local state="RUNNING/WAITING"
  if [ "${done}" -ge "${expected}" ]; then state="COMPLETE"; fi
  printf '%-44s %8s %8s  %s\n' "${label}" "${done}" "${expected}" "${state}"
}

expected_scaffolds="${HOMOLOGY_EXPECTED_SCAFFOLDS:-59}"
expected_genes="${HOMOLOGY_EXPECTED_GENES:-948}"

echo "Focused scaffold BLAST/DIAMOND audit progress"
echo "Run ID: ${RUN_ID}"
echo "Run directory: ${RUN_DIR}"
echo
printf '%-44s %8s %8s  %s\n' "Stage" "Done" "Expected" "Status"
printf '%-44s %8s %8s  %s\n' "-----" "----" "--------" "------"
status_row "Staged candidate scaffolds" \
  "$(awk 'NF && $1 !~ /^#/ {n++} END {print n + 0}' "${RUN_DIR}/00-input/candidate_scaffolds.txt" 2>/dev/null || echo 0)" \
  "${expected_scaffolds}"
status_row "Staged candidate genes" \
  "$(count_rows "${RUN_DIR}/00-input/candidate_genes.tsv")" \
  "${expected_genes}"
status_row "Extracted scaffold FASTA records" \
  "$(count_fasta "${RUN_DIR}/01-sequences/candidate_scaffolds.fna")" \
  "${expected_scaffolds}"
status_row "Nucleotide gene queries" \
  "$(count_fasta "${RUN_DIR}/01-sequences/candidate_gene_regions.fna")" \
  "${expected_genes}"
status_row "Protein query file produced" \
  "$(file_done "${RUN_DIR}/01-sequences/candidate_proteins.faa")" 1
status_row "BLASTn core_nt search produced" \
  "$(file_done "${RUN_DIR}/02-blastn-core-nt/gene_regions_vs_core_nt.tsv")" 1
status_row "DIAMOND nr_cluster_seq search produced" \
  "$(file_done "${RUN_DIR}/03-diamond-nr-cluster-seq/proteins_vs_nr_cluster_seq.tsv")" 1
status_row "Gene homology summary rows" \
  "$(count_rows "${RUN_DIR}/04-evidence/gene_homology_summary.tsv")" \
  "${expected_genes}"
status_row "Scaffold homology summary rows" \
  "$(count_rows "${RUN_DIR}/04-evidence/scaffold_homology_summary.tsv")" \
  "${expected_scaffolds}"

transfer_done=0
if [ -s "${RUN_DIR}/05-local-transfer/transfer_manifest.tsv" ] \
  && [ -s "${RUN_DIR}/05-local-transfer/scaffold_homology_for_local.tar.gz" ]; then
  transfer_done=1
fi
status_row "Transfer manifest and archive" "${transfer_done}" 1

echo
echo "Current SLURM jobs for ${USER:-current user}:"
if command -v squeue >/dev/null 2>&1; then
  queue_output=$(squeue -h -u "${USER}" -o "%.18i %.32j %.10T %R" 2>&1 || true)
  if [ -n "${queue_output}" ]; then
    printf '%s\n' "${queue_output}" | sed 's/^/  /'
  else
    echo "  None."
  fi
else
  echo "  squeue is unavailable."
fi

echo
echo "Worker log directory: ${LOG_DIR}"
echo "This report is read-only."
