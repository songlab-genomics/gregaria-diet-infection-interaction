#!/usr/bin/env bash
set -euo pipefail

# Read-only status report. This script never submits, unlocks, or removes files.

PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
RUN_ID="${1:-${ORIGIN_RUN_ID:-}}"
if [ -z "${RUN_ID}" ]; then
  echo "ERROR: provide ORIGIN_RUN_ID as the first argument." >&2
  exit 2
fi

RUN_DIR="${ORIGIN_DIR:-${PROJECT_ROOT}/output/runs/${RUN_ID}}"
LOG_DIR="${CLUSTER_LOG_DIR:-${PROJECT_ROOT}/logs/slurm/${RUN_ID}}"
if [ ! -d "${RUN_DIR}" ]; then
  echo "ERROR: scaffold-origin run directory not found: ${RUN_DIR}" >&2
  exit 2
fi

count_data_rows() {
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

status_row() {
  local label="$1"
  local done="$2"
  local expected="$3"
  local state="RUNNING/WAITING"
  if [ "${done}" -ge "${expected}" ]; then
    state="COMPLETE"
  fi
  printf '%-42s %8s %8s  %s\n' "${label}" "${done}" "${expected}" "${state}"
}

expected=$(
  awk 'NF && $1 !~ /^#/ {n++} END {print n + 0}' \
    "${RUN_DIR}/00-input/candidate_scaffolds.txt" 2>/dev/null || echo 0
)
if [ "${expected}" -eq 0 ]; then
  expected="${ORIGIN_EXPECTED_SCAFFOLDS:-77}"
fi

echo "Scaffold-origin audit progress"
echo "Run ID: ${RUN_ID}"
echo "Run directory: ${RUN_DIR}"
echo
printf '%-42s %8s %8s  %s\n' "Stage" "Done" "Expected" "Status"
printf '%-42s %8s %8s  %s\n' "-----" "----" "--------" "------"
status_row "Staged candidate scaffolds" \
  "$(awk 'NF && $1 !~ /^#/ {n++} END {print n + 0}' "${RUN_DIR}/00-input/candidate_scaffolds.txt" 2>/dev/null || echo 0)" \
  "${expected}"
status_row "Extracted scaffold FASTA records" \
  "$(count_fasta "${RUN_DIR}/01-sequences/candidate_scaffolds.fna")" \
  "${expected}"
status_row "Sequence metric rows" \
  "$(count_data_rows "${RUN_DIR}/01-sequences/candidate_scaffold_sequence_metrics.tsv")" \
  "${expected}"
status_row "Annotation summary rows" \
  "$(count_data_rows "${RUN_DIR}/02-annotation/candidate_scaffold_annotation_summary.tsv")" \
  "${expected}"
status_row "Kraken scaffold calls" \
  "$(awk 'NF {n++} END {print n + 0}' "${RUN_DIR}/03-taxonomy/candidate_scaffolds.kraken2.output" 2>/dev/null || echo 0)" \
  "${expected}"
cross_pafs=0
if [ -d "${RUN_DIR}/04-cross-species" ]; then
  cross_pafs=$(
    find "${RUN_DIR}/04-cross-species" -maxdepth 1 -type f -name '*.paf' -size +0c 2>/dev/null \
      | wc -l \
      | awk '{print $1}'
  )
fi
if [[ "${ORIGIN_RUN_CROSS_SPECIES:-true}" =~ ^(0|false|FALSE|no|NO|off|OFF)$ ]]; then
  cross_expected=0
else
  cross_expected="${ORIGIN_EXPECTED_CROSS_SPECIES:-3}"
fi
status_row "Cross-species PAF files" "${cross_pafs}" "${cross_expected}"
status_row "Cross-species summary rows" \
  "$(count_data_rows "${RUN_DIR}/04-cross-species/cross_species_scaffold_homology.tsv")" \
  "${expected}"
status_row "Integrated evidence rows" \
  "$(count_data_rows "${RUN_DIR}/05-evidence/scaffold_origin_evidence.tsv")" \
  "${expected}"

transfer_done=0
if [ -s "${RUN_DIR}/06-local-transfer/transfer_manifest.tsv" ] \
  && [ -s "${RUN_DIR}/06-local-transfer/scaffold_origin_audit_for_local.tar.gz" ]; then
  transfer_done=1
fi
status_row "Transfer manifest and archive" "${transfer_done}" 1

echo
echo "Current SLURM jobs for ${USER:-current user}:"
if command -v squeue >/dev/null 2>&1; then
  if queue_output=$(squeue -h -u "${USER}" -o "%.18i %.32j %.10T %R" 2>&1); then
    if [ -n "${queue_output}" ]; then
      printf '%s\n' "${queue_output}" | sed 's/^/  /'
    else
      echo "  None."
    fi
  else
    echo "  Unable to query SLURM."
    printf '%s\n' "${queue_output}" | sed 's/^/  /'
  fi
else
  echo "  squeue is unavailable."
fi

echo
echo "Worker log directory: ${LOG_DIR}"
echo "This report is read-only."
