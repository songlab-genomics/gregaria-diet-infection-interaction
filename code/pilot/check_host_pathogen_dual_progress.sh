#!/usr/bin/env bash
set -euo pipefail

# Read-only progress report for one host-pathogen dual RNA-seq run.
# It counts complete, non-empty outputs but never submits, unlocks, or deletes.

PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
SAMPLE_TABLE="${DUAL_SAMPLE_TABLE:-${PROJECT_ROOT}/code/config/fatbody_dual_rnaseq_samples.tsv}"
RUN_ID="${1:-${DUAL_RUN_ID:-}}"

if [ -z "${RUN_ID}" ]; then
  echo "ERROR: provide the run ID as the first argument or export DUAL_RUN_ID." >&2
  echo "Example: bash pilot/check_host_pathogen_dual_progress.sh host_pathogen_dual_20260727_000952" >&2
  exit 2
fi

RUN_DIR="${DUAL_DIR:-${PROJECT_ROOT}/output/runs/${RUN_ID}}"
LOG_DIR="${CLUSTER_LOG_DIR:-${PROJECT_ROOT}/logs/slurm/${RUN_ID}}"

if [ ! -d "${RUN_DIR}" ]; then
  echo "ERROR: run directory not found: ${RUN_DIR}" >&2
  exit 2
fi
if [ ! -s "${SAMPLE_TABLE}" ]; then
  echo "ERROR: sample table not found: ${SAMPLE_TABLE}" >&2
  exit 2
fi

SAMPLES=$(awk -F '\t' 'NR > 1 && $1 != "" {n++} END {print n + 0}' "${SAMPLE_TABLE}")
KRAKEN_SOURCES="${DUAL_KRAKEN_SOURCES:-all_trimmed,host_unmapped,competitive_unmapped}"
SOURCE_COUNT=$(awk -F ',' '{print NF}' <<< "${KRAKEN_SOURCES}")
TAXONOMY_EXPECTED=$((SAMPLES * SOURCE_COUNT))
BRACKEN_ENABLED="${DUAL_RUN_BRACKEN:-true}"

count_nonempty() {
  local directory="$1"
  local pattern="$2"
  if [ ! -d "${directory}" ]; then
    printf '0\n'
    return
  fi
  find "${directory}" -type f -name "${pattern}" -size +0c 2>/dev/null \
    | wc -l \
    | awk '{print $1}'
}

count_paths() {
  local count=0
  local path
  for path in "$@"; do
    if [ -s "${path}" ]; then
      count=$((count + 1))
    fi
  done
  printf '%s\n' "${count}"
}

report_stage() {
  local label="$1"
  local complete="$2"
  local expected="$3"
  local remaining=$((expected - complete))
  local state="RUNNING/WAITING"
  if [ "${remaining}" -lt 0 ]; then
    remaining=0
  fi
  if [ "${complete}" -ge "${expected}" ]; then
    state="COMPLETE"
  fi
  printf '%-48s %8s %8s %10s  %s\n' \
    "${label}" "${complete}" "${expected}" "${remaining}" "${state}"
}

host_index=0
if [ -s "${RUN_DIR}/00-reference/STAR-host-only/Genome" ] \
  && [ -s "${RUN_DIR}/00-reference/STAR-host-only/SA" ]; then
  host_index=1
fi

competitive_index=0
if [ -s "${RUN_DIR}/00-reference/STAR-competitive/Genome" ] \
  && [ -s "${RUN_DIR}/00-reference/STAR-competitive/SA" ]; then
  competitive_index=1
fi

matrix_count=$(count_paths \
  "${RUN_DIR}/04-count-matrices/host-only/host_exon_counts.tsv" \
  "${RUN_DIR}/04-count-matrices/host-only/host_transcript_exon_sensitivity_counts.tsv" \
  "${RUN_DIR}/04-count-matrices/competitive-host/host_exon_counts.tsv" \
  "${RUN_DIR}/04-count-matrices/competitive-host/host_transcript_exon_sensitivity_counts.tsv" \
  "${RUN_DIR}/04-count-matrices/competitive-fungus/metarhizium_exon_counts.tsv")

mapping_count=$(count_paths \
  "${RUN_DIR}/05-mapping-comparison/host_only_mapping_summary.tsv" \
  "${RUN_DIR}/05-mapping-comparison/competitive_mapping_summary.tsv" \
  "${RUN_DIR}/05-mapping-comparison/host_only_vs_competitive_sample_counts.tsv" \
  "${RUN_DIR}/05-mapping-comparison/host_only_vs_competitive_gene_counts.tsv")

transfer_count=$(count_paths \
  "${RUN_DIR}/11-local-transfer/transfer_manifest.tsv" \
  "${RUN_DIR}/11-local-transfer/fatbody_hpc_results_for_local.tar.gz")

echo "Host-pathogen dual RNA-seq progress"
echo "Run ID: ${RUN_ID}"
echo "Run directory: ${RUN_DIR}"
echo "Samples: ${SAMPLES}"
echo "Kraken sources: ${KRAKEN_SOURCES}"
echo
printf '%-48s %8s %8s %10s  %s\n' "Stage" "Done" "Expected" "Remaining" "Status"
printf '%-48s %8s %8s %10s  %s\n' "-----" "----" "--------" "---------" "------"
report_stage "fastp sample reports" \
  "$(count_nonempty "${RUN_DIR}/01-trimmed-fastp/qc" "*.fastp.json")" "${SAMPLES}"
report_stage "STAR host-only index" "${host_index}" 1
report_stage "STAR competitive index" "${competitive_index}" 1
report_stage "STAR host-only BAMs" \
  "$(count_nonempty "${RUN_DIR}/02-host-only-star" "*_Aligned.sortedByCoord.out.bam")" "${SAMPLES}"
report_stage "STAR competitive BAMs" \
  "$(count_nonempty "${RUN_DIR}/02-competitive-star" "*_Aligned.sortedByCoord.out.bam")" "${SAMPLES}"
report_stage "featureCounts: host-only exon" \
  "$(count_nonempty "${RUN_DIR}/03-featurecounts/host-only-exon" "*.featureCounts.txt.summary")" "${SAMPLES}"
report_stage "featureCounts: host-only transcript+exon" \
  "$(count_nonempty "${RUN_DIR}/03-featurecounts/host-only-transcript-exon" "*.featureCounts.txt.summary")" "${SAMPLES}"
report_stage "featureCounts: competitive host exon" \
  "$(count_nonempty "${RUN_DIR}/03-featurecounts/competitive-host-exon" "*.featureCounts.txt.summary")" "${SAMPLES}"
report_stage "featureCounts: competitive host transcript+exon" \
  "$(count_nonempty "${RUN_DIR}/03-featurecounts/competitive-host-transcript-exon" "*.featureCounts.txt.summary")" "${SAMPLES}"
report_stage "featureCounts: competitive fungus exon" \
  "$(count_nonempty "${RUN_DIR}/03-featurecounts/competitive-fungus-exon" "*.featureCounts.txt.summary")" "${SAMPLES}"
report_stage "Merged count matrices" "${matrix_count}" 5
report_stage "Mapping comparison tables" "${mapping_count}" 4
report_stage "Kraken2 reports" \
  "$(count_nonempty "${RUN_DIR}/08-taxonomy" "*.kraken2.report")" "${TAXONOMY_EXPECTED}"
case "${BRACKEN_ENABLED}" in
1|true|TRUE|yes|YES|y|Y|on|ON)
  report_stage "Bracken family tables" \
    "$(count_nonempty "${RUN_DIR}/08-taxonomy" "*.bracken.family.tsv")" "${TAXONOMY_EXPECTED}"
  ;;
*)
  printf '%-48s %8s %8s %10s  %s\n' "Bracken family tables" "-" "-" "-" "DISABLED"
  ;;
esac
report_stage "Local-transfer manifest and archive" "${transfer_count}" 2

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
  echo "  squeue is unavailable on this host."
fi

echo
echo "Worker log directory: ${LOG_DIR}"
if [ -d "${LOG_DIR}" ]; then
  echo "Non-empty worker logs: $(count_nonempty "${LOG_DIR}" "*.out") out, $(count_nonempty "${LOG_DIR}" "*.err") err"
  echo "Recent non-empty error logs:"
  find "${LOG_DIR}" -type f -name "*.err" -size +0c -print 2>/dev/null \
    | sort \
    | tail -n 10 \
    | sed 's/^/  /' || true
else
  echo "  Not created yet."
fi

echo
echo "This report is read-only. Do not submit a second controller while any"
echo "controller or worker from this run is still active."
