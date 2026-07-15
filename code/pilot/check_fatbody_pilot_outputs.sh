#!/usr/bin/env bash
set -euo pipefail

# Check whether a fat body pilot run finished all expression outputs.
# Run from the cluster after the SLURM job exits:
#   cd /scratch/mtecher/gregaria-diet-infection-interaction/code
#   bash pilot/check_fatbody_pilot_outputs.sh

PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
SAMPLE_TABLE="${SAMPLE_TABLE:-${PROJECT_ROOT}/code/pilot/fatbody_pilot_samples.tsv}"
TISSUE="${TISSUE:-fatbody}"

latest_run_id() {
  local candidate_file
  candidate_file="${TMPDIR:-/tmp}/fatbody_pilot_candidates.$$"
  : > "${candidate_file}"

  # Completed or partially completed output folders.
  if [ -d "${PROJECT_ROOT}/output/runs" ]; then
    find "${PROJECT_ROOT}/output/runs" -maxdepth 1 -type d -name "fatbody_pilot_*" \
      -exec basename {} \; >> "${candidate_file}" 2>/dev/null || true
  fi

  # Symlink folders are created before Snakemake submits trimming jobs, so they
  # are useful evidence when a run fails before any output folder is made.
  if [ -d "${PROJECT_ROOT}/raw_links" ]; then
    find "${PROJECT_ROOT}/raw_links" -maxdepth 1 -type d -name "fatbody_pilot_*" \
      -exec basename {} \; >> "${candidate_file}" 2>/dev/null || true
  fi

  # Newer pilot launches write worker logs in logs/slurm/{RUN_ID}.
  if [ -d "${PROJECT_ROOT}/logs/slurm" ]; then
    find "${PROJECT_ROOT}/logs/slurm" -maxdepth 1 -type d -name "fatbody_pilot_*" \
      -exec basename {} \; >> "${candidate_file}" 2>/dev/null || true
  fi

  # Older pilot launches wrote only controller logs in the code folder.
  if [ -d "${PROJECT_ROOT}/code" ]; then
    find "${PROJECT_ROOT}/code" -maxdepth 1 -type f -name "slurm-fatbody-pilot-*.out" \
      -exec awk '/^Run ID:/ {print $3}' {} \; >> "${candidate_file}" 2>/dev/null || true
  fi

  sort -u "${candidate_file}" | tail -n 1
  rm -f "${candidate_file}"
}

RUN_ID="${1:-${RUN_ID:-$(latest_run_id)}}"

if [ -z "${RUN_ID}" ]; then
  echo "ERROR: no fatbody_pilot_* run evidence found under ${PROJECT_ROOT}" >&2
  echo "Look for controller logs with:" >&2
  echo "  ls -lh ${PROJECT_ROOT}/code/slurm-fatbody-pilot-*.out ${PROJECT_ROOT}/code/slurm-fatbody-pilot-*.err" >&2
  exit 1
fi

WORKDIR="${WORKDIR:-${PROJECT_ROOT}/output/runs/${RUN_ID}}"
LINK_DATADIR="${LINK_DATADIR:-${PROJECT_ROOT}/raw_links/${RUN_ID}}"
LOG_DIR="${LOG_DIR:-${PROJECT_ROOT}/logs/slurm/${RUN_ID}}"

failures=0
warnings=0

note_missing() {
  echo "MISSING: $1"
  failures=$((failures + 1))
}

note_warn() {
  echo "WARNING: $1"
  warnings=$((warnings + 1))
}

check_file() {
  local path="$1"
  if [ ! -s "${path}" ]; then
    note_missing "${path}"
  fi
}

echo "Fat body pilot output check"
echo "Run ID: ${RUN_ID}"
echo "Project root: ${PROJECT_ROOT}"
echo "Workdir: ${WORKDIR}"
echo "Raw-link dir: ${LINK_DATADIR}"
echo "Sample table: ${SAMPLE_TABLE}"
echo

[ -d "${WORKDIR}" ] || note_missing "${WORKDIR}"
[ -d "${LINK_DATADIR}" ] || note_missing "${LINK_DATADIR}"
[ -s "${SAMPLE_TABLE}" ] || note_missing "${SAMPLE_TABLE}"

if [ ! -s "${SAMPLE_TABLE}" ]; then
  echo
  echo "Stopping early because the sample table is missing."
  exit 1
fi

sample_count=0
while IFS=$'\t' read -r sample_id raw_prefix tissue treatment diet note; do
  if [ "${sample_id}" = "sample_id" ]; then
    continue
  fi
  if [ "${tissue}" != "${TISSUE}" ]; then
    continue
  fi

  sample_count=$((sample_count + 1))
  echo "Checking sample ${sample_id} (${treatment}, diet ${diet})"

  check_file "${LINK_DATADIR}/00-reads-${TISSUE}/${sample_id}_MERGE_1.fastq.gz"
  check_file "${LINK_DATADIR}/00-reads-${TISSUE}/${sample_id}_MERGE_2.fastq.gz"

  check_file "${WORKDIR}/01-${TISSUE}-trimmed-fastp/TrimQC/${sample_id}_fastp.json"
  check_file "${WORKDIR}/01-${TISSUE}-trimmed-fastp/TrimQC/${sample_id}_fastp.html"
  check_file "${WORKDIR}/01-${TISSUE}-trimmed-fastp/${sample_id}_1.trimmed.fastq.gz"
  check_file "${WORKDIR}/01-${TISSUE}-trimmed-fastp/${sample_id}_2.trimmed.fastq.gz"

  check_file "${WORKDIR}/02-${TISSUE}-star/${sample_id}_Aligned.sortedByCoord.out.bam"
  check_file "${WORKDIR}/02-${TISSUE}-star/${sample_id}_Aligned.sortedByCoord.out.bam.csi"
  check_file "${WORKDIR}/02-${TISSUE}-star/${sample_id}_ReadsPerGene.out.tab"
  check_file "${WORKDIR}/03-${TISSUE}-DESeq2/${sample_id}_counts.txt"

  check_file "${WORKDIR}/04-${TISSUE}-featurecounts/${sample_id}_counts.txt"
  check_file "${WORKDIR}/03-${TISSUE}-DESeq2/${sample_id}_featurecounts.txt"
done < "${SAMPLE_TABLE}"

if [ "${sample_count}" -eq 0 ]; then
  note_missing "no ${TISSUE} samples found in ${SAMPLE_TABLE}"
fi

echo
echo "Output counts:"
find "${WORKDIR}/01-${TISSUE}-trimmed-fastp" -type f -name "*.trimmed.fastq.gz" 2>/dev/null | wc -l | awk '{print "  trimmed FASTQs:", $1}'
find "${WORKDIR}/02-${TISSUE}-star" -type f -name "*_Aligned.sortedByCoord.out.bam" 2>/dev/null | wc -l | awk '{print "  STAR BAMs:", $1}'
find "${WORKDIR}/03-${TISSUE}-DESeq2" -type f -name "*_counts.txt" 2>/dev/null | wc -l | awk '{print "  STAR gene-count files:", $1}'
find "${WORKDIR}/03-${TISSUE}-DESeq2" -type f -name "*_featurecounts.txt" 2>/dev/null | wc -l | awk '{print "  featureCounts gene-count files:", $1}'

echo
echo "Controller logs:"
if find "${PROJECT_ROOT}/code" -maxdepth 1 -type f \( -name "slurm-fatbody-pilot-*.out" -o -name "slurm-fatbody-pilot-*.err" \) -print -quit 2>/dev/null | grep -q .; then
  find "${PROJECT_ROOT}/code" -maxdepth 1 -type f \( -name "slurm-fatbody-pilot-*.out" -o -name "slurm-fatbody-pilot-*.err" \) \
    -printf "  %p\n" 2>/dev/null || true
else
  echo "  none found in ${PROJECT_ROOT}/code"
fi

echo
if [ -d "${LOG_DIR}" ]; then
  echo "Cluster worker logs: ${LOG_DIR}"
  log_hits="${TMPDIR:-/tmp}/fatbody_pilot_log_hits.$$"
  if grep -R -n -E "Error|ERROR|Traceback|command not found|exited with non-zero" "${LOG_DIR}" >"${log_hits}" 2>/dev/null; then
    note_warn "possible errors found in worker logs:"
    cat "${log_hits}"
  fi
  rm -f "${log_hits}"
else
  note_warn "worker log directory not found: ${LOG_DIR}; older pilots may have logs in the code/ folder"
fi

echo
if [ "${failures}" -eq 0 ]; then
  echo "PASS: all expected pilot expression outputs are present for ${sample_count} samples."
  if [ "${warnings}" -gt 0 ]; then
    echo "There were ${warnings} warning(s); inspect the messages above before scaling up."
  fi
else
  echo "FAIL: ${failures} expected file(s) are missing."
  echo "The run did not complete the full fastp -> STAR -> featureCounts pilot."
  exit 1
fi
