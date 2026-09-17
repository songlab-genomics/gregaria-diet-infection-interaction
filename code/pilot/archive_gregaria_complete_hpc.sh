#!/usr/bin/env bash

# Move the complete authoritative Gregaria analysis from scratch to durable
# Song Lab storage. Copy and verification are non-destructive. Whole-project
# scratch removal is a separate, strongly gated final action.

set -euo pipefail

export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin${PATH:+:${PATH}}"

MODE="${1:-audit}"
PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
DEST_ROOT="${DEST_ROOT:-/data/songlab/maeva/gregaria-diet-infection-interaction}"
RAW_READS_DIR="${RAW_READS_DIR:-/data/songlab/sequencing_data/RNAseq/mehreen}"
HOST_REFERENCE_FNA="${HOST_REFERENCE_FNA:-/data/songlab/maeva/gregaria-timecourse/reference/GCF_023897955.1_iqSchGreg1.2_genomic.fna}"
HOST_REFERENCE_GTF="${HOST_REFERENCE_GTF:-/data/songlab/maeva/gregaria-timecourse/reference/GCF_023897955.1_iqSchGreg1.2_genomic.gtf}"
REQUIRE_CORRECTED_RUN_ON_DATA="${REQUIRE_CORRECTED_RUN_ON_DATA:-false}"

CANONICAL_RUN="host_pathogen_dual_20260727_000952"
ADDON_RUN="host_pathogen_dual_1044_20260808_110441"
CORRECTED_RUN="host_pathogen_dual_corrected45_20260809-153146"
ORIGIN_RUN="scaffold_origin_20260729_030308"
HOMOLOGY_RUN="scaffold_homology_20260803_042409"

# These are the only full HPC runs retained in the durable archive.
ACTIVE_SOURCE_RUNS=(
  "${CANONICAL_RUN}"
  "${ADDON_RUN}"
  "${ORIGIN_RUN}"
  "${HOMOLOGY_RUN}"
)

# The corrected integration was made locally. It must be transferred directly
# to DEST_ROOT before complete verification, but it is not required on scratch.
DESTINATION_ONLY_RUNS=(
  "${CORRECTED_RUN}"
)

# These partial or superseded runs are neither copied nor retained. Most were
# already removed by archive_gregaria_publication_hpc.sh.
OBSOLETE_RUNS=(
  host_pathogen_dual_20260726_005105
  host_pathogen_dual_20260726_191843
  host_pathogen_dual_20260726_193749
  host_pathogen_dual_20260726_214209
  host_pathogen_dual_20260727_002338
  host_pathogen_dual_20260727_010428
  fatbody_pilot_20260713_031027
  11-local-transfer
)

# Keep biological and workflow outputs, but omit replaceable bookkeeping,
# scheduler logs, editor files, and interrupted temporary files.
COMMON_EXCLUDES=(
  --exclude='.snakemake/'
  --exclude='__pycache__/'
  --exclude='*.py[cod]'
  --exclude='.DS_Store'
  --exclude='slurm*.out'
  --exclude='slurm*.err'
  --exclude='*.sh-*.out'
  --exclude='*.sh-*.err'
  --exclude='*.partial'
  --exclude='*.tmp'
  --exclude='*_STARtmp/'
)

PROJECT_EXCLUDES=(
  "${COMMON_EXCLUDES[@]}"
  --exclude='/.git/'
  --exclude='/.Rproj.user/'
  --exclude='/logs/'
  --exclude='/output/runs/'
)

require_commands() {
  local command_name missing=0
  for command_name in awk df du find grep rsync sha256sum sort tar wc; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
      echo "Missing required command: ${command_name}" >&2
      missing=1
    fi
  done
  [[ "${missing}" -eq 0 ]] || exit 127
}

validate_safe_roots() {
  if [[ "${PROJECT_ROOT}" != "/scratch/mtecher/gregaria-diet-infection-interaction" ]]; then
    echo "Refusing unexpected scratch root: ${PROJECT_ROOT}" >&2
    return 1
  fi
  if [[ "${DEST_ROOT}" != "/data/songlab/maeva/gregaria-diet-infection-interaction" ]]; then
    echo "Refusing unexpected durable root: ${DEST_ROOT}" >&2
    return 1
  fi
  if [[ "${PROJECT_ROOT}" == "${DEST_ROOT}" ]]; then
    echo "Source and destination cannot be the same path." >&2
    return 1
  fi
}

array_contains() {
  local query="$1"
  shift
  local item
  for item in "$@"; do
    [[ "${item}" == "${query}" ]] && return 0
  done
  return 1
}

count_nonempty() {
  local directory="$1"
  local pattern="$2"
  if [[ ! -d "${directory}" ]]; then
    printf '0\n'
    return
  fi
  find "${directory}" -type f -name "${pattern}" -size +0c | wc -l | awk '{print $1}'
}

assert_count() {
  local label="$1"
  local observed="$2"
  local expected="$3"
  if [[ "${observed}" -ne "${expected}" ]]; then
    echo "FAIL  ${label}: observed ${observed}, expected ${expected}" >&2
    return 1
  fi
  printf 'PASS  %-58s %s\n' "${label}" "${observed}"
}

require_nonempty_file() {
  local file="$1"
  if [[ ! -s "${file}" ]]; then
    echo "Missing or empty required file: ${file}" >&2
    return 1
  fi
}

classify_source_runs() {
  local runs_root="${PROJECT_ROOT}/output/runs"
  local run_path run_name unknown=0

  [[ -d "${runs_root}" ]] || {
    echo "Missing source run directory: ${runs_root}" >&2
    return 1
  }

  echo "=== Source run classification ==="
  for run_path in "${runs_root}"/*; do
    [[ -e "${run_path}" ]] || continue
    run_name=$(basename "${run_path}")
    if array_contains "${run_name}" "${ACTIVE_SOURCE_RUNS[@]}"; then
      printf 'RETAIN    %-55s %s\n' "${run_name}" "$(du -sh "${run_path}" | awk '{print $1}')"
    elif [[ "${run_name}" == "${CORRECTED_RUN}" ]]; then
      printf 'RETAIN    %-55s %s (source copy available)\n' "${run_name}" "$(du -sh "${run_path}" | awk '{print $1}')"
    elif array_contains "${run_name}" "${OBSOLETE_RUNS[@]}"; then
      printf 'OBSOLETE  %-55s %s\n' "${run_name}" "$(du -sh "${run_path}" | awk '{print $1}')"
    else
      printf 'UNKNOWN   %s\n' "${run_name}" >&2
      unknown=1
    fi
  done

  local active
  for active in "${ACTIVE_SOURCE_RUNS[@]}"; do
    if [[ ! -d "${runs_root}/${active}" ]]; then
      echo "Missing authoritative source run: ${runs_root}/${active}" >&2
      unknown=1
    fi
  done

  if [[ "${unknown}" -ne 0 ]]; then
    echo "Stop: classify every UNKNOWN run before copying or deleting scratch." >&2
    return 1
  fi
}

classify_destination_runs() {
  local runs_root="${DEST_ROOT}/output/runs"
  local run_path run_name invalid=0

  [[ -d "${runs_root}" ]] || {
    echo "Missing destination run directory: ${runs_root}" >&2
    return 1
  }

  echo "=== Durable run classification ==="
  for run_path in "${runs_root}"/*; do
    [[ -e "${run_path}" ]] || continue
    run_name=$(basename "${run_path}")
    if array_contains "${run_name}" "${ACTIVE_SOURCE_RUNS[@]}" || \
       array_contains "${run_name}" "${DESTINATION_ONLY_RUNS[@]}"; then
      printf 'RETAIN    %-55s %s\n' "${run_name}" "$(du -sh "${run_path}" | awk '{print $1}')"
    elif array_contains "${run_name}" "${OBSOLETE_RUNS[@]}"; then
      echo "OBSOLETE destination run must be reviewed and removed: ${run_path}" >&2
      invalid=1
    else
      echo "UNKNOWN destination run must be classified: ${run_path}" >&2
      invalid=1
    fi
  done
  [[ "${invalid}" -eq 0 ]]
}

validate_dual_run() {
  local root="$1"
  local samples="$2"
  local taxonomy_expected=$((samples * 3))
  local featurecounts_expected=$((samples * 5))
  local label
  label=$(basename "${root}")

  echo "=== Output counts: ${label} ==="
  assert_count "fastp JSON reports" \
    "$(count_nonempty "${root}/01-trimmed-fastp/qc" '*.fastp.json')" "${samples}"
  assert_count "trimmed paired FASTQ files" \
    "$(count_nonempty "${root}/01-trimmed-fastp" '*.trimmed.fastq.gz')" "$((samples * 2))"
  assert_count "host-only coordinate-sorted BAMs" \
    "$(count_nonempty "${root}/02-host-only-star" '*_Aligned.sortedByCoord.out.bam')" "${samples}"
  assert_count "host-only BAM indexes" \
    "$(count_nonempty "${root}/02-host-only-star" '*_Aligned.sortedByCoord.out.bam.csi')" "${samples}"
  assert_count "competitive coordinate-sorted BAMs" \
    "$(count_nonempty "${root}/02-competitive-star" '*_Aligned.sortedByCoord.out.bam')" "${samples}"
  assert_count "competitive BAM indexes" \
    "$(count_nonempty "${root}/02-competitive-star" '*_Aligned.sortedByCoord.out.bam.csi')" "${samples}"
  assert_count "host-unmapped paired FASTQ files" \
    "$(count_nonempty "${root}/02-host-only-star/unmapped" '*_host_unmapped_R*.fastq.gz')" "$((samples * 2))"
  assert_count "competitive-unmapped paired FASTQ files" \
    "$(count_nonempty "${root}/02-competitive-star/unmapped" '*_competitive_unmapped_R*.fastq.gz')" "$((samples * 2))"
  assert_count "featureCounts per-sample tables" \
    "$(count_nonempty "${root}/03-featurecounts" '*.featureCounts.txt')" "${featurecounts_expected}"
  assert_count "featureCounts summaries across five definitions" \
    "$(count_nonempty "${root}/03-featurecounts" '*.featureCounts.txt.summary')" "${featurecounts_expected}"
  assert_count "merged count matrices" \
    "$(count_nonempty "${root}/04-count-matrices" '*.tsv')" 5
  assert_count "top-level mapping comparison tables" \
    "$(find "${root}/05-mapping-comparison" -maxdepth 1 -type f -name '*.tsv' -size +0c 2>/dev/null | wc -l | awk '{print $1}')" 4
  assert_count "Kraken2 reports" \
    "$(count_nonempty "${root}/08-taxonomy" '*.kraken2.report')" "${taxonomy_expected}"
  assert_count "compressed Kraken2 classifications" \
    "$(count_nonempty "${root}/08-taxonomy" '*.kraken2.output.gz')" "${taxonomy_expected}"
  assert_count "Bracken family tables" \
    "$(count_nonempty "${root}/08-taxonomy" '*.bracken.family.tsv')" "${taxonomy_expected}"
  require_nonempty_file "${root}/11-local-transfer/transfer_manifest.tsv"
  require_nonempty_file "${root}/11-local-transfer/fatbody_hpc_results_for_local.tar.gz"
  tar -tzf "${root}/11-local-transfer/fatbody_hpc_results_for_local.tar.gz" >/dev/null
  echo "PASS  transfer archive is readable"
}

validate_audit_runs() {
  local root="$1"
  require_nonempty_file "${root}/output/runs/${ORIGIN_RUN}/06-local-transfer/transfer_manifest.tsv"
  require_nonempty_file "${root}/output/runs/${ORIGIN_RUN}/06-local-transfer/scaffold_origin_audit_for_local.tar.gz"
  tar -tzf "${root}/output/runs/${ORIGIN_RUN}/06-local-transfer/scaffold_origin_audit_for_local.tar.gz" >/dev/null
  require_nonempty_file "${root}/output/runs/${HOMOLOGY_RUN}/05-local-transfer/transfer_manifest.tsv"
  require_nonempty_file "${root}/output/runs/${HOMOLOGY_RUN}/05-local-transfer/scaffold_homology_for_local.tar.gz"
  tar -tzf "${root}/output/runs/${HOMOLOGY_RUN}/05-local-transfer/scaffold_homology_for_local.tar.gz" >/dev/null
  echo "PASS  scaffold-origin and homology transfer archives are readable"
}

validate_corrected_run() {
  local root="${DEST_ROOT}/output/runs/${CORRECTED_RUN}"
  local file sample_columns relative_output expected_hash observed_hash
  local matrices=0

  if [[ ! -d "${root}" ]]; then
    if [[ "${REQUIRE_CORRECTED_RUN_ON_DATA}" == "true" ]]; then
      echo "Missing corrected 45-sample integration at durable destination: ${root}" >&2
      return 1
    fi
    echo "SKIP  corrected 45-sample integration is retained as a local derivative."
    echo "      Its complete 44-sample source run and sample-1044 add-on are verified on durable storage."
    return 0
  fi

  echo "=== Corrected 45-sample integration ==="
  while IFS= read -r -d '' file; do
    sample_columns=$(awk -F '\t' 'NR == 1 {print NF - 1; exit}' "${file}")
    assert_count "$(basename "$(dirname "${file}")")/$(basename "${file}") sample columns" \
      "${sample_columns}" 45
    matrices=$((matrices + 1))
  done < <(find "${root}/04-count-matrices" -type f -name '*.tsv' -print0)
  assert_count "corrected host count matrices" "${matrices}" 4

  require_nonempty_file "${root}/matrix_integration_audit.tsv"
  require_nonempty_file "${root}/sample_assignment_audit.tsv"
  require_nonempty_file "${root}/run_provenance.tsv"

  while IFS=$'\t' read -r relative_output expected_hash; do
    file="${DEST_ROOT}/${relative_output}"
    require_nonempty_file "${file}"
    observed_hash=$(sha256sum "${file}" | awk '{print $1}')
    if [[ "${observed_hash}" != "${expected_hash}" ]]; then
      echo "Checksum mismatch for corrected matrix: ${file}" >&2
      echo "Expected: ${expected_hash}" >&2
      echo "Observed: ${observed_hash}" >&2
      return 1
    fi
  done < <(
    awk -F '\t' 'NR > 1 {print $5 "\t" $10}' \
      "${root}/matrix_integration_audit.tsv"
  )
  echo "PASS  corrected matrix SHA-256 checksums match the integration audit"

  grep -q $'^1024\tInfected\t33\t' "${root}/sample_assignment_audit.tsv"
  grep -q $'^1044\tInfected\t83\t' "${root}/sample_assignment_audit.tsv"
  echo "PASS  corrected assignments include 1024 = infected diet 33 and 1044 = infected diet 83"
}

validate_external_inputs() {
  [[ -d "${RAW_READS_DIR}" ]] || {
    echo "Missing durable raw-read directory: ${RAW_READS_DIR}" >&2
    return 1
  }
  require_nonempty_file "${HOST_REFERENCE_FNA}"
  require_nonempty_file "${HOST_REFERENCE_GTF}"
  echo "PASS  durable raw-read directory and host reference are present"
}

find_samtools() {
  local candidate
  if [[ -n "${SAMTOOLS_BIN:-}" && -x "${SAMTOOLS_BIN}" ]]; then
    printf '%s\n' "${SAMTOOLS_BIN}"
    return
  fi
  if command -v samtools >/dev/null 2>&1; then
    command -v samtools
    return
  fi
  for candidate in \
    /scratch/mtecher/conda-envs/metatranscriptomics-qc/bin/samtools \
    /scratch/mtecher/conda-envs/myENV/bin/samtools \
    /home/mtecher/.conda/envs/myENV/bin/samtools
  do
    if [[ -x "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return
    fi
  done
  return 1
}

validate_destination_bams() {
  local samtools_bin bam run count=0
  if ! samtools_bin=$(find_samtools); then
    echo "samtools is required for final BAM validation." >&2
    echo "Activate the mapping environment or export SAMTOOLS_BIN before verification." >&2
    return 1
  fi

  echo "=== Durable BAM readability and index checks ==="
  for run in "${CANONICAL_RUN}" "${ADDON_RUN}"; do
    while IFS= read -r -d '' bam; do
      "${samtools_bin}" quickcheck -v "${bam}"
      "${samtools_bin}" idxstats "${bam}" >/dev/null
      count=$((count + 1))
    done < <(
      find "${DEST_ROOT}/output/runs/${run}" \
        -type f -name '*_Aligned.sortedByCoord.out.bam' -print0
    )
  done
  assert_count "readable indexed BAM files" "${count}" 90
}

audit_complete_archive() {
  validate_safe_roots
  [[ -d "${PROJECT_ROOT}" ]] || {
    echo "Missing scratch project: ${PROJECT_ROOT}" >&2
    return 1
  }

  echo "Gregaria complete HPC archive audit"
  echo "Scratch source:    ${PROJECT_ROOT}"
  echo "Durable target:    ${DEST_ROOT}"
  echo "Durable raw reads: ${RAW_READS_DIR}"
  echo
  classify_source_runs
  echo
  validate_dual_run "${PROJECT_ROOT}/output/runs/${CANONICAL_RUN}" 44
  echo
  validate_dual_run "${PROJECT_ROOT}/output/runs/${ADDON_RUN}" 1
  echo
  validate_audit_runs "${PROJECT_ROOT}"
  echo
  echo "=== Storage ==="
  du -sh "${PROJECT_ROOT}"
  df -h "${PROJECT_ROOT}" "${DEST_ROOT}"
  echo
  echo "AUDIT ONLY: no files were copied or deleted."
}

copy_non_run_project_files() {
  mkdir -p "${DEST_ROOT}"
  rsync -aH --partial --info=progress2 "${PROJECT_EXCLUDES[@]}" \
    "${PROJECT_ROOT}/" "${DEST_ROOT}/"
}

copy_full_run() {
  local run="$1"
  local source="${PROJECT_ROOT}/output/runs/${run}"
  local destination="${DEST_ROOT}/output/runs/${run}"
  mkdir -p "${destination}"
  rsync -aH --partial --info=progress2 "${COMMON_EXCLUDES[@]}" \
    "${source}/" "${destination}/"
}

copy_complete_archive() {
  audit_complete_archive

  echo
  echo "=== Copying non-run project files ==="
  copy_non_run_project_files

  local run
  for run in "${ACTIVE_SOURCE_RUNS[@]}"; do
    echo
    echo "=== Copying complete run: ${run} ==="
    copy_full_run "${run}"
  done
  if [[ -d "${PROJECT_ROOT}/output/runs/${CORRECTED_RUN}" ]]; then
    echo
    echo "=== Copying source-available corrected integration ==="
    copy_full_run "${CORRECTED_RUN}"
  else
    echo
    echo "NOTICE: corrected 45-sample integration remains a local derivative."
    echo "        Set REQUIRE_CORRECTED_RUN_ON_DATA=true only if it must also be retained on /data."
  fi

  printf 'copied_at\t%s\nsource_root\t%s\ndestination_root\t%s\narchive_scope\tcomplete_authoritative_runs\n' \
    "$(date --iso-8601=seconds)" "${PROJECT_ROOT}" "${DEST_ROOT}" \
    > "${DEST_ROOT}/complete_archive_provenance.tsv"

  echo
  echo "Copy finished. Scratch remains unchanged. Run '$0 verify'."
}

verify_pair() {
  local label="$1"
  local source="$2"
  local destination="$3"
  shift 3
  local differences

  [[ -e "${source}" ]] || {
    echo "Missing verification source: ${source}" >&2
    return 1
  }
  [[ -e "${destination}" ]] || {
    echo "Missing verification destination: ${destination}" >&2
    return 1
  }

  differences=$(rsync -rlHnci --no-times "$@" "${source}/" "${destination}/")
  if [[ -n "${differences}" ]]; then
    echo "FAIL  ${label}" >&2
    printf '%s\n' "${differences}" >&2
    return 1
  fi
  echo "PASS  ${label}"
}

verify_complete_archive() {
  validate_safe_roots
  classify_source_runs

  echo
  echo "=== Checksum comparison: non-run project files ==="
  verify_pair "project files outside output/runs" \
    "${PROJECT_ROOT}" "${DEST_ROOT}" "${PROJECT_EXCLUDES[@]}"

  echo
  echo "=== Checksum comparison: complete authoritative runs ==="
  local run
  for run in "${ACTIVE_SOURCE_RUNS[@]}"; do
    verify_pair "output/runs/${run}" \
      "${PROJECT_ROOT}/output/runs/${run}" \
      "${DEST_ROOT}/output/runs/${run}" \
      "${COMMON_EXCLUDES[@]}"
  done

  echo
  validate_dual_run "${DEST_ROOT}/output/runs/${CANONICAL_RUN}" 44
  echo
  validate_dual_run "${DEST_ROOT}/output/runs/${ADDON_RUN}" 1
  echo
  validate_audit_runs "${DEST_ROOT}"
  echo
  validate_corrected_run
  echo
  validate_external_inputs
  echo
  classify_destination_runs
  echo
  validate_destination_bams

  printf 'verified_at\t%s\nsource_root\t%s\ndestination_root\t%s\nstatus\tPASS\n' \
    "$(date --iso-8601=seconds)" "${PROJECT_ROOT}" "${DEST_ROOT}" \
    > "${DEST_ROOT}/complete_archive_verification.tsv"

  echo
  echo "VERIFIED: complete authoritative HPC outputs are present on durable storage."
  echo "Scratch remains unchanged."
}

preview_complete_cleanup() {
  audit_complete_archive
  echo
  echo "=== Complete scratch-removal preview ==="
  du -sh "${PROJECT_ROOT}"
  echo "The entire directory below would be removed only after full verification:"
  echo "  ${PROJECT_ROOT}"
  echo
  echo "The Kraken and homology databases outside this project root are not removed:"
  echo "  /scratch/mtecher/kraken2"
  echo "  /scratch/mtecher/scaffold_homology_databases"
  echo
  echo "PREVIEW ONLY: nothing was deleted."
}

cleanup_complete_archive() {
  if [[ "${CONFIRM_COMPLETE_SCRATCH_REMOVAL:-}" != "DELETE_VERIFIED_GREGARIA_PROJECT" ]]; then
    echo "Complete scratch removal was not authorized." >&2
    echo "First run: $0 cleanup-preview" >&2
    echo "Then set CONFIRM_COMPLETE_SCRATCH_REMOVAL=DELETE_VERIFIED_GREGARIA_PROJECT" >&2
    return 2
  fi

  verify_complete_archive

  printf 'scratch_removal_started_at\t%s\nsource_root\t%s\n' \
    "$(date --iso-8601=seconds)" "${PROJECT_ROOT}" \
    > "${DEST_ROOT}/scratch_removal_provenance.tsv"

  echo
  echo "Removing fully verified scratch project: ${PROJECT_ROOT}"
  rm -rf -- "${PROJECT_ROOT}"
  if [[ -e "${PROJECT_ROOT}" ]]; then
    echo "Scratch project still exists after removal attempt: ${PROJECT_ROOT}" >&2
    return 1
  fi

  printf 'scratch_removal_completed_at\t%s\nsource_root\t%s\nstatus\tREMOVED\n' \
    "$(date --iso-8601=seconds)" "${PROJECT_ROOT}" \
    >> "${DEST_ROOT}/scratch_removal_provenance.tsv"
  echo "Complete: the verified project was removed from scratch."
  echo "Durable archive: ${DEST_ROOT}"
}

case "${MODE}" in
  audit)
    require_commands
    audit_complete_archive
    ;;
  copy)
    require_commands
    copy_complete_archive
    ;;
  verify)
    require_commands
    verify_complete_archive
    ;;
  cleanup-preview)
    require_commands
    preview_complete_cleanup
    ;;
  cleanup)
    require_commands
    cleanup_complete_archive
    ;;
  *)
    echo "Usage: $0 {audit|copy|verify|cleanup-preview|cleanup}" >&2
    exit 2
    ;;
esac
