#!/usr/bin/env bash

# Preserve the compact, publication-relevant products on Song Lab storage.
# Audit, copy, and verify never delete anything. Cleanup is separately gated.

set -euo pipefail

# The archive helper needs only standard system tools. Prepend their canonical
# locations so a stale Conda/module environment cannot hide chmod, find, tar,
# rsync, or the other utilities used below.
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin${PATH:+:${PATH}}"

MODE="${1:-audit}"
PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
DEST_ROOT="${DEST_ROOT:-/data/songlab/maeva/gregaria-diet-infection-interaction}"

CANONICAL_RUN="host_pathogen_dual_20260727_000952"
ADDON_RUN="host_pathogen_dual_1044_20260808_110441"
CORRECTED_RUN="host_pathogen_dual_corrected45_20260809-153146"
ORIGIN_RUN="scaffold_origin_20260729_030308"
HOMOLOGY_RUN="scaffold_homology_20260803_042409"

# These runs were superseded by the complete 44-library run and 1044 add-on.
# Cleanup is restricted to this explicit allowlist.
OBSOLETE_RUNS=(
  host_pathogen_dual_20260726_005105
  host_pathogen_dual_20260726_191843
  host_pathogen_dual_20260726_193749
  host_pathogen_dual_20260726_214209
  fatbody_pilot_20260713_031027
)

# Preserve the exact code/configuration plus compact, checksummed result bundles.
# The corrected 45-sample integration is copied when it exists on the HPC.
PRESERVE_PATHS=(
  code
  data/metadata
  data/excluded_loci
  reference/metarhizium/GCF_000187425.2
  "output/runs/${CANONICAL_RUN}/11-local-transfer"
  "output/runs/${ADDON_RUN}/11-local-transfer"
  "output/runs/${CORRECTED_RUN}"
  "output/runs/${ORIGIN_RUN}/06-local-transfer"
  "output/runs/${HOMOLOGY_RUN}/05-local-transfer"
)

RSYNC_EXCLUDES=(
  --exclude='.snakemake/'
  --exclude='logs/'
  --exclude='slurm*.out'
  --exclude='slurm*.err'
  --exclude='*.sh-*.out'
  --exclude='*.sh-*.err'
)

require_commands() {
  local command_name
  local missing=0
  for command_name in awk du find rsync sort tar wc; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
      echo "Missing required command: ${command_name}" >&2
      missing=1
    fi
  done
  [[ "${missing}" -eq 0 ]] || exit 127
}

relative_file_count() {
  local path="$1"
  if [[ -f "${path}" ]]; then
    printf '1\n'
  elif [[ -d "${path}" ]]; then
    find "${path}" -type f \
      ! -name 'slurm*.out' ! -name 'slurm*.err' \
      ! -name '*.sh-*.out' ! -name '*.sh-*.err' | wc -l | awk '{print $1}'
  else
    printf '0\n'
  fi
}

show_preserve_inventory() {
  local relative source destination state

  printf '%-88s %10s %10s %s\n' "Publication archive item" "Files" "Size" "Destination"
  printf '%-88s %10s %10s %s\n' "------------------------" "-----" "----" "-----------"
  for relative in "${PRESERVE_PATHS[@]}"; do
    source="${PROJECT_ROOT}/${relative}"
    destination="${DEST_ROOT}/${relative}"
    if [[ -e "${source}" ]]; then
      state="missing"
      [[ -e "${destination}" ]] && state="present"
      printf '%-88s %10s %10s %s\n' \
        "${relative}" \
        "$(relative_file_count "${source}")" \
        "$(du -sh "${source}" | awk '{print $1}')" \
        "${state}"
    else
      printf '%-88s %10s %10s %s\n' "${relative}" "-" "-" "SOURCE MISSING"
    fi
  done
}

find_slurm_logs() {
  local print_action="${1:--print}"
  find "${PROJECT_ROOT}" \
    -path "${PROJECT_ROOT}/.git" -prune -o \
    -type f \( \
      -name 'slurm*.out' -o -name 'slurm*.err' -o \
      -name '*.sh-*.out' -o -name '*.sh-*.err' \
    \) "${print_action}"
}

show_cleanup_inventory() {
  local run path

  echo "=== Explicitly obsolete run directories ==="
  for run in "${OBSOLETE_RUNS[@]}"; do
    path="${PROJECT_ROOT}/output/runs/${run}"
    if [[ -d "${path}" ]]; then
      du -sh "${path}"
    else
      echo "MISSING  ${path}"
    fi
  done

  echo
  echo "=== SLURM output/error files ==="
  local log_count log_bytes
  log_count=$(find_slurm_logs | wc -l | awk '{print $1}')
  if [[ "${log_count}" -gt 0 ]]; then
    log_bytes=$(find_slurm_logs -print0 2>/dev/null \
      | du --files0-from=- -ch 2>/dev/null \
      | tail -n 1 \
      | awk '{print $1}')
  else
    log_bytes="0"
  fi
  echo "Files: ${log_count}"
  echo "Size:  ${log_bytes:-0}"
  find_slurm_logs | sort
}

validate_transfer_archives() {
  local root="$1"
  local archive
  local archives=(
    "${root}/output/runs/${CANONICAL_RUN}/11-local-transfer/fatbody_hpc_results_for_local.tar.gz"
    "${root}/output/runs/${ADDON_RUN}/11-local-transfer/fatbody_hpc_results_for_local.tar.gz"
    "${root}/output/runs/${ORIGIN_RUN}/06-local-transfer/scaffold_origin_audit_for_local.tar.gz"
    "${root}/output/runs/${HOMOLOGY_RUN}/05-local-transfer/scaffold_homology_for_local.tar.gz"
  )

  for archive in "${archives[@]}"; do
    if [[ ! -s "${archive}" ]]; then
      echo "Missing required transfer archive: ${archive}" >&2
      return 1
    fi
    tar -tzf "${archive}" >/dev/null
    echo "PASS  readable archive: ${archive}"
  done
}

audit_archive() {
  if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "Missing scratch project: ${PROJECT_ROOT}" >&2
    exit 1
  fi

  echo "Gregaria publication archive audit"
  echo "Scratch source: ${PROJECT_ROOT}"
  echo "Durable target: ${DEST_ROOT}"
  echo
  show_preserve_inventory
  echo
  validate_transfer_archives "${PROJECT_ROOT}"
  echo
  show_cleanup_inventory
  echo
  echo "AUDIT ONLY: no files were copied or deleted."
}

copy_item() {
  local relative="$1"
  local source="${PROJECT_ROOT}/${relative}"
  local destination="${DEST_ROOT}/${relative}"

  if [[ ! -e "${source}" ]]; then
    if [[ "${relative}" == "output/runs/${CORRECTED_RUN}" ]]; then
      echo "SKIP  optional corrected-45 integration is not present on this HPC."
      return
    fi
    echo "Missing required archive input: ${source}" >&2
    return 1
  fi

  mkdir -p "$(dirname "${destination}")"
  if [[ -d "${source}" ]]; then
    mkdir -p "${destination}"
    rsync -aH --partial --info=progress2 "${RSYNC_EXCLUDES[@]}" \
      "${source}/" "${destination}/"
  else
    rsync -aH --partial --info=progress2 "${source}" "${destination}"
  fi
}

copy_archive() {
  audit_archive
  mkdir -p "${DEST_ROOT}"

  local relative
  for relative in "${PRESERVE_PATHS[@]}"; do
    echo
    echo "=== Copying ${relative} ==="
    copy_item "${relative}"
  done

  printf 'archived_at\t%s\nsource_root\t%s\ndestination_root\t%s\n' \
    "$(date --iso-8601=seconds)" "${PROJECT_ROOT}" "${DEST_ROOT}" \
    > "${DEST_ROOT}/archive_provenance.tsv"

  echo
  echo "Copy finished. Run '$0 verify' before cleanup."
}

verify_item() {
  local relative="$1"
  local source="${PROJECT_ROOT}/${relative}"
  local destination="${DEST_ROOT}/${relative}"
  local differences

  if [[ ! -e "${source}" ]]; then
    if [[ "${relative}" == "output/runs/${CORRECTED_RUN}" ]]; then
      echo "SKIP  ${relative} was not present at the source."
      return
    fi
    echo "Missing required source: ${source}" >&2
    return 1
  fi
  if [[ ! -e "${destination}" ]]; then
    echo "Missing destination copy: ${destination}" >&2
    return 1
  fi

  if [[ -d "${source}" ]]; then
    # Scratch and project storage can assign different mtimes during staging.
    # Verify file content by checksum while ignoring timestamp-only drift.
    differences=$(rsync -aHnci --no-times "${RSYNC_EXCLUDES[@]}" \
      "${source}/" "${destination}/")
  else
    differences=$(rsync -aHnci --no-times "${source}" "${destination}")
  fi

  if [[ -n "${differences}" ]]; then
    echo "FAIL  ${relative}" >&2
    printf '%s\n' "${differences}" >&2
    return 1
  fi
  echo "PASS  ${relative}"
}

verify_archive() {
  if [[ ! -d "${DEST_ROOT}" ]]; then
    echo "Missing durable archive: ${DEST_ROOT}" >&2
    exit 1
  fi

  echo "=== Checksum-level source/destination comparison ==="
  local relative
  for relative in "${PRESERVE_PATHS[@]}"; do
    verify_item "${relative}"
  done

  echo
  echo "=== Durable transfer-archive integrity ==="
  validate_transfer_archives "${DEST_ROOT}"
  echo
  echo "VERIFIED: publication essentials are present in durable storage."
  echo "No scratch files were deleted."
}

cleanup_archive() {
  if [[ "${CONFIRM_SCRATCH_CLEANUP:-}" != "YES" ]]; then
    echo "Cleanup was not authorized." >&2
    echo "First run: $0 cleanup-preview" >&2
    echo "Then export CONFIRM_SCRATCH_CLEANUP=YES and run: $0 cleanup" >&2
    exit 2
  fi

  # Cleanup cannot start unless the compact durable copy passes verification.
  verify_archive

  local run path
  for run in "${OBSOLETE_RUNS[@]}"; do
    path="${PROJECT_ROOT}/output/runs/${run}"
    if [[ -d "${path}" ]]; then
      echo "Removing explicitly obsolete run: ${path}"
      rm -rf -- "${path}"
    fi
  done

  echo "Removing SLURM output/error files listed by cleanup-preview."
  local log_file log_count=0
  while IFS= read -r -d '' log_file; do
    rm -f -- "${log_file}"
    log_count=$((log_count + 1))
  done < <(find_slurm_logs -print0)
  echo "Removed ${log_count} SLURM output/error files."

  echo "Cleanup complete. Canonical, add-on, and audit run directories remain."
}

case "${MODE}" in
  audit)
    require_commands
    audit_archive
    ;;
  copy)
    require_commands
    copy_archive
    ;;
  verify)
    require_commands
    verify_archive
    ;;
  cleanup-preview)
    require_commands
    show_cleanup_inventory
    echo
    echo "PREVIEW ONLY: no files were deleted."
    ;;
  cleanup)
    require_commands
    cleanup_archive
    ;;
  *)
    echo "Usage: $0 {audit|copy|verify|cleanup-preview|cleanup}" >&2
    exit 2
    ;;
esac
