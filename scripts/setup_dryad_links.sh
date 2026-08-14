#!/usr/bin/env bash
set -euo pipefail

# Connect project-relative analysis paths to the separately deposited data.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
DEFAULT_DRYAD_DIR="/Users/maevatecher/Dropbox/3. Research Papers/Locusts/Ongoing/gregaria-diet-infection-interaction-dryad"
DRYAD_DIR="${GREGARIA_DIET_DRYAD_DIR:-$DEFAULT_DRYAD_DIR}"

LINK_PATHS=(
  "data/GO_Annotations"
  "data/raw_read_counts"
  "data/reference"
  "data/external"
  "data/legacy_csv"
  "data/scaffold_origin"
  "output/rmd_runs"
  "output/runs"
  "outputs"
)

if [[ ! -d "$DRYAD_DIR" ]]; then
  echo "ERROR: Dryad data directory not found: $DRYAD_DIR" >&2
  echo "Set GREGARIA_DIET_DRYAD_DIR to the downloaded package location." >&2
  exit 1
fi

for rel_path in "${LINK_PATHS[@]}"; do
  source_path="$DRYAD_DIR/$rel_path"
  target_path="$PROJECT_DIR/$rel_path"

  if [[ ! -e "$source_path" ]]; then
    echo "ERROR: expected Dryad path is missing: $source_path" >&2
    exit 1
  fi

  if [[ -L "$target_path" ]]; then
    current_target=$(readlink "$target_path")
    if [[ "$current_target" == "$source_path" ]]; then
      printf 'Ready: %s -> %s\n' "$rel_path" "$source_path"
      continue
    fi
    echo "ERROR: $target_path links to a different location: $current_target" >&2
    exit 1
  fi

  if [[ -e "$target_path" ]]; then
    echo "ERROR: refusing to replace existing path: $target_path" >&2
    echo "Move or reconcile it before rerunning this setup." >&2
    exit 1
  fi

  mkdir -p "$(dirname "$target_path")"
  ln -s "$source_path" "$target_path"
  printf 'Linked: %s -> %s\n' "$rel_path" "$source_path"
done

echo "Dryad data links are ready. R Markdown project paths remain unchanged."
