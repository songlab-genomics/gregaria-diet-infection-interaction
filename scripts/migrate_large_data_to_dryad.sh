#!/usr/bin/env bash
set -euo pipefail

# One-time migration of large inputs/results out of GitHub's working tree.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
DEFAULT_DRYAD_DIR="/Users/maevatecher/Dropbox/3. Research Papers/Locusts/Ongoing/gregaria-diet-infection-interaction-dryad"
DRYAD_DIR="${GREGARIA_DIET_DRYAD_DIR:-$DEFAULT_DRYAD_DIR}"

MOVE_PATHS=(
  "data/GO_Annotations"
  "data/raw_read_counts"
  "data/reference"
  "data/external"
  "data/scaffold_origin"
  "output/rmd_runs"
  "output/runs"
  "outputs"
)

mkdir -p "$DRYAD_DIR"

for rel_path in "${MOVE_PATHS[@]}"; do
  source_path="$PROJECT_DIR/$rel_path"
  target_path="$DRYAD_DIR/$rel_path"

  if [[ -L "$source_path" ]]; then
    printf 'Already linked: %s\n' "$rel_path"
    continue
  fi
  if [[ ! -e "$source_path" ]]; then
    echo "ERROR: migration source is missing: $source_path" >&2
    exit 1
  fi
  if [[ -e "$target_path" ]]; then
    echo "ERROR: refusing to overwrite Dryad path: $target_path" >&2
    exit 1
  fi

  mkdir -p "$(dirname "$target_path")"
  mv "$source_path" "$target_path"
  printf 'Moved: %s\n' "$rel_path"
done

cp "$PROJECT_DIR/DRYAD_PACKAGE.md" "$DRYAD_DIR/README.md"
echo "Large-data migration complete: $DRYAD_DIR"
