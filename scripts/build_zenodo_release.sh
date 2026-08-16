#!/usr/bin/env bash
set -euo pipefail

# Build a portable code-and-website archive from committed files only.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
PROJECT_NAME=$(basename "$PROJECT_DIR")

cd "$PROJECT_DIR"

if [[ -n "$(git status --porcelain)" ]]; then
  echo "ERROR: the Git working tree is not clean." >&2
  echo "Commit or intentionally exclude pending changes before archiving." >&2
  git status --short >&2
  exit 1
fi

release_label="${1:-$(git describe --tags --exact-match 2>/dev/null || git rev-parse --short HEAD)}"
release_label=$(printf '%s' "$release_label" | tr -cs 'A-Za-z0-9._-' '-')
output_dir="${2:-$PROJECT_DIR/../zenodo-release}"
archive_root="${PROJECT_NAME}-${release_label}"
archive_path="$output_dir/${archive_root}.zip"

mkdir -p "$output_dir"
git archive \
  --format=zip \
  --prefix="${archive_root}/" \
  --output="$archive_path" \
  HEAD

# Confirm that every archived member can be read before reporting success.
unzip -tq "$archive_path"

echo "Zenodo code archive: $archive_path"
echo "Archive size: $(du -h "$archive_path" | awk '{print $1}')"
echo "SHA-256: $(shasum -a 256 "$archive_path" | awk '{print $1}')"
echo "Dryad data are intentionally excluded; deposit that package separately."
