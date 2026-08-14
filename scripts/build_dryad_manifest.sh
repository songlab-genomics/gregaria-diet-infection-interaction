#!/usr/bin/env bash
set -euo pipefail

# Add checksums for new files; set DRYAD_MANIFEST_FULL=true for a full refresh.
DEFAULT_DRYAD_DIR="/Users/maevatecher/Dropbox/3. Research Papers/Locusts/Ongoing/gregaria-diet-infection-interaction-dryad"
DRYAD_DIR="${GREGARIA_DIET_DRYAD_DIR:-$DEFAULT_DRYAD_DIR}"
MANIFEST="$DRYAD_DIR/dryad_file_manifest.tsv"
TMP_MANIFEST="$DRYAD_DIR/.dryad_file_manifest.tsv.tmp"
JOBS="${DRYAD_MANIFEST_JOBS:-4}"
FULL_REFRESH="${DRYAD_MANIFEST_FULL:-false}"

if [[ ! -d "$DRYAD_DIR" ]]; then
  echo "ERROR: Dryad data directory not found: $DRYAD_DIR" >&2
  exit 1
fi

if [[ "$FULL_REFRESH" == "true" ]]; then
  printf 'relative_path\tbytes\tsha256\n' > "$TMP_MANIFEST"
elif [[ ! -f "$TMP_MANIFEST" ]]; then
  if [[ -f "$MANIFEST" ]]; then
    cp "$MANIFEST" "$TMP_MANIFEST"
  else
    printf 'relative_path\tbytes\tsha256\n' > "$TMP_MANIFEST"
  fi
fi

all_paths=$(mktemp)
done_paths=$(mktemp)
remaining_paths=$(mktemp)
new_rows=$(mktemp)
trap 'rm -f "$all_paths" "$done_paths" "$remaining_paths" "$new_rows"' EXIT

find "$DRYAD_DIR" -type f \
  ! -name 'dryad_file_manifest.tsv' \
  ! -name '.dryad_file_manifest.tsv.tmp' \
  -print | sed "s#^$DRYAD_DIR/##" | LC_ALL=C sort > "$all_paths"
tail -n +2 "$TMP_MANIFEST" | cut -f1 | LC_ALL=C sort -u > "$done_paths"
comm -23 "$all_paths" "$done_paths" > "$remaining_paths"

export DRYAD_DIR
xargs -P "$JOBS" -I '{}' bash -c '
  rel_path=$1
  file_path="$DRYAD_DIR/$rel_path"
  bytes=$(stat -f "%z" "$file_path")
  checksum=$(shasum -a 256 "$file_path" | awk "{print \$1}")
  printf "%s\t%s\t%s\n" "$rel_path" "$bytes" "$checksum"
' _ '{}' < "$remaining_paths" > "$new_rows"

{
  printf 'relative_path\tbytes\tsha256\n'
  {
    tail -n +2 "$TMP_MANIFEST"
    cat "$new_rows"
  } | LC_ALL=C sort -t $'\t' -k1,1
} > "$MANIFEST"

rm -f "$TMP_MANIFEST"
echo "Wrote: $MANIFEST"
