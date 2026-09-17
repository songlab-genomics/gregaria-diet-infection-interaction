#!/usr/bin/env bash
set -euo pipefail

# Expand the versioned NCBI annotation without modifying the committed archive.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
REFERENCE_DIR="$PROJECT_DIR/data/reference"
GFF="$REFERENCE_DIR/GCF_023897955.1_iqSchGreg1.2_genomic.gff"
GFF_GZ="${GFF}.gz"

if [[ ! -s "$GFF_GZ" ]]; then
  echo "ERROR: missing compressed reference annotation: $GFF_GZ" >&2
  exit 1
fi

gzip -t "$GFF_GZ"

if [[ -s "$GFF" ]]; then
  if gzip -dc "$GFF_GZ" | cmp -s - "$GFF"; then
    echo "Reference annotation is already prepared: $GFF"
    exit 0
  fi
  echo "ERROR: existing GFF differs from the committed compressed copy: $GFF" >&2
  exit 1
fi

tmp="${GFF}.tmp.$$"
trap 'rm -f "$tmp"' EXIT
gzip -dc "$GFF_GZ" > "$tmp"
mv "$tmp" "$GFF"
trap - EXIT

echo "Reference annotation prepared: $GFF"
