#!/usr/bin/env python3
"""Combine one-row TSV files while preserving their shared header."""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    rows = []
    fieldnames = []
    for path in args.input:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for name in reader.fieldnames or []:
                if name not in fieldnames:
                    fieldnames.append(name)
            rows.extend(reader)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
