#!/usr/bin/env python3
"""Convert one STAR Log.final.out file into a single-row QC table."""

import argparse
import csv
from pathlib import Path


FIELDS = {
    "Number of input reads": "input_reads",
    "Uniquely mapped reads number": "uniquely_mapped_reads",
    "Uniquely mapped reads %": "uniquely_mapped_percent",
    "Number of reads mapped to multiple loci": "multimapped_reads",
    "% of reads mapped to multiple loci": "multimapped_percent",
    "Number of reads mapped to too many loci": "too_many_loci_reads",
    "% of reads mapped to too many loci": "too_many_loci_percent",
    "% of reads unmapped: too many mismatches": "unmapped_mismatches_percent",
    "% of reads unmapped: too short": "unmapped_too_short_percent",
    "% of reads unmapped: other": "unmapped_other_percent",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--branch", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def clean_value(value):
    value = value.strip()
    if value.endswith("%"):
        return value[:-1]
    return value


def main():
    args = parse_args()
    values = {}
    with open(args.input) as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = (part.strip() for part in line.split("|", 1))
            if key in FIELDS:
                values[FIELDS[key]] = clean_value(value)

    missing = [output_name for output_name in FIELDS.values() if output_name not in values]
    if missing:
        raise ValueError(
            f"STAR log {args.input} is missing expected fields: {', '.join(missing)}"
        )

    row = {"sample_id": args.sample, "mapping_branch": args.branch, **values}
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
