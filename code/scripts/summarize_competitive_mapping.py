#!/usr/bin/env python3
"""Summarize primary STAR alignments from a host-fungus competitive BAM.

The script reads SAM records from stdin. Feed it primary mapped records with:
samtools view -F 2308 combined.bam | summarize_competitive_mapping.py ...
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--host-prefix", default="HOST__")
    parser.add_argument("--fungus-prefix", default="MR__")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def alignment_category(reference, host_prefix, fungus_prefix):
    if reference.startswith(fungus_prefix):
        return "fungus"
    if reference.startswith(host_prefix + "NC_"):
        return "host_chromosome"
    if reference.startswith(host_prefix + "NW_"):
        return "host_unplaced"
    if reference.startswith(host_prefix):
        return "host_other"
    return "other_reference"


def nh_value(fields):
    for field in fields[11:]:
        if field.startswith("NH:i:"):
            return int(field.rsplit(":", 1)[1])
    return 1


def main():
    args = parse_args()
    counts = Counter()
    for line in sys.stdin:
        if not line.strip() or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        category = alignment_category(
            fields[2],
            args.host_prefix,
            args.fungus_prefix,
        )
        counts["primary_mapped_read_alignments"] += 1
        if nh_value(fields) == 1:
            counts["unique_read_alignments"] += 1
            counts[f"unique_{category}"] += 1
        else:
            counts["multimapping_primary_read_alignments"] += 1
            counts[f"multimapping_{category}"] += 1

    total_unique = counts["unique_read_alignments"]
    row = {"sample_id": args.sample}
    for key in [
        "primary_mapped_read_alignments",
        "unique_read_alignments",
        "multimapping_primary_read_alignments",
        "unique_host_chromosome",
        "unique_host_unplaced",
        "unique_host_other",
        "unique_fungus",
        "unique_other_reference",
    ]:
        row[key] = counts[key]
    for category in ["host_chromosome", "host_unplaced", "host_other", "fungus"]:
        numerator = counts[f"unique_{category}"]
        row[f"pct_unique_{category}"] = (
            100.0 * numerator / total_unique if total_unique else 0.0
        )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
