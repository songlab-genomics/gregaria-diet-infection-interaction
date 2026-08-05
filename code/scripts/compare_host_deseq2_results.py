#!/usr/bin/env python3
"""Reconcile host DESeq2 results from host-only and competitive mappings."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host-only", required=True)
    parser.add_argument("--competitive", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    return parser.parse_args()


def read_results(path):
    rows = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            rows[(row["contrast_id"], row["gene_id"])] = row
    return rows


def as_bool(value):
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def direction(value):
    if math.isnan(value) or value == 0:
        return "zero_or_missing"
    return "positive" if value > 0 else "negative"


def main():
    args = parse_args()
    host_only = read_results(args.host_only)
    competitive = read_results(args.competitive)
    keys = sorted(set(host_only) | set(competitive))

    rows = []
    summaries = defaultdict(
        lambda: {
            "host_only_significant": 0,
            "competitive_significant": 0,
            "shared_same_direction": 0,
            "shared_opposite_direction": 0,
            "host_only_only": 0,
            "competitive_only": 0,
        }
    )
    for contrast_id in sorted({key[0] for key in keys}):
        summaries[contrast_id]

    for contrast_id, gene_id in keys:
        host = host_only.get((contrast_id, gene_id), {})
        comp = competitive.get((contrast_id, gene_id), {})
        host_sig = as_bool(host.get("significant", False))
        comp_sig = as_bool(comp.get("significant", False))
        if not host_sig and not comp_sig:
            continue

        host_lfc = as_float(host.get("log2FoldChange"))
        comp_lfc = as_float(comp.get("log2FoldChange"))
        if host_sig and comp_sig:
            if direction(host_lfc) == direction(comp_lfc):
                classification = "shared_same_direction"
            else:
                classification = "shared_opposite_direction"
        elif host_sig:
            classification = "host_only_only"
        else:
            classification = "competitive_only"

        summary = summaries[contrast_id]
        summary["host_only_significant"] += int(host_sig)
        summary["competitive_significant"] += int(comp_sig)
        summary[classification] += 1
        rows.append(
            {
                "contrast_id": contrast_id,
                "contrast_family": host.get(
                    "contrast_family", comp.get("contrast_family", "")
                ),
                "gene_id": gene_id,
                "host_only_log2FoldChange": host.get("log2FoldChange", ""),
                "competitive_log2FoldChange": comp.get("log2FoldChange", ""),
                "host_only_padj": host.get("padj", ""),
                "competitive_padj": comp.get("padj", ""),
                "host_only_significant": host_sig,
                "competitive_significant": comp_sig,
                "mapping_support_class": classification,
            }
        )

    output_fields = [
        "contrast_id",
        "contrast_family",
        "gene_id",
        "host_only_log2FoldChange",
        "competitive_log2FoldChange",
        "host_only_padj",
        "competitive_padj",
        "host_only_significant",
        "competitive_significant",
        "mapping_support_class",
    ]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields)
        writer.writeheader()
        writer.writerows(rows)

    summary_fields = [
        "contrast_id",
        "host_only_significant",
        "competitive_significant",
        "shared_same_direction",
        "shared_opposite_direction",
        "host_only_only",
        "competitive_only",
        "significant_union",
        "significant_jaccard",
    ]
    with open(args.summary, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=summary_fields, delimiter="\t"
        )
        writer.writeheader()
        for contrast_id in sorted(summaries):
            summary = summaries[contrast_id]
            shared = (
                summary["shared_same_direction"]
                + summary["shared_opposite_direction"]
            )
            union = shared + summary["host_only_only"] + summary["competitive_only"]
            writer.writerow(
                {
                    "contrast_id": contrast_id,
                    **summary,
                    "significant_union": union,
                    "significant_jaccard": shared / union if union else "",
                }
            )


if __name__ == "__main__":
    main()
