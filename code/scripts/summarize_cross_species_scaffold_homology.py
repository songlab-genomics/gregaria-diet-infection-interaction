#!/usr/bin/env python3
"""Summarize minimap2 PAF coverage of candidate scaffolds in related genomes."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-list", required=True)
    parser.add_argument(
        "--alignment",
        action="append",
        default=[],
        help="Repeated species=/path/to/alignment.paf values.",
    )
    parser.add_argument(
        "--minimum-mapq",
        type=int,
        default=20,
        help="Minimum MAPQ for confident cross-species coverage (default: 20).",
    )
    parser.add_argument(
        "--minimum-identity",
        type=float,
        default=70,
        help="Minimum minimap2 divergence-derived identity (default: 70).",
    )
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def merge_intervals(intervals):
    if not intervals:
        return 0
    merged_bp = 0
    current_start, current_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged_bp += current_end - current_start
            current_start, current_end = start, end
    return merged_bp + current_end - current_start


def summarize_paf(path):
    records = defaultdict(list)
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            tags = {}
            for value in fields[12:]:
                parts = value.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            try:
                divergence_identity = 100 * (1 - float(tags["dv"]))
            except (KeyError, ValueError):
                divergence_identity = (
                    100
                    * int(fields[9])
                    / max(int(fields[10]), 1)
                )
            records[fields[0]].append(
                {
                    "query_length": int(fields[1]),
                    "query_start": int(fields[2]),
                    "query_end": int(fields[3]),
                    "target": fields[5],
                    "matching_bases": int(fields[9]),
                    "alignment_block_length": int(fields[10]),
                    "mapq": int(fields[11]),
                    "divergence_identity_percent": divergence_identity,
                }
            )
    return records


def main():
    args = parse_args()
    with open(args.candidate_list, encoding="utf-8") as handle:
        candidates = [
            line.strip()
            for line in handle
            if line.strip() and not line.startswith("#")
        ]

    alignments = []
    for value in args.alignment:
        if "=" not in value:
            raise ValueError("--alignment must use species=/path/to/file.paf")
        species, path = value.split("=", 1)
        alignments.append((species, path))

    rows = []
    if not alignments:
        for scaffold in candidates:
            rows.append(
                {
                    "scaffold": scaffold,
                    "comparison_species": "not_run",
                    "query_length_bp": "",
                    "aligned_query_bp": "",
                    "query_coverage_percent": "",
                    "weighted_identity_percent": "",
                    "best_target_sequence": "",
                    "best_alignment_identity_percent": "",
                    "alignment_records": 0,
                    "high_conf_aligned_query_bp": "",
                    "high_conf_query_coverage_percent": "",
                    "high_conf_weighted_identity_percent": "",
                    "high_conf_alignment_records": "",
                }
            )
    else:
        for species, path in alignments:
            records = summarize_paf(path)
            for scaffold in candidates:
                hits = records.get(scaffold, [])
                if not hits:
                    rows.append(
                        {
                            "scaffold": scaffold,
                            "comparison_species": species,
                            "query_length_bp": "",
                            "aligned_query_bp": 0,
                            "query_coverage_percent": 0,
                            "weighted_identity_percent": 0,
                            "best_target_sequence": "",
                            "best_alignment_identity_percent": 0,
                            "alignment_records": 0,
                            "high_conf_aligned_query_bp": 0,
                            "high_conf_query_coverage_percent": 0,
                            "high_conf_weighted_identity_percent": 0,
                            "high_conf_alignment_records": 0,
                        }
                    )
                    continue

                query_length = hits[0]["query_length"]
                aligned_bp = merge_intervals(
                    [(hit["query_start"], hit["query_end"]) for hit in hits]
                )
                total_query_span = sum(
                    hit["query_end"] - hit["query_start"] for hit in hits
                )
                weighted_identity = (
                    sum(
                        (hit["query_end"] - hit["query_start"])
                        * hit["divergence_identity_percent"]
                        for hit in hits
                    )
                    / total_query_span
                    if total_query_span
                    else 0
                )
                high_conf = [
                    hit
                    for hit in hits
                    if hit["mapq"] >= args.minimum_mapq
                    and hit["divergence_identity_percent"]
                    >= args.minimum_identity
                ]
                high_conf_aligned_bp = merge_intervals(
                    [
                        (hit["query_start"], hit["query_end"])
                        for hit in high_conf
                    ]
                )
                high_conf_query_span = sum(
                    hit["query_end"] - hit["query_start"]
                    for hit in high_conf
                )
                high_conf_identity = (
                    sum(
                        (hit["query_end"] - hit["query_start"])
                        * hit["divergence_identity_percent"]
                        for hit in high_conf
                    )
                    / high_conf_query_span
                    if high_conf_query_span
                    else 0
                )
                best = max(
                    hits,
                    key=lambda hit: (
                        hit["mapq"],
                        hit["divergence_identity_percent"],
                        hit["query_end"] - hit["query_start"],
                    ),
                )
                rows.append(
                    {
                        "scaffold": scaffold,
                        "comparison_species": species,
                        "query_length_bp": query_length,
                        "aligned_query_bp": aligned_bp,
                        "query_coverage_percent": 100 * aligned_bp / query_length,
                        "weighted_identity_percent": weighted_identity,
                        "best_target_sequence": best["target"],
                        "best_alignment_identity_percent": best[
                            "divergence_identity_percent"
                        ],
                        "alignment_records": len(hits),
                        "high_conf_aligned_query_bp": high_conf_aligned_bp,
                        "high_conf_query_coverage_percent": (
                            100 * high_conf_aligned_bp / query_length
                        ),
                        "high_conf_weighted_identity_percent": high_conf_identity,
                        "high_conf_alignment_records": len(high_conf),
                    }
                )

    fields = [
        "scaffold",
        "comparison_species",
        "query_length_bp",
        "aligned_query_bp",
        "query_coverage_percent",
        "weighted_identity_percent",
        "best_target_sequence",
        "best_alignment_identity_percent",
        "alignment_records",
        "high_conf_aligned_query_bp",
        "high_conf_query_coverage_percent",
        "high_conf_weighted_identity_percent",
        "high_conf_alignment_records",
    ]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
