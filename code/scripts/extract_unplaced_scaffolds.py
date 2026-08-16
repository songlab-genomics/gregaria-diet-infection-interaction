#!/usr/bin/env python3
"""Extract selected unplaced scaffolds and report sequence properties."""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument(
        "--sequence-list",
        default="",
        help="Optional newline-delimited scaffold IDs; otherwise extract every NW_ sequence.",
    )
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-table", required=True)
    return parser.parse_args()


def fasta_records(path):
    name = None
    description = ""
    chunks = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, description, "".join(chunks)
                header = line[1:].rstrip("\n")
                name = header.split(maxsplit=1)[0]
                description = header
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, description, "".join(chunks)


def main():
    args = parse_args()
    selected = None
    if args.sequence_list:
        selected = {
            line.strip()
            for line in open(args.sequence_list)
            if line.strip() and not line.startswith("#")
        }
    Path(args.output_fasta).parent.mkdir(parents=True, exist_ok=True)
    rows = []
    with open(args.output_fasta, "w") as fasta_out:
        for name, description, sequence in fasta_records(args.fasta):
            if selected is None and not name.startswith("NW_"):
                continue
            if selected is not None and name not in selected:
                continue
            clean = sequence.upper()
            acgt = sum(clean.count(base) for base in "ACGT")
            gc = clean.count("G") + clean.count("C")
            fasta_out.write(f">{description}\n")
            for start in range(0, len(clean), 80):
                fasta_out.write(clean[start:start + 80] + "\n")
            rows.append(
                {
                    "scaffold": name,
                    "length_bp": len(clean),
                    "gc_percent": 100.0 * gc / acgt if acgt else 0.0,
                    "n_percent": 100.0 * clean.count("N") / len(clean) if clean else 0.0,
                }
            )

    if selected is not None:
        found = {row["scaffold"] for row in rows}
        missing = sorted(selected - found)
        if missing:
            raise ValueError(
                "Scaffolds listed but absent from FASTA: " + ", ".join(missing)
            )

    with open(args.output_table, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["scaffold", "length_bp", "gc_percent", "n_percent"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
