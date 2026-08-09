#!/usr/bin/env python3
"""Combine the frozen 44-sample matrices with the independently run sample 1044.

The script never edits source runs. It validates gene order, count integrity,
sample identities, corrected treatment assignments, and rRNA-list provenance
before writing a new run-scoped 45-sample dataset.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from collections import Counter
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path


MATRIX_SPECS = (
    (
        "competitive-host/host_transcript_exon_counts.tsv",
        "competitive-host/host_transcript_exon_sensitivity_counts.tsv",
        "competitive-host/host_transcript_exon_counts.tsv",
        "Transcript + exon",
        "Competitive host",
    ),
    (
        "host-only/host_transcript_exon_counts.tsv",
        "host-only/host_transcript_exon_sensitivity_counts.tsv",
        "host-only/host_transcript_exon_counts.tsv",
        "Transcript + exon",
        "Host only",
    ),
    (
        "competitive-host/host_exon_counts.tsv",
        "competitive-host/host_exon_counts.tsv",
        "competitive-host/host_exon_counts.tsv",
        "Exon only",
        "Competitive host",
    ),
    (
        "host-only/host_exon_counts.tsv",
        "host-only/host_exon_counts.tsv",
        "host-only/host_exon_counts.tsv",
        "Exon only",
        "Host only",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old-run", type=Path, required=True)
    parser.add_argument("--addon-run", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--rrna-list", type=Path, required=True)
    parser.add_argument("--output-run", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_metadata(path: Path) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"label", "Treatment", "Diet"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Metadata lacks required columns: {path}")
    by_sample = {row["label"]: row for row in rows}
    if len(rows) != 45 or len(by_sample) != 45:
        raise ValueError("Corrected metadata must contain 45 unique samples.")
    if (by_sample["1024"]["Treatment"], by_sample["1024"]["Diet"]) != ("Infected", "33"):
        raise ValueError("Sample 1024 must be assigned to infected diet 33.")
    if (by_sample["1044"]["Treatment"], by_sample["1044"]["Diet"]) != ("Infected", "83"):
        raise ValueError("Sample 1044 must be assigned to infected diet 83.")
    observed = Counter((row["Treatment"], row["Diet"]) for row in rows)
    expected = Counter(
        {
            ("Control", "33"): 7,
            ("Control", "50"): 6,
            ("Control", "83"): 7,
            ("Infected", "33"): 10,
            ("Infected", "50"): 9,
            ("Infected", "83"): 6,
        }
    )
    if observed != expected:
        raise ValueError(f"Unexpected corrected group sizes: {dict(observed)}")
    return rows, by_sample


def read_header(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def combine_matrix(old_path: Path, addon_path: Path, output_path: Path) -> dict[str, str | int]:
    old_header = read_header(old_path)
    addon_header = read_header(addon_path)
    if old_header[0] != "gene_id" or addon_header != ["gene_id", "1044"]:
        raise ValueError(f"Unexpected matrix headers: {old_path} / {addon_path}")
    old_samples = old_header[1:]
    if len(old_samples) != 44 or "1044" in old_samples:
        raise ValueError(f"Expected a frozen 44-sample matrix: {old_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    genes = 0
    seen: set[str] = set()
    with (
        old_path.open(newline="") as old_handle,
        addon_path.open(newline="") as addon_handle,
        output_path.open("w", newline="") as output_handle,
    ):
        old_reader = csv.reader(old_handle, delimiter="\t")
        addon_reader = csv.reader(addon_handle, delimiter="\t")
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        next(old_reader)
        next(addon_reader)
        writer.writerow(["gene_id", *old_samples, "1044"])
        for old_row, addon_row in zip_longest(old_reader, addon_reader):
            if old_row is None or addon_row is None:
                raise ValueError(f"Matrix row counts differ: {old_path} / {addon_path}")
            gene_id = old_row[0]
            if addon_row[0] != gene_id:
                raise ValueError(f"Gene order mismatch at {gene_id}: {addon_path}")
            if gene_id in seen:
                raise ValueError(f"Duplicate gene ID in matrix: {gene_id}")
            seen.add(gene_id)
            values = old_row[1:] + addon_row[1:]
            if any(not value.isdigit() for value in values):
                raise ValueError(f"Non-integer or negative count for {gene_id}")
            writer.writerow([gene_id, *values])
            genes += 1

    return {
        "genes": genes,
        "old_samples": len(old_samples),
        "addon_samples": 1,
        "combined_samples": len(old_samples) + 1,
        "sha256": sha256(output_path),
    }


def combine_row_tables(old_path: Path, addon_path: Path, output_path: Path) -> int:
    with old_path.open(newline="") as old_handle, addon_path.open(newline="") as addon_handle:
        old_reader = csv.DictReader(old_handle, delimiter="\t")
        addon_reader = csv.DictReader(addon_handle, delimiter="\t")
        if old_reader.fieldnames != addon_reader.fieldnames:
            raise ValueError(f"Table schema mismatch: {old_path} / {addon_path}")
        rows = [*old_reader, *addon_reader]
        fieldnames = list(old_reader.fieldnames or [])
    if "sample_id" in fieldnames:
        rows.sort(key=lambda row: int(row["sample_id"]))
        if len({row["sample_id"] for row in rows}) != 45:
            raise ValueError(f"Combined table does not contain 45 unique samples: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def combine_taxonomy_manifest(
    old_run: Path,
    addon_run: Path,
    metadata: dict[str, dict[str, str]],
    output_path: Path,
) -> int:
    rows: list[dict[str, str]] = []
    for source_run in (old_run, addon_run):
        manifest = source_run / "08-taxonomy/family_abundance_manifest.tsv"
        with manifest.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample = row["sample_id"]
                source = row["source"]
                row["treatment"] = metadata[sample]["Treatment"]
                row["diet"] = metadata[sample]["Diet"]
                row["abundance_file"] = str(
                    source_run / f"08-taxonomy/{source}/{sample}.bracken.family.tsv"
                )
                row["kraken_report"] = str(
                    source_run / f"08-taxonomy/{source}/{sample}.kraken2.report"
                )
                if not Path(row["abundance_file"]).is_file() or not Path(row["kraken_report"]).is_file():
                    raise FileNotFoundError(f"Missing local taxonomy output for {sample} / {source}")
                rows.append(row)
    rows.sort(key=lambda row: (int(row["sample_id"]), row["source"]))
    if len(rows) != 45 * 3:
        raise ValueError(f"Expected 135 taxonomy-manifest rows, found {len(rows)}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main() -> None:
    args = parse_args()
    if args.output_run.exists():
        raise FileExistsError(f"Refusing to overwrite existing run: {args.output_run}")

    metadata_rows, metadata_by_sample = read_metadata(args.metadata)
    old_rrna = args.old_run / "inputs/gregaria_rrna_list.txt"
    addon_rrna = args.addon_run / "inputs/gregaria_rrna_list.txt"
    rrna_hashes = {sha256(args.rrna_list), sha256(old_rrna), sha256(addon_rrna)}
    if len(rrna_hashes) != 1:
        raise ValueError("The project, frozen-run, and add-on rRNA lists differ.")

    matrix_audit: list[dict[str, str | int]] = []
    for output_rel, old_rel, addon_rel, definition, mapping_branch in MATRIX_SPECS:
        old_path = args.old_run / "04-count-matrices" / old_rel
        addon_path = args.addon_run / "04-count-matrices" / addon_rel
        output_path = args.output_run / "04-count-matrices" / output_rel
        audit = combine_matrix(old_path, addon_path, output_path)
        matrix_audit.append(
            {
                "mapping_branch": mapping_branch,
                "count_definition": definition,
                "source_44": str(old_path),
                "source_1044": str(addon_path),
                "output_45": str(output_path),
                **audit,
            }
        )

    for filename in ("host_only_mapping_summary.tsv", "competitive_mapping_summary.tsv"):
        combine_row_tables(
            args.old_run / "05-mapping-comparison" / filename,
            args.addon_run / "05-mapping-comparison" / filename,
            args.output_run / "05-mapping-comparison" / filename,
        )

    taxonomy_rows = combine_taxonomy_manifest(
        args.old_run,
        args.addon_run,
        metadata_by_sample,
        args.output_run / "08-taxonomy/family_abundance_manifest.tsv",
    )

    inputs = args.output_run / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.metadata, inputs / "mehreen_metadata_corrected_45.txt")
    shutil.copy2(args.rrna_list, inputs / "gregaria_rrna_list.txt")

    with (args.output_run / "matrix_integration_audit.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(matrix_audit[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(matrix_audit)

    with (args.output_run / "sample_assignment_audit.tsv").open("w", newline="") as handle:
        fieldnames = ["sample_id", "treatment", "diet", "assignment_note"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(metadata_rows, key=lambda item: int(item["label"])):
            note = ""
            if row["label"] == "1024":
                note = "Corrected from infected diet 83 to infected diet 33"
            elif row["label"] == "1044":
                note = "Missing library added as infected diet 83"
            writer.writerow(
                {
                    "sample_id": row["label"],
                    "treatment": row["Treatment"],
                    "diet": row["Diet"],
                    "assignment_note": note,
                }
            )

    provenance = [
        ("created_utc", datetime.now(timezone.utc).isoformat()),
        ("source_44_sample_run", str(args.old_run)),
        ("source_1044_addon_run", str(args.addon_run)),
        ("corrected_metadata", str(args.metadata)),
        ("samples", "45"),
        ("sample_1024_assignment", "Infected diet 33"),
        ("sample_1044_assignment", "Infected diet 83"),
        ("primary_count_definition", "featureCounts -t transcript,exon -g gene_id"),
        ("sensitivity_count_definition", "featureCounts -t exon -g gene_id"),
        ("rrna_list_sha256", next(iter(rrna_hashes))),
        ("taxonomy_manifest_rows", str(taxonomy_rows)),
        ("fungal_count_note", "Not merged: old fungal matrix is exon-only; add-on fungal matrix is transcript+exon"),
    ]
    with (args.output_run / "run_provenance.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["parameter", "value"])
        writer.writerows(provenance)

    print(f"Corrected 45-sample run created: {args.output_run}")
    print("Validated assignments: 1024 = infected diet 33; 1044 = infected diet 83")
    print("Host matrices: transcript+exon and exon-only, host-only and competitive")


if __name__ == "__main__":
    main()
