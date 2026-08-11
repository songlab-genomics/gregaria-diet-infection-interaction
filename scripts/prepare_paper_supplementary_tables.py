#!/usr/bin/env python3
"""Prepare auditable JSON tables for the paper supplementary workbook."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def latest_run(root: Path, prefix: str, required_file: str) -> Path:
    candidates = sorted(
        path
        for path in root.glob(f"{prefix}_*")
        if path.is_dir() and (path / required_file).is_file()
    )
    if not candidates:
        raise FileNotFoundError(
            f"No completed {prefix} run contains {required_file} under {root}"
        )
    return candidates[-1]


def clean_for_json(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame.columns = [str(column) for column in frame.columns]
    # Object dtype allows missing numeric values to become JSON null, not NaN.
    return frame.astype(object).where(pd.notna(frame), None)


def frame_payload(frame: pd.DataFrame) -> dict:
    clean = clean_for_json(frame)
    return {
        "columns": clean.columns.tolist(),
        "rows": clean.values.tolist(),
        "row_count": int(len(clean)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    project = args.project.resolve()
    output_path = args.output.resolve()
    run_root = project / "output" / "rmd_runs"

    catalogue_run = latest_run(
        run_root,
        "main-chromosome-all-deg-catalogue",
        "all_samples_deg_catalogue.csv",
    )
    enrichment_run = latest_run(
        run_root, "go-term-enrichment", "GO_enrichment_all_sets.csv"
    )
    publication_run = latest_run(
        run_root,
        "publication-figures",
        "figureS3_candidate_evidence_master.csv",
    )
    qc_run = latest_run(
        run_root,
        "count-qc-exploration",
        "mapping_quality_and_individual_mass_45_samples.csv",
    )
    burden_run = latest_run(
        run_root,
        "metarhizium-read-burden",
        "bracken_host_unmapped_burden_by_sample.csv",
    )

    deg = pd.read_csv(catalogue_run / "all_samples_deg_catalogue.csv")
    deg.insert(0, "count_definition", "Transcript + exon")
    deg.insert(1, "analysis_scope", "Nuclear chromosomes only")
    deg.insert(2, "analysis_cohort", "44 QC-passed libraries; 1044 excluded")

    rrna_file = project / "data" / "excluded_loci" / "gregaria_rrna_list.txt"
    rrna_ids = {
        line.strip()
        for line in rrna_file.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }
    rrna_overlap = sorted(set(deg["gene_id"].astype(str)) & rrna_ids)
    if rrna_overlap:
        raise ValueError(
            f"The primary DEG catalogue unexpectedly contains {len(rrna_overlap)} rRNA IDs"
        )

    go = pd.read_csv(enrichment_run / "GO_enrichment_all_sets.csv")
    go.insert(0, "count_definition", "Transcript + exon")
    go.insert(1, "analysis_scope", "Nuclear chromosomes only")

    kegg = pd.read_csv(enrichment_run / "KEGG_pathway_enrichment_all_sets.csv")
    kegg.insert(0, "count_definition", "Transcript + exon")
    kegg.insert(1, "analysis_scope", "Nuclear chromosomes only")

    candidates = pd.read_csv(
        publication_run / "figureS3_candidate_evidence_master.csv"
    )
    candidate_results = pd.read_csv(
        publication_run / "figureS3_candidate_contrast_results.csv"
    )

    manual_n = int(
        (candidates["evidence_tier"] == "Manual curation candidate").sum()
    )
    eggnog_supported_n = int(
        (candidates["evidence_tier"] == "eggNOG-supported candidate").sum()
    )
    if manual_n != 115:
        raise ValueError(f"Expected 115 manual curation candidates, found {manual_n}")

    mapping = pd.read_csv(
        qc_run / "mapping_quality_and_individual_mass_45_samples.csv",
        dtype={"label": str},
    )
    competitive = pd.read_csv(
        project
        / "output"
        / "runs"
        / "host_pathogen_dual_corrected45_20260809-153146"
        / "05-mapping-comparison"
        / "competitive_mapping_summary.tsv",
        sep="\t",
        dtype={"sample_id": str},
    )
    microbial = pd.read_csv(
        burden_run / "bracken_host_unmapped_burden_by_sample.csv",
        dtype={"label": str},
    )

    sample = (
        mapping.merge(
            competitive[
                [
                    "sample_id",
                    "primary_mapped_read_alignments",
                    "unique_read_alignments",
                    "multimapping_primary_read_alignments",
                    "unique_host_chromosome",
                    "unique_host_unplaced",
                    "unique_fungus",
                    "pct_unique_host_chromosome",
                    "pct_unique_host_unplaced",
                    "pct_unique_fungus",
                ]
            ],
            left_on="label",
            right_on="sample_id",
            how="left",
            validate="one_to_one",
        )
        .merge(
            microbial[
                [
                    "label",
                    "microbial_family_reads",
                    "microbial_reads_per_million",
                    "log10_microbial_rpm",
                ]
            ],
            on="label",
            how="left",
            validate="one_to_one",
        )
        .drop(columns=["sample_id"])
    )
    sample.insert(1, "analysis_included", sample["label"] != "1044")
    sample.insert(
        2,
        "analysis_status",
        sample["analysis_included"].map(
            {
                True: "Included in primary 44-library analysis",
                False: "QC excluded: exceptionally low host alignment",
            }
        ),
    )
    sample = sample.sort_values(["Treatment", "Diet", "label"], kind="stable")

    manifest = pd.DataFrame(
        [
            {
                "table": "Sample metadata and mapping",
                "rows": len(sample),
                "source": str(qc_run),
                "notes": "All 45 sequenced libraries; 1044 is retained only for QC documentation.",
            },
            {
                "table": "Transcript-plus-exon DEG catalogue",
                "rows": len(deg),
                "source": str(catalogue_run),
                "notes": "Significant non-rRNA DEGs on nuclear chromosomes; 44 QC-passed libraries.",
            },
            {
                "table": "GO enrichment",
                "rows": len(go),
                "source": str(enrichment_run),
                "notes": "eggNOG-derived GO over-representation results with contributing gene IDs.",
            },
            {
                "table": "KEGG enrichment",
                "rows": len(kegg),
                "source": str(enrichment_run),
                "notes": "eggNOG-derived KEGG over-representation results with pathway scope.",
            },
            {
                "table": "Physiological candidate master",
                "rows": len(candidates),
                "source": str(publication_run),
                "notes": (
                    f"{manual_n} manual curation candidates plus "
                    f"{eggnog_supported_n} eggNOG-supported candidate DEGs."
                ),
            },
            {
                "table": "Candidate contrast results",
                "rows": len(candidate_results),
                "source": str(publication_run),
                "notes": "All displayed contrast results, including non-significant cells, for selected candidates.",
            },
        ]
    )

    payload = {
        "metadata": {
            "title": "Paper supplementary data tables",
            "analysis": "Primary transcript-plus-exon analysis",
            "cohort": "44 QC-passed libraries; sample 1044 shown only in the QC sheet",
            "gene_scope": "Non-rRNA genes on nuclear chromosome pseudomolecules",
            "candidate_rule": (
                "The 115 manual curation candidates form the primary set; "
                "eggNOG-supported candidates must be absent from that set and "
                "significant in at least one displayed contrast"
            ),
            "rrna_overlap_in_deg_catalogue": len(rrna_overlap),
        },
        "tables": {
            "Sources": frame_payload(manifest),
            "Samples_mapping": frame_payload(sample),
            "DEG_transcript_exon": frame_payload(deg),
            "GO_enrichment": frame_payload(go),
            "KEGG_enrichment": frame_payload(kegg),
            "Physiology_candidates": frame_payload(candidates),
            "Candidate_DEG_results": frame_payload(candidate_results),
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, ensure_ascii=True))


if __name__ == "__main__":
    main()
