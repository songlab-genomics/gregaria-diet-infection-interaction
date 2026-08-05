ORIGIN_INPUT_DIR = ORIGIN_DIR + "/00-input"
ORIGIN_SEQUENCE_DIR = ORIGIN_DIR + "/01-sequences"
ORIGIN_ANNOTATION_DIR = ORIGIN_DIR + "/02-annotation"
ORIGIN_TAXONOMY_DIR = ORIGIN_DIR + "/03-taxonomy"
ORIGIN_CROSS_SPECIES_DIR = ORIGIN_DIR + "/04-cross-species"
ORIGIN_EVIDENCE_DIR = ORIGIN_DIR + "/05-evidence"
ORIGIN_TRANSFER_DIR = ORIGIN_DIR + "/06-local-transfer"

ORIGIN_STAGED_LIST = ORIGIN_INPUT_DIR + "/candidate_scaffolds.txt"
ORIGIN_STAGED_GENES = ORIGIN_INPUT_DIR + "/candidate_deg_genes.tsv"
ORIGIN_STAGED_PROVENANCE = ORIGIN_INPUT_DIR + "/candidate_provenance.tsv"
ORIGIN_STAGED_FCS = ORIGIN_INPUT_DIR + "/ncbi_fcs_report.txt"
ORIGIN_INPUT_MANIFEST = ORIGIN_INPUT_DIR + "/input_manifest.tsv"

ORIGIN_FASTA = ORIGIN_SEQUENCE_DIR + "/candidate_scaffolds.fna"
ORIGIN_SEQUENCE_METRICS = ORIGIN_SEQUENCE_DIR + "/candidate_scaffold_sequence_metrics.tsv"
ORIGIN_ALL_GENES = ORIGIN_ANNOTATION_DIR + "/candidate_scaffold_all_annotated_genes.tsv"
ORIGIN_ANNOTATION_SUMMARY = (
    ORIGIN_ANNOTATION_DIR + "/candidate_scaffold_annotation_summary.tsv"
)
ORIGIN_KRAKEN_REPORT = ORIGIN_TAXONOMY_DIR + "/candidate_scaffolds.kraken2.report"
ORIGIN_KRAKEN_OUTPUT = ORIGIN_TAXONOMY_DIR + "/candidate_scaffolds.kraken2.output"
ORIGIN_TAXONOMY = ORIGIN_TAXONOMY_DIR + "/candidate_scaffold_taxonomy.tsv"
ORIGIN_CROSS_SPECIES_SUMMARY = (
    ORIGIN_CROSS_SPECIES_DIR + "/cross_species_scaffold_homology.tsv"
)
ORIGIN_EVIDENCE = ORIGIN_EVIDENCE_DIR + "/scaffold_origin_evidence.tsv"
ORIGIN_FOLLOW_UP = ORIGIN_EVIDENCE_DIR + "/scaffold_origin_follow_up.tsv"
ORIGIN_TRANSFER_MANIFEST = ORIGIN_TRANSFER_DIR + "/transfer_manifest.tsv"
ORIGIN_TRANSFER_BUNDLE = (
    ORIGIN_TRANSFER_DIR + "/scaffold_origin_audit_for_local.tar.gz"
)


def origin_cross_species_pafs():
    if not ORIGIN_RUN_CROSS_SPECIES:
        return []
    return [
        ORIGIN_CROSS_SPECIES_DIR + f"/{species}.paf"
        for species in sorted(ORIGIN_CROSS_SPECIES_REFERENCES)
    ]


def scaffold_origin_transfer_files():
    return [
        ORIGIN_INPUT_MANIFEST,
        ORIGIN_STAGED_LIST,
        ORIGIN_STAGED_GENES,
        ORIGIN_STAGED_PROVENANCE,
        ORIGIN_STAGED_FCS,
        ORIGIN_FASTA,
        ORIGIN_SEQUENCE_METRICS,
        ORIGIN_ALL_GENES,
        ORIGIN_ANNOTATION_SUMMARY,
        ORIGIN_KRAKEN_REPORT,
        ORIGIN_KRAKEN_OUTPUT,
        ORIGIN_TAXONOMY,
        ORIGIN_CROSS_SPECIES_SUMMARY,
        *origin_cross_species_pafs(),
        ORIGIN_EVIDENCE,
        ORIGIN_FOLLOW_UP,
    ]


def scaffold_origin_audit_outputs():
    return [
        ORIGIN_EVIDENCE,
        ORIGIN_FOLLOW_UP,
        ORIGIN_TRANSFER_MANIFEST,
        ORIGIN_TRANSFER_BUNDLE,
    ]


rule scaffold_origin_audit_hpc:
    input:
        scaffold_origin_audit_outputs()


rule origin_stage_inputs:
    input:
        candidates=ORIGIN_CANDIDATE_LIST,
        genes=ORIGIN_CANDIDATE_GENES,
        provenance=ORIGIN_CANDIDATE_PROVENANCE,
        fcs=ORIGIN_NCBI_FCS_REPORT
    output:
        candidates=ORIGIN_STAGED_LIST,
        genes=ORIGIN_STAGED_GENES,
        provenance=ORIGIN_STAGED_PROVENANCE,
        fcs=ORIGIN_STAGED_FCS,
        manifest=ORIGIN_INPUT_MANIFEST
    shell:
        r"""
        python3 scripts/stage_scaffold_origin_inputs.py \
          --candidate-list {input.candidates} \
          --candidate-genes {input.genes} \
          --provenance {input.provenance} \
          --fcs-report {input.fcs} \
          --output-list {output.candidates} \
          --output-genes {output.genes} \
          --output-provenance {output.provenance} \
          --output-fcs-report {output.fcs} \
          --manifest {output.manifest}
        """


rule origin_extract_candidate_scaffolds:
    input:
        fasta=ORIGIN_HOST_FASTA,
        candidates=ORIGIN_STAGED_LIST
    output:
        fasta=ORIGIN_FASTA,
        metrics=ORIGIN_SEQUENCE_METRICS
    shell:
        r"""
        python3 scripts/extract_unplaced_scaffolds.py \
          --fasta {input.fasta} \
          --sequence-list {input.candidates} \
          --output-fasta {output.fasta} \
          --output-table {output.metrics}
        """


rule origin_summarize_candidate_annotation:
    input:
        candidates=ORIGIN_STAGED_LIST,
        candidate_genes=ORIGIN_STAGED_GENES,
        annotation=ORIGIN_HOST_ANNOTATION,
        metrics=ORIGIN_SEQUENCE_METRICS
    output:
        genes=ORIGIN_ALL_GENES,
        scaffolds=ORIGIN_ANNOTATION_SUMMARY
    shell:
        r"""
        python3 scripts/summarize_candidate_scaffold_annotation.py \
          --candidate-list {input.candidates} \
          --candidate-genes {input.candidate_genes} \
          --annotation {input.annotation} \
          --metrics {input.metrics} \
          --output-genes {output.genes} \
          --output-scaffolds {output.scaffolds}
        """


rule origin_kraken_candidate_scaffolds:
    input:
        fasta=ORIGIN_FASTA
    output:
        report=ORIGIN_KRAKEN_REPORT,
        classified=ORIGIN_KRAKEN_OUTPUT
    params:
        db=ORIGIN_KRAKEN_DB,
        confidence=ORIGIN_KRAKEN_CONFIDENCE,
        min_hit_groups=ORIGIN_KRAKEN_MIN_HIT_GROUPS,
        conda_module=ORIGIN_CONDA_MODULE,
        conda_env=ORIGIN_CONDA_ENV,
        cmd=ORIGIN_KRAKEN_CMD
    threads: 16
    shell:
        r"""
        set -euo pipefail
        if command -v module >/dev/null 2>&1; then
          module purge
          if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
        fi
        if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
        mkdir -p $(dirname {output.report})
        {params.cmd} \
          --db {params.db} \
          --threads {threads} \
          --confidence {params.confidence} \
          --minimum-hit-groups {params.min_hit_groups} \
          --report {output.report} \
          --output {output.classified} \
          {input.fasta}
        """


rule origin_summarize_scaffold_taxonomy:
    input:
        metrics=ORIGIN_SEQUENCE_METRICS,
        classified=ORIGIN_KRAKEN_OUTPUT,
        report=ORIGIN_KRAKEN_REPORT,
        nodes=ORIGIN_KRAKEN_DB + "/nodes.dmp",
        names=ORIGIN_KRAKEN_DB + "/names.dmp"
    output:
        ORIGIN_TAXONOMY
    shell:
        r"""
        python3 scripts/summarize_scaffold_kraken.py \
          --metrics {input.metrics} \
          --kraken-output {input.classified} \
          --kraken-report {input.report} \
          --nodes {input.nodes} \
          --names {input.names} \
          --output {output}
        """


rule origin_align_cross_species:
    input:
        query=ORIGIN_FASTA,
        target=lambda wildcards: ORIGIN_CROSS_SPECIES_REFERENCES[wildcards.species]
    output:
        ORIGIN_CROSS_SPECIES_DIR + "/{species}.paf"
    params:
        conda_module=ORIGIN_CONDA_MODULE,
        conda_env=ORIGIN_CONDA_ENV,
        cmd=ORIGIN_MINIMAP2_CMD
    threads: 16
    shell:
        r"""
        set -euo pipefail
        if command -v module >/dev/null 2>&1; then
          module purge
          if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
        fi
        if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
        mkdir -p $(dirname {output})
        {params.cmd} \
          -x asm20 \
          --secondary=no \
          -t {threads} \
          {input.target} \
          {input.query} \
          > {output}
        if [ ! -s {output} ]; then
          printf '# no alignments\n' > {output}
        fi
        """


rule origin_summarize_cross_species_homology:
    input:
        candidates=ORIGIN_STAGED_LIST,
        alignments=origin_cross_species_pafs()
    output:
        ORIGIN_CROSS_SPECIES_SUMMARY
    params:
        alignment_args=(
            " ".join(
                f"--alignment {species}={ORIGIN_CROSS_SPECIES_DIR}/{species}.paf"
                for species in sorted(ORIGIN_CROSS_SPECIES_REFERENCES)
            )
            if ORIGIN_RUN_CROSS_SPECIES
            else ""
        )
    shell:
        r"""
        python3 scripts/summarize_cross_species_scaffold_homology.py \
          --candidate-list {input.candidates} \
          {params.alignment_args} \
          --output {output}
        """


rule origin_build_evidence_table:
    input:
        provenance=ORIGIN_STAGED_PROVENANCE,
        annotation=ORIGIN_ANNOTATION_SUMMARY,
        taxonomy=ORIGIN_TAXONOMY,
        cross_species=ORIGIN_CROSS_SPECIES_SUMMARY
    output:
        evidence=ORIGIN_EVIDENCE,
        follow_up=ORIGIN_FOLLOW_UP
    shell:
        r"""
        python3 scripts/build_scaffold_origin_audit.py \
          --provenance {input.provenance} \
          --annotation-summary {input.annotation} \
          --taxonomy {input.taxonomy} \
          --cross-species {input.cross_species} \
          --output {output.evidence} \
          --follow-up {output.follow_up}
        """


rule origin_prepare_local_transfer_bundle:
    input:
        scaffold_origin_transfer_files()
    output:
        manifest=ORIGIN_TRANSFER_MANIFEST,
        archive=ORIGIN_TRANSFER_BUNDLE
    params:
        run_dir=ORIGIN_DIR,
        input_args=" ".join(
            f"--input {path}" for path in scaffold_origin_transfer_files()
        )
    shell:
        r"""
        python3 scripts/make_scaffold_origin_transfer_bundle.py \
          --run-dir {params.run_dir} \
          {params.input_args} \
          --manifest {output.manifest} \
          --archive {output.archive}
        """
