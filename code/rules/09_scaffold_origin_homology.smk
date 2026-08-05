HOMOLOGY_INPUT_DIR = HOMOLOGY_DIR + "/00-input"
HOMOLOGY_SEQUENCE_DIR = HOMOLOGY_DIR + "/01-sequences"
HOMOLOGY_BLASTN_DIR = HOMOLOGY_DIR + "/02-blastn-core-nt"
HOMOLOGY_DIAMOND_DIR = HOMOLOGY_DIR + "/03-diamond-nr-cluster-seq"
HOMOLOGY_EVIDENCE_DIR = HOMOLOGY_DIR + "/04-evidence"
HOMOLOGY_TRANSFER_DIR = HOMOLOGY_DIR + "/05-local-transfer"

HOMOLOGY_STAGED_SCAFFOLDS = HOMOLOGY_INPUT_DIR + "/candidate_scaffolds.txt"
HOMOLOGY_STAGED_GENES = HOMOLOGY_INPUT_DIR + "/candidate_genes.tsv"
HOMOLOGY_STAGED_PROVENANCE = HOMOLOGY_INPUT_DIR + "/candidate_provenance.tsv"
HOMOLOGY_INPUT_MANIFEST = HOMOLOGY_INPUT_DIR + "/input_manifest.tsv"

HOMOLOGY_SCAFFOLD_FASTA = HOMOLOGY_SEQUENCE_DIR + "/candidate_scaffolds.fna"
HOMOLOGY_SCAFFOLD_METRICS = (
    HOMOLOGY_SEQUENCE_DIR + "/candidate_scaffold_sequence_metrics.tsv"
)
HOMOLOGY_NUCLEOTIDE_QUERIES = (
    HOMOLOGY_SEQUENCE_DIR + "/candidate_gene_regions.fna"
)
HOMOLOGY_PROTEIN_QUERIES = HOMOLOGY_SEQUENCE_DIR + "/candidate_proteins.faa"
HOMOLOGY_QUERY_MANIFEST = HOMOLOGY_SEQUENCE_DIR + "/query_manifest.tsv"
HOMOLOGY_MISSING_QUERIES = HOMOLOGY_SEQUENCE_DIR + "/missing_protein_queries.tsv"

HOMOLOGY_BLASTN_HITS = HOMOLOGY_BLASTN_DIR + "/gene_regions_vs_core_nt.tsv"
HOMOLOGY_DIAMOND_HITS = (
    HOMOLOGY_DIAMOND_DIR + "/proteins_vs_nr_cluster_seq.tsv"
)

HOMOLOGY_TOP_HITS = HOMOLOGY_EVIDENCE_DIR + "/top_tier_taxonomic_hits.tsv"
HOMOLOGY_QUERY_SUMMARY = HOMOLOGY_EVIDENCE_DIR + "/query_homology_summary.tsv"
HOMOLOGY_GENE_SUMMARY = HOMOLOGY_EVIDENCE_DIR + "/gene_homology_summary.tsv"
HOMOLOGY_SCAFFOLD_SUMMARY = (
    HOMOLOGY_EVIDENCE_DIR + "/scaffold_homology_summary.tsv"
)
HOMOLOGY_NON_METAZOAN = (
    HOMOLOGY_EVIDENCE_DIR + "/non_metazoan_candidate_genes.tsv"
)
HOMOLOGY_TRANSFER_MANIFEST = HOMOLOGY_TRANSFER_DIR + "/transfer_manifest.tsv"
HOMOLOGY_TRANSFER_BUNDLE = (
    HOMOLOGY_TRANSFER_DIR + "/scaffold_homology_for_local.tar.gz"
)


def scaffold_origin_homology_transfer_files():
    return [
        HOMOLOGY_INPUT_MANIFEST,
        HOMOLOGY_STAGED_SCAFFOLDS,
        HOMOLOGY_STAGED_GENES,
        HOMOLOGY_STAGED_PROVENANCE,
        HOMOLOGY_SCAFFOLD_FASTA,
        HOMOLOGY_SCAFFOLD_METRICS,
        HOMOLOGY_NUCLEOTIDE_QUERIES,
        HOMOLOGY_PROTEIN_QUERIES,
        HOMOLOGY_QUERY_MANIFEST,
        HOMOLOGY_MISSING_QUERIES,
        HOMOLOGY_BLASTN_HITS,
        HOMOLOGY_DIAMOND_HITS,
        HOMOLOGY_TOP_HITS,
        HOMOLOGY_QUERY_SUMMARY,
        HOMOLOGY_GENE_SUMMARY,
        HOMOLOGY_SCAFFOLD_SUMMARY,
        HOMOLOGY_NON_METAZOAN,
    ]


def scaffold_origin_homology_outputs():
    return [
        HOMOLOGY_GENE_SUMMARY,
        HOMOLOGY_SCAFFOLD_SUMMARY,
        HOMOLOGY_NON_METAZOAN,
        HOMOLOGY_TRANSFER_MANIFEST,
        HOMOLOGY_TRANSFER_BUNDLE,
    ]


rule scaffold_origin_homology_hpc:
    input:
        scaffold_origin_homology_outputs()


rule homology_stage_inputs:
    input:
        scaffolds=HOMOLOGY_CANDIDATE_SCAFFOLDS,
        genes=HOMOLOGY_CANDIDATE_GENES,
        provenance=HOMOLOGY_CANDIDATE_PROVENANCE
    output:
        scaffolds=HOMOLOGY_STAGED_SCAFFOLDS,
        genes=HOMOLOGY_STAGED_GENES,
        provenance=HOMOLOGY_STAGED_PROVENANCE,
        manifest=HOMOLOGY_INPUT_MANIFEST
    shell:
        r"""
        python3 scripts/stage_scaffold_homology_inputs.py \
          --scaffolds {input.scaffolds} \
          --genes {input.genes} \
          --provenance {input.provenance} \
          --output-scaffolds {output.scaffolds} \
          --output-genes {output.genes} \
          --output-provenance {output.provenance} \
          --manifest {output.manifest}
        """


rule homology_extract_candidate_scaffolds:
    input:
        fasta=HOMOLOGY_HOST_FASTA,
        candidates=HOMOLOGY_STAGED_SCAFFOLDS
    output:
        fasta=HOMOLOGY_SCAFFOLD_FASTA,
        metrics=HOMOLOGY_SCAFFOLD_METRICS
    shell:
        r"""
        python3 scripts/extract_unplaced_scaffolds.py \
          --fasta {input.fasta} \
          --sequence-list {input.candidates} \
          --output-fasta {output.fasta} \
          --output-table {output.metrics}
        """


rule homology_prepare_gene_queries:
    input:
        genes=HOMOLOGY_STAGED_GENES,
        scaffolds=HOMOLOGY_SCAFFOLD_FASTA,
        annotation=HOMOLOGY_HOST_ANNOTATION,
        proteins=HOMOLOGY_HOST_PROTEIN_FASTA
    output:
        nucleotide=HOMOLOGY_NUCLEOTIDE_QUERIES,
        protein=HOMOLOGY_PROTEIN_QUERIES,
        manifest=HOMOLOGY_QUERY_MANIFEST,
        missing=HOMOLOGY_MISSING_QUERIES
    shell:
        r"""
        python3 scripts/prepare_scaffold_homology_queries.py \
          --candidate-genes {input.genes} \
          --candidate-scaffolds-fasta {input.scaffolds} \
          --annotation {input.annotation} \
          --protein-fasta {input.proteins} \
          --output-nucleotide {output.nucleotide} \
          --output-protein {output.protein} \
          --output-manifest {output.manifest} \
          --output-missing {output.missing}
        """


rule homology_blastn_core_nt:
    input:
        query=HOMOLOGY_NUCLEOTIDE_QUERIES,
        databases_ready=HOMOLOGY_DATABASE_READY
    output:
        HOMOLOGY_BLASTN_HITS
    params:
        database=HOMOLOGY_CORE_NT_DB,
        conda_module=HOMOLOGY_CONDA_MODULE,
        conda_env=HOMOLOGY_CONDA_ENV,
        cmd=HOMOLOGY_BLASTN_CMD,
        self_taxid=HOMOLOGY_SELF_TAXID,
        max_targets=HOMOLOGY_MAX_TARGET_SEQS
    threads: 32
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
          -query {input.query} \
          -db {params.database} \
          -task dc-megablast \
          -negative_taxids {params.self_taxid} \
          -evalue 1e-10 \
          -perc_identity 60 \
          -qcov_hsp_perc 20 \
          -max_target_seqs {params.max_targets} \
          -max_hsps 5 \
          -num_threads {threads} \
          -outfmt '6 qseqid qlen sseqid staxids pident length qcovhsp evalue bitscore stitle' \
          -out {output}
        touch {output}
        """


rule homology_diamond_nr_cluster_seq:
    input:
        query=HOMOLOGY_PROTEIN_QUERIES,
        databases_ready=HOMOLOGY_DATABASE_READY
    output:
        HOMOLOGY_DIAMOND_HITS
    params:
        database=HOMOLOGY_NR_CLUSTER_SEQ_DB,
        conda_module=HOMOLOGY_CONDA_MODULE,
        conda_env=HOMOLOGY_CONDA_ENV,
        cmd=HOMOLOGY_DIAMOND_CMD,
        self_taxid=HOMOLOGY_SELF_TAXID,
        max_targets=HOMOLOGY_MAX_TARGET_SEQS
    threads: 32
    shell:
        r"""
        set -euo pipefail
        if command -v module >/dev/null 2>&1; then
          module purge
          if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
        fi
        if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
        mkdir -p $(dirname {output})
        {params.cmd} blastp \
          --query {input.query} \
          --db {params.database} \
          --taxon-exclude {params.self_taxid} \
          --very-sensitive \
          --evalue 1e-5 \
          --query-cover 20 \
          --max-target-seqs {params.max_targets} \
          --threads {threads} \
          --outfmt 6 qseqid qlen sseqid staxids pident length qcovhsp evalue bitscore stitle \
          --out {output}
        touch {output}
        """


rule homology_summarize_taxonomic_hits:
    input:
        manifest=HOMOLOGY_QUERY_MANIFEST,
        genes=HOMOLOGY_STAGED_GENES,
        blastn=HOMOLOGY_BLASTN_HITS,
        diamond=HOMOLOGY_DIAMOND_HITS,
        nodes=HOMOLOGY_TAXONOMY_NODES,
        names=HOMOLOGY_TAXONOMY_NAMES
    output:
        top_hits=HOMOLOGY_TOP_HITS,
        queries=HOMOLOGY_QUERY_SUMMARY,
        genes=HOMOLOGY_GENE_SUMMARY,
        scaffolds=HOMOLOGY_SCAFFOLD_SUMMARY,
        non_metazoan=HOMOLOGY_NON_METAZOAN
    params:
        top_fraction=HOMOLOGY_TOP_FRACTION
    shell:
        r"""
        python3 scripts/summarize_scaffold_homology_hits.py \
          --manifest {input.manifest} \
          --candidate-genes {input.genes} \
          --blastn-hits {input.blastn} \
          --diamond-hits {input.diamond} \
          --nodes {input.nodes} \
          --names {input.names} \
          --exclude-taxid {params.self_taxid} \
          --top-fraction {params.top_fraction} \
          --output-top-hits {output.top_hits} \
          --output-query-summary {output.queries} \
          --output-gene-summary {output.genes} \
          --output-scaffold-summary {output.scaffolds} \
          --output-non-metazoan {output.non_metazoan}
        """


rule homology_prepare_local_transfer_bundle:
    input:
        scaffold_origin_homology_transfer_files()
    output:
        manifest=HOMOLOGY_TRANSFER_MANIFEST,
        archive=HOMOLOGY_TRANSFER_BUNDLE
    params:
        run_dir=HOMOLOGY_DIR,
        input_args=" ".join(
            f"--input {path}"
            for path in scaffold_origin_homology_transfer_files()
        )
    shell:
        r"""
        python3 scripts/make_scaffold_origin_transfer_bundle.py \
          --run-dir {params.run_dir} \
          {params.input_args} \
          --manifest {output.manifest} \
          --archive {output.archive}
        """
