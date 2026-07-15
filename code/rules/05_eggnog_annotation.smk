### =================================================================
### EGGNOG-MAPPER FUNCTIONAL ANNOTATION
### =================================================================
###
### Purpose:
###   Annotate the S. gregaria reference proteins from scratch with eggNOG-
###   mapper, then convert protein-level annotations back to LOCID/GeneID so
###   the workflowR DEG enrichment pages can use GO and KEGG/KO gene sets.
###
### Output:
###   New files are written under EGGNOG_DIR. Existing read counts, STAR BAMs,
###   and DESeq2 outputs are not touched.

EGGNOG_ANNOTATIONS = EGGNOG_DIR + f"/{EGGNOG_PREFIX}.emapper.annotations"
EGGNOG_SEED_ORTHOLOGS = EGGNOG_DIR + f"/{EGGNOG_PREFIX}.emapper.seed_orthologs"
EGGNOG_GENE_ANNOTATIONS = EGGNOG_DIR + f"/postprocess/{EGGNOG_PREFIX}.emapper.gene_annotations.tsv"
EGGNOG_GENE_GO = EGGNOG_DIR + f"/postprocess/{EGGNOG_PREFIX}.emapper.gene_go.tsv"
EGGNOG_GENE_KEGG = EGGNOG_DIR + f"/postprocess/{EGGNOG_PREFIX}.emapper.gene_kegg.tsv"
EGGNOG_MANIFEST = EGGNOG_DIR + f"/postprocess/{EGGNOG_PREFIX}.emapper.manifest.tsv"
EGGNOG_DESERTLOCUSTR_DIR = EGGNOG_DIR + "/desertlocustr"
EGGNOG_DESERTLOCUSTR_MANIFEST = EGGNOG_DESERTLOCUSTR_DIR + f"/{EGGNOG_PREFIX}.desertlocustr_manifest.tsv"
EGGNOG_DESERTLOCUSTR_GO_TERM2GENE = EGGNOG_DESERTLOCUSTR_DIR + f"/{EGGNOG_PREFIX}.GO_TERM2GENE.tsv"
EGGNOG_DESERTLOCUSTR_GO_TERM2NAME = EGGNOG_DESERTLOCUSTR_DIR + f"/{EGGNOG_PREFIX}.GO_TERM2NAME.tsv"
EGGNOG_DESERTLOCUSTR_KEGG_KO = EGGNOG_DESERTLOCUSTR_DIR + f"/{EGGNOG_PREFIX}.KEGGTERM2LOC.tsv"
EGGNOG_DB = EGGNOG_DATA_DIR + "/eggnog.db"
EGGNOG_DMND = EGGNOG_DATA_DIR + "/eggnog_proteins.dmnd"
EGGNOG_TAXA_DB = EGGNOG_DATA_DIR + "/eggnog.taxa.db"


def eggnog_annotation_outputs():
    """Final files needed by the enrichment pages and manual annotation audit."""
    return [
        EGGNOG_ANNOTATIONS,
        EGGNOG_GENE_ANNOTATIONS,
        EGGNOG_GENE_GO,
        EGGNOG_GENE_KEGG,
        EGGNOG_MANIFEST,
        EGGNOG_DESERTLOCUSTR_MANIFEST,
        EGGNOG_DESERTLOCUSTR_GO_TERM2GENE,
        EGGNOG_DESERTLOCUSTR_GO_TERM2NAME,
        EGGNOG_DESERTLOCUSTR_KEGG_KO,
    ]


rule eggnog_annotations_all:
    input:
        eggnog_annotation_outputs()


rule eggnog_desertlocustr_all:
    input:
        eggnog_annotation_outputs()


rule eggnog_download_database:
    output:
        db = EGGNOG_DB,
        dmnd = EGGNOG_DMND,
        taxa = EGGNOG_TAXA_DB,
        marker = touch(EGGNOG_DATA_DIR + "/.eggnog_download_complete")
    threads: 8
    params:
        eggnog_module = EGGNOG_MODULE,
        eggnog_conda_module = EGGNOG_CONDA_MODULE,
        eggnog_conda_env = EGGNOG_CONDA_ENV,
        download_cmd = EGGNOG_DOWNLOAD_CMD,
        data_dir = EGGNOG_DATA_DIR
    log:
        stdout = EGGNOG_DIR + "/logs/eggnog_download_database.out",
        stderr = EGGNOG_DIR + "/logs/eggnog_download_database.err"
    shell:
        r"""
        mkdir -p {params.data_dir} $(dirname {log.stdout})
        {{
          if command -v module >/dev/null 2>&1; then
            module purge
            if [ -n "{params.eggnog_module}" ]; then module load {params.eggnog_module}; fi
            if [ -n "{params.eggnog_conda_module}" ]; then module load {params.eggnog_conda_module}; fi
          fi
          if [ -n "{params.eggnog_conda_env}" ]; then
            if command -v conda >/dev/null 2>&1; then
              eval "$(conda shell.bash hook)"
              conda activate {params.eggnog_conda_env}
            elif command -v mamba >/dev/null 2>&1; then
              eval "$(mamba shell hook --shell bash)"
              mamba activate {params.eggnog_conda_env}
            else
              source activate {params.eggnog_conda_env}
            fi
          fi

          command -v {params.download_cmd}

          if [ ! -s "{output.db}" ] || [ ! -s "{output.dmnd}" ] || [ ! -s "{output.taxa}" ]; then
            stale_dir="{params.data_dir}/incomplete_download_$(date +%Y%m%d_%H%M%S)"
            mkdir -p "$stale_dir"
            # Keep partial archives for audit, but move them out of the active DB directory.
            mv {params.data_dir}/*.gz "$stale_dir"/ 2>/dev/null || true
            mv {params.data_dir}/*.tar.gz "$stale_dir"/ 2>/dev/null || true
            mv {params.data_dir}/.eggnog_download_complete* "$stale_dir"/ 2>/dev/null || true
            echo "Moved incomplete eggNOG archives, if any, to: $stale_dir"
          fi

          # This is a large one-time download. Keep it outside the repo.
          {params.download_cmd} -P -y --data_dir {params.data_dir}
          test -s {output.db}
          test -s {output.dmnd}
          test -s {output.taxa}
        }} > {log.stdout} 2> {log.stderr}
        """


rule eggnog_preflight:
    input:
        protein = EGGNOG_PROTEIN_FASTA,
        protein2gene = EGGNOG_PROTEIN2GENE,
        eggnog_db = EGGNOG_DB,
        eggnog_dmnd = EGGNOG_DMND,
        eggnog_taxa = EGGNOG_TAXA_DB
    output:
        ok = EGGNOG_DIR + "/logs/eggnog_preflight.ok"
    params:
        eggnog_module = EGGNOG_MODULE,
        eggnog_conda_module = EGGNOG_CONDA_MODULE,
        eggnog_conda_env = EGGNOG_CONDA_ENV,
        eggnog_cmd = EGGNOG_CMD,
        data_dir = EGGNOG_DATA_DIR
    log:
        stdout = EGGNOG_DIR + "/logs/eggnog_preflight.out",
        stderr = EGGNOG_DIR + "/logs/eggnog_preflight.err"
    shell:
        r"""
        mkdir -p $(dirname {output.ok}) $(dirname {log.stdout})
        {{
          if command -v module >/dev/null 2>&1; then
            module purge
            if [ -n "{params.eggnog_module}" ]; then module load {params.eggnog_module}; fi
            if [ -n "{params.eggnog_conda_module}" ]; then module load {params.eggnog_conda_module}; fi
          fi
          if [ -n "{params.eggnog_conda_env}" ]; then
            if command -v conda >/dev/null 2>&1; then
              eval "$(conda shell.bash hook)"
              conda activate {params.eggnog_conda_env}
            elif command -v mamba >/dev/null 2>&1; then
              eval "$(mamba shell hook --shell bash)"
              mamba activate {params.eggnog_conda_env}
            else
              source activate {params.eggnog_conda_env}
            fi
          fi

          echo "Checking eggNOG-mapper command"
          command -v {params.eggnog_cmd}

          echo "Checking protein FASTA: {input.protein}"
          test -s {input.protein}

          echo "Checking protein-to-gene map: {input.protein2gene}"
          test -s {input.protein2gene}

          echo "Checking eggNOG database directory: {params.data_dir}"
          test -d {params.data_dir}
          test -s {input.eggnog_db}
          test -s {input.eggnog_dmnd}
          test -s {input.eggnog_taxa}

          date > {output.ok}
        }} > {log.stdout} 2> {log.stderr}
        """


rule eggnog_run_emapper:
    input:
        protein = EGGNOG_PROTEIN_FASTA,
        ok = EGGNOG_DIR + "/logs/eggnog_preflight.ok"
    output:
        annotations = EGGNOG_ANNOTATIONS
    threads: 16
    params:
        eggnog_module = EGGNOG_MODULE,
        eggnog_conda_module = EGGNOG_CONDA_MODULE,
        eggnog_conda_env = EGGNOG_CONDA_ENV,
        eggnog_cmd = EGGNOG_CMD,
        data_dir = EGGNOG_DATA_DIR,
        prefix = EGGNOG_PREFIX,
        out_dir = EGGNOG_DIR,
        go_evidence = EGGNOG_GO_EVIDENCE,
        target_orthologs = EGGNOG_TARGET_ORTHOLOGS,
        tax_scope = EGGNOG_TAX_SCOPE,
        tmp_dir = EGGNOG_DIR + "/tmp"
    log:
        stdout = EGGNOG_DIR + "/logs/eggnog_run_emapper.out",
        stderr = EGGNOG_DIR + "/logs/eggnog_run_emapper.err"
    shell:
        r"""
        mkdir -p {params.out_dir} {params.tmp_dir} $(dirname {log.stdout})
        {{
          if command -v module >/dev/null 2>&1; then
            module purge
            if [ -n "{params.eggnog_module}" ]; then module load {params.eggnog_module}; fi
            if [ -n "{params.eggnog_conda_module}" ]; then module load {params.eggnog_conda_module}; fi
          fi
          if [ -n "{params.eggnog_conda_env}" ]; then
            if command -v conda >/dev/null 2>&1; then
              eval "$(conda shell.bash hook)"
              conda activate {params.eggnog_conda_env}
            elif command -v mamba >/dev/null 2>&1; then
              eval "$(mamba shell hook --shell bash)"
              mamba activate {params.eggnog_conda_env}
            else
              source activate {params.eggnog_conda_env}
            fi
          fi

          tax_scope_arg=""
          if [ -n "{params.tax_scope}" ]; then tax_scope_arg="--tax_scope {params.tax_scope}"; fi

          # --override is safe here because EGGNOG_DIR contains the run ID.
          {params.eggnog_cmd} \
            -i {input.protein} \
            --itype proteins \
            --data_dir {params.data_dir} \
            --cpu {threads} \
            --go_evidence {params.go_evidence} \
            --target_orthologs {params.target_orthologs} \
            $tax_scope_arg \
            --output {params.prefix} \
            --output_dir {params.out_dir} \
            --temp_dir {params.tmp_dir} \
            --override

          test -s {output.annotations}
        }} > {log.stdout} 2> {log.stderr}
        """


rule eggnog_build_desertlocustr_universe:
    input:
        gene_annotations = EGGNOG_GENE_ANNOTATIONS,
        gene_go = EGGNOG_GENE_GO,
        gene_kegg = EGGNOG_GENE_KEGG
    output:
        manifest = EGGNOG_DESERTLOCUSTR_MANIFEST,
        go_term2gene = EGGNOG_DESERTLOCUSTR_GO_TERM2GENE,
        go_term2name = EGGNOG_DESERTLOCUSTR_GO_TERM2NAME,
        kegg_ko = EGGNOG_DESERTLOCUSTR_KEGG_KO
    params:
        builder = "../code/scripts/build_desertlocustr_universe.py",
        out_dir = EGGNOG_DESERTLOCUSTR_DIR,
        prefix = EGGNOG_PREFIX,
        species_label = EGGNOG_SPECIES_LABEL,
        run_id = EGGNOG_RUN_ID,
        go_name_table = EGGNOG_GO_NAME_TABLE
    log:
        stdout = EGGNOG_DIR + "/logs/eggnog_build_desertlocustr_universe.out",
        stderr = EGGNOG_DIR + "/logs/eggnog_build_desertlocustr_universe.err"
    shell:
        r"""
        mkdir -p {params.out_dir} $(dirname {log.stdout})
        {{
          python {params.builder} \
            --gene-annotations {input.gene_annotations} \
            --gene-go {input.gene_go} \
            --gene-kegg {input.gene_kegg} \
            --go-name-table "{params.go_name_table}" \
            --out-dir {params.out_dir} \
            --prefix {params.prefix} \
            --species "{params.species_label}" \
            --run-id {params.run_id}
        }} > {log.stdout} 2> {log.stderr}
        """


rule eggnog_parse_annotations:
    input:
        annotations = EGGNOG_ANNOTATIONS,
        protein2gene = EGGNOG_PROTEIN2GENE
    output:
        gene_annotations = EGGNOG_GENE_ANNOTATIONS,
        gene_go = EGGNOG_GENE_GO,
        gene_kegg = EGGNOG_GENE_KEGG,
        manifest = EGGNOG_MANIFEST
    params:
        parser = "../code/scripts/parse_emapper_annotations.py",
        species_label = EGGNOG_SPECIES_LABEL,
        run_id = EGGNOG_RUN_ID,
        protein_fasta = EGGNOG_PROTEIN_FASTA,
        data_dir = EGGNOG_DATA_DIR
    log:
        stdout = EGGNOG_DIR + "/logs/eggnog_parse_annotations.out",
        stderr = EGGNOG_DIR + "/logs/eggnog_parse_annotations.err"
    shell:
        r"""
        mkdir -p $(dirname {output.gene_annotations}) $(dirname {log.stdout})
        {{
          python {params.parser} \
            --emapper {input.annotations} \
            --protein2gene {input.protein2gene} \
            --species "{params.species_label}" \
            --gene-annotations {output.gene_annotations} \
            --gene-go {output.gene_go} \
            --gene-kegg {output.gene_kegg} \
            --manifest {output.manifest} \
            --run-id {params.run_id} \
            --protein-fasta {params.protein_fasta} \
            --eggnog-data-dir {params.data_dir}
        }} > {log.stdout} 2> {log.stderr}
        """
