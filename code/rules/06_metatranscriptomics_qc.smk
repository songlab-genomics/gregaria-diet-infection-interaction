### =================================================================
### METATRANSCRIPTOMICS + CROSS-SAMPLE RNA-SEQ QC
### =================================================================
###
### Purpose:
###   Reuse completed STAR alignments to summarize host RNA-seq QC and classify
###   host-unmapped read pairs for bacteria, fungi, viruses, and other taxa.
###
### Output:
###   New files are written under METATX_DIR. Existing FASTQs, STAR BAMs,
###   featureCounts tables, and DESeq2 outputs are not modified.

import glob
import os


def metatx_project_dir(species):
    return f"{METATX_PROJECT_BASE}/{species}-timecourse"


def metatx_workdir(species):
    return f"{metatx_project_dir(species)}/data"


def metatx_metadata_file(species):
    candidates = sorted(glob.glob(f"../data/list/{species}/alltissues_*_super_metadata2026.csv"))
    if not candidates:
        candidates = sorted(glob.glob(f"../data/list/{species}/alltissues_*_super_metadata.csv"))
    return candidates[0] if candidates else ""


def metatx_discover_tissues(species):
    star_dirs = glob.glob(f"{metatx_workdir(species)}/02-*-star")
    tissues = []
    for path in star_dirs:
        name = os.path.basename(path)
        if name.startswith("02-") and name.endswith("-star"):
            tissues.append(name.replace("02-", "").replace("-star", ""))
    tissues = sorted(set(tissues))
    if METATX_TISSUES:
        tissues = [t for t in tissues if t in set(METATX_TISSUES)]
    return tissues


def metatx_sample_records():
    records = []
    for species in METATX_SPECIES:
        metadata_file = metatx_metadata_file(species)
        rrna_list = f"../data/list/excluded_loci/{species}_rrna_list.txt"
        for tissue in metatx_discover_tissues(species):
            for bam in sorted(glob.glob(f"{metatx_workdir(species)}/02-{tissue}-star/*_Aligned.sortedByCoord.out.bam")):
                sample = os.path.basename(bam).replace("_Aligned.sortedByCoord.out.bam", "")
                records.append({
                    "species": species,
                    "tissue": tissue,
                    "sample": sample,
                    "workdir": metatx_workdir(species),
                    "trimmed_r1": f"{metatx_workdir(species)}/01-{tissue}-trimmed-fastp/{sample}_1.trimmed.fastq.gz",
                    "trimmed_r2": f"{metatx_workdir(species)}/01-{tissue}-trimmed-fastp/{sample}_2.trimmed.fastq.gz",
                    "fastp_json": f"{metatx_workdir(species)}/01-{tissue}-trimmed-fastp/TrimQC/{sample}_fastp.json",
                    "bam": bam,
                    "bam_csi": bam + ".csi",
                    "star_log": f"{metatx_workdir(species)}/02-{tissue}-star/{sample}_Log.final.out",
                    "reads_per_gene": f"{metatx_workdir(species)}/02-{tissue}-star/{sample}_ReadsPerGene.out.tab",
                    "featurecounts_counts": f"{metatx_workdir(species)}/03-{tissue}-DESeq2/{sample}_featurecounts.txt",
                    "featurecounts_summary": f"{metatx_workdir(species)}/04-{tissue}-featurecounts/{sample}_counts.txt.summary",
                    "rrna_list": rrna_list,
                    "metadata_file": metadata_file,
                })
    return records


METATX_RECORDS = metatx_sample_records()
METATX_SPECIES_WC = [r["species"] for r in METATX_RECORDS]
METATX_TISSUES_WC = [r["tissue"] for r in METATX_RECORDS]
METATX_SAMPLES_WC = [r["sample"] for r in METATX_RECORDS]
METATX_RECORD_BY_KEY = {(r["species"], r["tissue"], r["sample"]): r for r in METATX_RECORDS}

print("Detected metatranscriptomics samples:")
for sp in sorted(set(METATX_SPECIES_WC)):
    n_sp = sum(1 for r in METATX_RECORDS if r["species"] == sp)
    tissues = sorted({r["tissue"] for r in METATX_RECORDS if r["species"] == sp})
    print(f"  {sp}: {n_sp} samples across {', '.join(tissues) if tissues else 'no tissues'}")


METATX_TABLE_DIR = METATX_DIR + "/tables"
METATX_MANIFEST = METATX_TABLE_DIR + "/metatx_sample_manifest.tsv"
METATX_SAMPLE_QC = METATX_TABLE_DIR + "/metatx_sample_qc.tsv"
METATX_TAXA_LONG = METATX_TABLE_DIR + "/metatx_taxa_long.tsv"
METATX_ORGANISM_SUMMARY = METATX_TABLE_DIR + "/metatx_organism_group_summary.tsv"
METATX_TOP_TAXA = METATX_TABLE_DIR + "/metatx_top_taxa.tsv"


def metatx_path(kind, wildcards):
    r = METATX_RECORD_BY_KEY[(wildcards.species, wildcards.tissue, wildcards.sample)]
    return r[kind]


def metatx_out(pattern):
    return expand(
        pattern,
        zip,
        species = METATX_SPECIES_WC,
        tissue = METATX_TISSUES_WC,
        sample = METATX_SAMPLES_WC,
    )


def metatx_host_unmapped_stats():
    return metatx_out(METATX_DIR + "/{species}/{tissue}/01-host-unmapped/{sample}.host_unmapped.stats.tsv")


def metatx_flagstat_outputs():
    return metatx_out(METATX_DIR + "/{species}/{tissue}/00-host-qc/{sample}.flagstat.txt")


def metatx_markdup_outputs():
    return metatx_out(METATX_DIR + "/{species}/{tissue}/00-host-qc/{sample}.markdup.metrics.txt")


def metatx_kraken2_reports():
    if not METATX_RUN_KRAKEN2:
        return []
    return metatx_out(METATX_DIR + "/{species}/{tissue}/02-kraken2/{sample}.kraken2.report")


def metatx_kaiju_summaries():
    if not METATX_RUN_KAIJU:
        return []
    return metatx_out(METATX_DIR + "/{species}/{tissue}/03-kaiju/{sample}.kaiju.summary.tsv")


def metatx_final_outputs():
    return [
        METATX_MANIFEST,
        METATX_SAMPLE_QC,
        METATX_TAXA_LONG,
        METATX_ORGANISM_SUMMARY,
        METATX_TOP_TAXA,
    ]


rule metatranscriptomics_qc_all:
    input:
        metatx_final_outputs()


rule metatx_sample_manifest:
    output:
        manifest = METATX_MANIFEST
    run:
        import csv
        os.makedirs(os.path.dirname(output.manifest), exist_ok=True)
        fieldnames = [
            "species", "tissue", "sample", "workdir", "trimmed_r1", "trimmed_r2",
            "fastp_json", "bam", "bam_csi", "star_log", "reads_per_gene",
            "featurecounts_counts", "featurecounts_summary", "rrna_list",
            "metadata_file", "flagstat", "samtools_stats", "markdup_metrics",
            "host_unmapped_r1", "host_unmapped_r2", "host_unmapped_stats",
            "kraken2_report", "kraken2_output", "kaiju_output", "kaiju_summary",
        ]
        with open(output.manifest, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for record in METATX_RECORDS:
                sp, tissue, sample = record["species"], record["tissue"], record["sample"]
                row = dict(record)
                row.update({
                    "flagstat": METATX_DIR + f"/{sp}/{tissue}/00-host-qc/{sample}.flagstat.txt",
                    "samtools_stats": METATX_DIR + f"/{sp}/{tissue}/00-host-qc/{sample}.samtools_stats.txt",
                    "markdup_metrics": METATX_DIR + f"/{sp}/{tissue}/00-host-qc/{sample}.markdup.metrics.txt",
                    "host_unmapped_r1": METATX_DIR + f"/{sp}/{tissue}/01-host-unmapped/{sample}_host_unmapped_R1.fastq.gz",
                    "host_unmapped_r2": METATX_DIR + f"/{sp}/{tissue}/01-host-unmapped/{sample}_host_unmapped_R2.fastq.gz",
                    "host_unmapped_stats": METATX_DIR + f"/{sp}/{tissue}/01-host-unmapped/{sample}.host_unmapped.stats.tsv",
                    "kraken2_report": METATX_DIR + f"/{sp}/{tissue}/02-kraken2/{sample}.kraken2.report",
                    "kraken2_output": METATX_DIR + f"/{sp}/{tissue}/02-kraken2/{sample}.kraken2.output",
                    "kaiju_output": METATX_DIR + f"/{sp}/{tissue}/03-kaiju/{sample}.kaiju.out",
                    "kaiju_summary": METATX_DIR + f"/{sp}/{tissue}/03-kaiju/{sample}.kaiju.summary.tsv",
                })
                writer.writerow({key: row.get(key, "") for key in fieldnames})


rule metatx_samtools_qc:
    input:
        bam = lambda wildcards: metatx_path("bam", wildcards)
    output:
        flagstat = METATX_DIR + "/{species}/{tissue}/00-host-qc/{sample}.flagstat.txt",
        stats = METATX_DIR + "/{species}/{tissue}/00-host-qc/{sample}.samtools_stats.txt"
    params:
        samtools_module = SAMTOOLS_MODULE,
        conda_module = METATX_CONDA_MODULE,
        conda_env = METATX_CONDA_ENV,
        tmp_prefix = METATX_DIR + "/{species}/{tissue}/01-host-unmapped/{sample}.collate_tmp"
    threads: 4
    log:
        stdout = METATX_DIR + "/{species}/{tissue}/logs/metatx_samtools_qc/{sample}.out",
        stderr = METATX_DIR + "/{species}/{tissue}/logs/metatx_samtools_qc/{sample}.err"
    shell:
        r"""
        mkdir -p $(dirname {output.flagstat}) $(dirname {log.stdout})
        {{
          if command -v module >/dev/null 2>&1; then
            module purge
            if [ -n "{params.samtools_module}" ]; then module load {params.samtools_module}; fi
            if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
          fi
          if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
          samtools flagstat -@ {threads} {input.bam} > {output.flagstat}
          samtools stats -@ {threads} {input.bam} > {output.stats}
        }} > {log.stdout} 2> {log.stderr}
        """


rule metatx_markdup_metrics:
    input:
        bam = lambda wildcards: metatx_path("bam", wildcards)
    output:
        metrics = METATX_DIR + "/{species}/{tissue}/00-host-qc/{sample}.markdup.metrics.txt"
    params:
        samtools_module = SAMTOOLS_MODULE,
        conda_module = METATX_CONDA_MODULE,
        conda_env = METATX_CONDA_ENV,
        tmp_prefix = METATX_DIR + "/{species}/{tissue}/00-host-qc/{sample}.collate_tmp"
    threads: 4
    log:
        stdout = METATX_DIR + "/{species}/{tissue}/logs/metatx_markdup_metrics/{sample}.out",
        stderr = METATX_DIR + "/{species}/{tissue}/logs/metatx_markdup_metrics/{sample}.err"
    shell:
        r"""
        mkdir -p $(dirname {output.metrics}) $(dirname {log.stdout})
        {{
          if command -v module >/dev/null 2>&1; then
            module purge
            if [ -n "{params.samtools_module}" ]; then module load {params.samtools_module}; fi
            if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
          fi
          if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
          samtools collate -@ {threads} -u -O {input.bam} {params.tmp_prefix} \
            | samtools fixmate -@ {threads} -m -u - - \
            | samtools sort -@ {threads} -u - \
            | samtools markdup -@ {threads} -s -f {output.metrics} - /dev/null
        }} > {log.stdout} 2> {log.stderr}
        """


rule metatx_export_host_unmapped_fastq:
    input:
        bam = lambda wildcards: metatx_path("bam", wildcards)
    output:
        r1 = METATX_DIR + "/{species}/{tissue}/01-host-unmapped/{sample}_host_unmapped_R1.fastq.gz",
        r2 = METATX_DIR + "/{species}/{tissue}/01-host-unmapped/{sample}_host_unmapped_R2.fastq.gz",
        stats = METATX_DIR + "/{species}/{tissue}/01-host-unmapped/{sample}.host_unmapped.stats.tsv"
    params:
        samtools_module = SAMTOOLS_MODULE,
        conda_module = METATX_CONDA_MODULE,
        conda_env = METATX_CONDA_ENV
    threads: 4
    log:
        stdout = METATX_DIR + "/{species}/{tissue}/logs/metatx_export_host_unmapped/{sample}.out",
        stderr = METATX_DIR + "/{species}/{tissue}/logs/metatx_export_host_unmapped/{sample}.err"
    shell:
        r"""
        mkdir -p $(dirname {output.r1}) $(dirname {log.stdout})
        {{
          if command -v module >/dev/null 2>&1; then
            module purge
            if [ -n "{params.samtools_module}" ]; then module load {params.samtools_module}; fi
            if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
          fi
          if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
          samtools collate -@ {threads} -u -O {input.bam} {params.tmp_prefix} \
            | samtools fastq -@ {threads} -f 12 -F 256 -n \
                -1 >(gzip -c > {output.r1}) \
                -2 >(gzip -c > {output.r2}) \
                -0 /dev/null -s /dev/null -
          pairs=$(gzip -cd {output.r1} | awk 'END {{print int(NR/4)}}')
          printf "SampleName\tHost_unmapped_read_pairs\n{wildcards.sample}\t%s\n" "$pairs" > {output.stats}
        }} > {log.stdout} 2> {log.stderr}
        """


rule metatx_kraken2_classify:
    input:
        r1 = METATX_DIR + "/{species}/{tissue}/01-host-unmapped/{sample}_host_unmapped_R1.fastq.gz",
        r2 = METATX_DIR + "/{species}/{tissue}/01-host-unmapped/{sample}_host_unmapped_R2.fastq.gz",
        stats = METATX_DIR + "/{species}/{tissue}/01-host-unmapped/{sample}.host_unmapped.stats.tsv"
    output:
        report = METATX_DIR + "/{species}/{tissue}/02-kraken2/{sample}.kraken2.report",
        classified = METATX_DIR + "/{species}/{tissue}/02-kraken2/{sample}.kraken2.output"
    params:
        module = KRAKEN2_MODULE,
        conda_module = METATX_CONDA_MODULE,
        conda_env = METATX_CONDA_ENV,
        cmd = KRAKEN2_CMD,
        db = KRAKEN2_DB
    threads: 16
    log:
        stdout = METATX_DIR + "/{species}/{tissue}/logs/metatx_kraken2/{sample}.out",
        stderr = METATX_DIR + "/{species}/{tissue}/logs/metatx_kraken2/{sample}.err"
    shell:
        r"""
        mkdir -p $(dirname {output.report}) $(dirname {log.stdout})
        {{
          if command -v module >/dev/null 2>&1; then
            module purge
            if [ -n "{params.module}" ]; then module load {params.module}; fi
            if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
          fi
          if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
          pairs=$(awk 'NR==2 {{print $2}}' {input.stats})
          if [ "${{pairs:-0}}" -eq 0 ]; then
            printf "100.00\t0\t0\tU\t0\tunclassified\n" > {output.report}
            : > {output.classified}
          else
            {params.cmd} \
              --db {params.db} \
              --threads {threads} \
              --paired \
              --gzip-compressed \
              --report {output.report} \
              --output {output.classified} \
              {input.r1} {input.r2}
          fi
        }} > {log.stdout} 2> {log.stderr}
        """


rule metatx_kaiju_classify:
    input:
        r1 = METATX_DIR + "/{species}/{tissue}/01-host-unmapped/{sample}_host_unmapped_R1.fastq.gz",
        r2 = METATX_DIR + "/{species}/{tissue}/01-host-unmapped/{sample}_host_unmapped_R2.fastq.gz",
        stats = METATX_DIR + "/{species}/{tissue}/01-host-unmapped/{sample}.host_unmapped.stats.tsv"
    output:
        out = METATX_DIR + "/{species}/{tissue}/03-kaiju/{sample}.kaiju.out",
        summary = METATX_DIR + "/{species}/{tissue}/03-kaiju/{sample}.kaiju.summary.tsv"
    params:
        module = KAIJU_MODULE,
        conda_module = METATX_CONDA_MODULE,
        conda_env = METATX_CONDA_ENV,
        kaiju = KAIJU_CMD,
        kaiju2table = KAIJU2TABLE_CMD,
        fmi = KAIJU_DB_FMI,
        nodes = KAIJU_DB_NODES,
        names = KAIJU_DB_NAMES,
        taxon_level = KAIJU_TAXON_LEVEL
    threads: 16
    log:
        stdout = METATX_DIR + "/{species}/{tissue}/logs/metatx_kaiju/{sample}.out",
        stderr = METATX_DIR + "/{species}/{tissue}/logs/metatx_kaiju/{sample}.err"
    shell:
        r"""
        mkdir -p $(dirname {output.out}) $(dirname {log.stdout})
        {{
          if command -v module >/dev/null 2>&1; then
            module purge
            if [ -n "{params.module}" ]; then module load {params.module}; fi
            if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
          fi
          if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi
          pairs=$(awk 'NR==2 {{print $2}}' {input.stats})
          if [ "${{pairs:-0}}" -eq 0 ]; then
            : > {output.out}
            printf "file\tpercent\treads\ttaxon_id\ttaxon_name\n" > {output.summary}
          else
            {params.kaiju} -z {threads} \
              -t {params.nodes} \
              -f {params.fmi} \
              -i {input.r1} \
              -j {input.r2} \
              -o {output.out}
            {params.kaiju2table} \
              -t {params.nodes} \
              -n {params.names} \
              -r {params.taxon_level} \
              -o {output.summary} \
              {output.out}
          fi
        }} > {log.stdout} 2> {log.stderr}
        """


rule metatx_summarize_all:
    input:
        manifest = METATX_MANIFEST,
        host_unmapped = metatx_host_unmapped_stats(),
        flagstat = metatx_flagstat_outputs(),
        markdup = metatx_markdup_outputs(),
        kraken2 = metatx_kraken2_reports(),
        kaiju = metatx_kaiju_summaries()
    output:
        sample_qc = METATX_SAMPLE_QC,
        taxa_long = METATX_TAXA_LONG,
        organism_summary = METATX_ORGANISM_SUMMARY,
        top_taxa = METATX_TOP_TAXA
    params:
        script = "../code/scripts/summarize_metatranscriptomics_qc.py"
    log:
        stdout = METATX_DIR + "/logs/metatx_summarize_all.out",
        stderr = METATX_DIR + "/logs/metatx_summarize_all.err"
    shell:
        r"""
        mkdir -p $(dirname {output.sample_qc}) $(dirname {log.stdout})
        {{
          python {params.script} \
            --manifest {input.manifest} \
            --sample-qc {output.sample_qc} \
            --taxa-long {output.taxa_long} \
            --organism-summary {output.organism_summary} \
            --top-taxa {output.top_taxa}
        }} > {log.stdout} 2> {log.stderr}
        """
