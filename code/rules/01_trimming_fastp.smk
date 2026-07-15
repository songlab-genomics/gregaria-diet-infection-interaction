### =================================================================
###  DATA PREPROCESSING: READS FILTERING AND TRIMMING
### =================================================================

rule trimming_fastp:
    input:
        read1 = DATADIR + "/00-reads-{tissue}/{locust}_MERGE_1.fastq.gz",
        read2 = DATADIR + "/00-reads-{tissue}/{locust}_MERGE_2.fastq.gz"
    output:
        reportjson = WORKDIR + "/01-{tissue}-trimmed-fastp/TrimQC/{locust}_fastp.json",
        reporthtml = WORKDIR + "/01-{tissue}-trimmed-fastp/TrimQC/{locust}_fastp.html",
        tread1 = WORKDIR + "/01-{tissue}-trimmed-fastp/{locust}_1.trimmed.fastq.gz",
        tread2 = WORKDIR + "/01-{tissue}-trimmed-fastp/{locust}_2.trimmed.fastq.gz"
    params:
        fastp_module = FASTP_MODULE,
        conda_module = FASTP_CONDA_MODULE,
        conda_env = FASTP_CONDA_ENV,
        fastp_cmd = FASTP_CMD
    threads: 16
    shell:
        r"""
        if command -v module >/dev/null 2>&1; then
          module purge
          if [ -n "{params.fastp_module}" ]; then module load {params.fastp_module}; fi
          if [ -n "{params.conda_module}" ]; then module load {params.conda_module}; fi
        fi
        if [ -n "{params.conda_env}" ]; then source activate {params.conda_env}; fi

        mkdir -p $(dirname {output.reportjson})
        mkdir -p $(dirname {output.tread1})

        echo "Using fastp command: $(command -v {params.fastp_cmd})"
        {params.fastp_cmd} --version || true

        {params.fastp_cmd} --thread {threads} \
          --in1 {input.read1} \
          --in2 {input.read2} \
          --out1 {output.tread1} \
          --out2 {output.tread2} \
          --trim_front1 2 \
          --trim_front2 2 \
          --detect_adapter_for_pe \
          -l 50 \
          --json {output.reportjson} \
          --html {output.reporthtml}
        """

rule trimming_tgalore:
        input:
                read1 = DATADIR + "/00-{species}-reads/{locust}_1.fastq.gz",
                read2 = DATADIR + "/00-{species}-reads/{locust}_2.fastq.gz"
        output: 
                tread1 = WORKDIR + "/01-{species}-trimmed-tgalore/{locust}_1_val_1.fq.gz",
                tread2 = WORKDIR + "/01-{species}-trimmed-tgalore/{locust}_2_val_2.fq.gz"
        shell:
                """
                module purge
                module load trimgalore-0.6.10-gcc-12.1.0

                trim_galore --trim-n \\
                --cores 16 \\
                --quality 20 \\
                --clip_R1 2 \\
                --clip_R2 2 \\
                --nextera \\
                --output_dir {WORKDIR}/01-{wildcards.species}-trimmed-tgalore/ \\
                --fastqc \\
                --paired {input.read1} {input.read2}
                """    
