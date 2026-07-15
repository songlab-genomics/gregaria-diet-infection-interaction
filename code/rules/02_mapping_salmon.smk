### =================================================================
### SALMON (index + quant)
### =================================================================

rule salmon_index:
    input:
        genome = REF_GENOME,
        transcriptome = REF_TRANSCRIPTOME
    output:
        decoy    = SALMON_DIR + "/decoys.txt",
        gentrome = SALMON_DIR + "/gentrome.fa",
        indexdir = directory(SALMON_INDEX_DIR)
    threads: 12
    shell:
        r"""
        module purge
        module load salmon-1.10.2-gcc-12.1.0

        mkdir -p {SALMON_DIR}
        mkdir -p {output.indexdir}

        # Decoys: contig names from genome fasta
        grep "^>" {input.genome} | cut -d " " -f 1 > {output.decoy}
        sed -i.bak -e 's/>//g' {output.decoy}

        # Gentrome: transcriptome + genome
        cat {input.transcriptome} {input.genome} > {output.gentrome}

        salmon index \
          -t {output.gentrome} \
          -d {output.decoy} \
          -p {threads} \
          -i {output.indexdir}
        """

rule salmon_quant:
    input:
        indexdir = SALMON_INDEX_DIR,
        trimmed_read1 = WORKDIR + "/01-{tissue}-trimmed-fastp/{locust}_1.trimmed.fastq.gz",
        trimmed_read2 = WORKDIR + "/01-{tissue}-trimmed-fastp/{locust}_2.trimmed.fastq.gz"
    output:
        outdir = directory(WORKDIR + "/02-{tissue}-salmon/{locust}_salmon_out")
    threads: 16
    shell:
        r"""
        module purge
        module load salmon-1.10.2-gcc-12.1.0

        mkdir -p $(dirname {output.outdir})

        salmon quant \
          -i {input.indexdir} \
          -l A \
          -p {threads} \
          -1 {input.trimmed_read1} \
          -2 {input.trimmed_read2} \
          -o {output.outdir}
        """
