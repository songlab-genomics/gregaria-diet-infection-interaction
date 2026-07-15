### =================================================================
### RNA-SEQ VARIANT CALLING
### =================================================================
###
### Purpose:
###   Call transcriptome-derived SNPs from the STAR genomic BAMs. These
###   variants are useful for exploratory population structure/genetic
###   variation when WGS is not available, but they only represent expressed
###   loci and can be influenced by tissue/time expression and RNA editing.
###
### Method:
###   This follows the GATK RNA-seq short-variant logic:
###     1. index the reference FASTA and create a GATK sequence dictionary
###     2. mark PCR/optical duplicates in coordinate-sorted STAR BAMs
###     3. split reads across intronic N-cigar operations
###     4. call one GVCF per library with HaplotypeCaller
###     5. jointly genotype samples per tissue
###     6. retain hard-filtered PASS SNPs for downstream diversity analyses

RNASEQ_VARIANT_REF = RNASEQ_VARIANT_DIR + "/reference/" + REF_GENOME.rsplit("/", 1)[-1]
RNASEQ_VARIANT_REF_DICT = RNASEQ_VARIANT_REF.rsplit(".", 1)[0] + ".dict"


def rnaseq_variant_sample_outputs():
    """Per-sample RNA-seq GVCFs expected from all discovered tissues."""
    return [
        RNASEQ_VARIANT_DIR + f"/{tissue}/03-gvcf/{locust}.rnaseq.g.vcf.gz"
        for tissue, locusts in LOCUSTS.items()
        for locust in locusts
    ]


def rnaseq_variant_tissue_outputs():
    """Final per-tissue PASS SNP VCFs."""
    return [
        RNASEQ_VARIANT_DIR + f"/{tissue}/06-filtered/{SPECIES}.{tissue}.rnaseq.pass.snps.vcf.gz"
        for tissue, locusts in LOCUSTS.items()
        if len(locusts) > 0
    ]


rule rnaseq_variants_all:
    input:
        rnaseq_variant_tissue_outputs()


rule rnaseq_variant_reference:
    input:
        ref = REF_GENOME
    output:
        ref = RNASEQ_VARIANT_REF,
        fai = RNASEQ_VARIANT_REF + ".fai",
        dict = RNASEQ_VARIANT_REF_DICT
    threads: 2
    params:
        samtools_module = SAMTOOLS_MODULE,
        gatk_module = GATK_MODULE,
        gatk_conda_module = GATK_CONDA_MODULE,
        gatk_conda_env = GATK_CONDA_ENV,
        gatk_cmd = GATK_CMD
    shell:
        r"""
        module purge
        if [ -n "{params.samtools_module}" ]; then module load {params.samtools_module}; fi
        if [ -n "{params.gatk_module}" ]; then module load {params.gatk_module}; fi
        if [ -n "{params.gatk_conda_module}" ]; then module load {params.gatk_conda_module}; fi
        if [ -n "{params.gatk_conda_env}" ]; then source activate {params.gatk_conda_env}; fi

        mkdir -p $(dirname {output.ref})

        # Stage a symlinked reference inside this run folder so /data reference files are not edited.
        ln -sf {input.ref} {output.ref}
        samtools faidx {output.ref}
        {params.gatk_cmd} CreateSequenceDictionary -R {output.ref} -O {output.dict}
        """


rule rnaseq_mark_duplicates:
    input:
        bam = WORKDIR + "/02-{tissue}-star/{locust}_Aligned.sortedByCoord.out.bam",
        csi = WORKDIR + "/02-{tissue}-star/{locust}_Aligned.sortedByCoord.out.bam.csi"
    output:
        bam = RNASEQ_VARIANT_DIR + "/{tissue}/01-markdup/{locust}.markdup.bam",
        csi = RNASEQ_VARIANT_DIR + "/{tissue}/01-markdup/{locust}.markdup.bam.csi",
        metrics = RNASEQ_VARIANT_DIR + "/{tissue}/01-markdup/{locust}.markdup.metrics.txt"
    log:
        stdout = RNASEQ_VARIANT_DIR + "/{tissue}/logs/rnaseq_mark_duplicates/{locust}.out",
        stderr = RNASEQ_VARIANT_DIR + "/{tissue}/logs/rnaseq_mark_duplicates/{locust}.err"
    threads: 4
    params:
        samtools_module = SAMTOOLS_MODULE,
        gatk_module = GATK_MODULE,
        gatk_conda_module = GATK_CONDA_MODULE,
        gatk_conda_env = GATK_CONDA_ENV,
        gatk_cmd = GATK_CMD,
        java_mem = "48g",
        tmp = RNASEQ_VARIANT_DIR + "/{tissue}/tmp/rnaseq_mark_duplicates/{locust}"
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) $(dirname {log.stdout}) {params.tmp}
        {{
          module purge
          if [ -n "{params.samtools_module}" ]; then module load {params.samtools_module}; fi
          if [ -n "{params.gatk_module}" ]; then module load {params.gatk_module}; fi
          if [ -n "{params.gatk_conda_module}" ]; then module load {params.gatk_conda_module}; fi
          if [ -n "{params.gatk_conda_env}" ]; then source activate {params.gatk_conda_env}; fi

          # Duplicates are marked, not removed, so later filtering can still inspect them.
          {params.gatk_cmd} --java-options "-Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp}" MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX false \
            --VALIDATION_STRINGENCY SILENT \
            --TMP_DIR {params.tmp}

          # Locust chromosomes/scaffolds can be large, so use CSI rather than BAI indexes.
          samtools index -@ {threads} -c {output.bam}
        }} > {log.stdout} 2> {log.stderr}
        """


rule rnaseq_split_ncigar:
    input:
        ref = RNASEQ_VARIANT_REF,
        dict = RNASEQ_VARIANT_REF_DICT,
        fai = RNASEQ_VARIANT_REF + ".fai",
        bam = RNASEQ_VARIANT_DIR + "/{tissue}/01-markdup/{locust}.markdup.bam",
        csi = RNASEQ_VARIANT_DIR + "/{tissue}/01-markdup/{locust}.markdup.bam.csi"
    output:
        bam = RNASEQ_VARIANT_DIR + "/{tissue}/02-split-ncigar/{locust}.splitncigar.bam",
        csi = RNASEQ_VARIANT_DIR + "/{tissue}/02-split-ncigar/{locust}.splitncigar.bam.csi"
    log:
        stdout = RNASEQ_VARIANT_DIR + "/{tissue}/logs/rnaseq_split_ncigar/{locust}.out",
        stderr = RNASEQ_VARIANT_DIR + "/{tissue}/logs/rnaseq_split_ncigar/{locust}.err"
    threads: 4
    params:
        samtools_module = SAMTOOLS_MODULE,
        gatk_module = GATK_MODULE,
        gatk_conda_module = GATK_CONDA_MODULE,
        gatk_conda_env = GATK_CONDA_ENV,
        gatk_cmd = GATK_CMD,
        java_mem = "48g",
        tmp = RNASEQ_VARIANT_DIR + "/{tissue}/tmp/rnaseq_split_ncigar/{locust}"
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) $(dirname {log.stdout}) {params.tmp}
        {{
          module purge
          if [ -n "{params.samtools_module}" ]; then module load {params.samtools_module}; fi
          if [ -n "{params.gatk_module}" ]; then module load {params.gatk_module}; fi
          if [ -n "{params.gatk_conda_module}" ]; then module load {params.gatk_conda_module}; fi
          if [ -n "{params.gatk_conda_env}" ]; then source activate {params.gatk_conda_env}; fi

          # SplitNCigarReads makes spliced RNA-seq alignments compatible with variant calling.
          # GATK's default BAM index can be BAI; we disable it and create CSI for large locust scaffolds.
          {params.gatk_cmd} --java-options "-Xmx{params.java_mem} -Djava.io.tmpdir={params.tmp}" SplitNCigarReads \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.bam} \
            --create-output-bam-index false \
            --tmp-dir {params.tmp}

          samtools index -@ {threads} -c {output.bam}
        }} > {log.stdout} 2> {log.stderr}
        """


rule rnaseq_haplotypecaller_gvcf:
    input:
        ref = RNASEQ_VARIANT_REF,
        dict = RNASEQ_VARIANT_REF_DICT,
        fai = RNASEQ_VARIANT_REF + ".fai",
        bam = RNASEQ_VARIANT_DIR + "/{tissue}/02-split-ncigar/{locust}.splitncigar.bam",
        csi = RNASEQ_VARIANT_DIR + "/{tissue}/02-split-ncigar/{locust}.splitncigar.bam.csi"
    output:
        gvcf = RNASEQ_VARIANT_DIR + "/{tissue}/03-gvcf/{locust}.rnaseq.g.vcf.gz",
        tbi = RNASEQ_VARIANT_DIR + "/{tissue}/03-gvcf/{locust}.rnaseq.g.vcf.gz.tbi"
    threads: 6
    params:
        gatk_module = GATK_MODULE,
        gatk_conda_module = GATK_CONDA_MODULE,
        gatk_conda_env = GATK_CONDA_ENV,
        gatk_cmd = GATK_CMD,
        java_mem = "24g"
    shell:
        r"""
        module purge
        if [ -n "{params.gatk_module}" ]; then module load {params.gatk_module}; fi
        if [ -n "{params.gatk_conda_module}" ]; then module load {params.gatk_conda_module}; fi
        if [ -n "{params.gatk_conda_env}" ]; then source activate {params.gatk_conda_env}; fi

        mkdir -p $(dirname {output.gvcf})

        # GVCF mode keeps per-sample evidence so samples can be jointly genotyped per tissue.
        {params.gatk_cmd} --java-options "-Xmx{params.java_mem}" HaplotypeCaller \
          -R {input.ref} \
          -I {input.bam} \
          -O {output.gvcf} \
          -ERC GVCF \
          --dont-use-soft-clipped-bases true \
          --standard-min-confidence-threshold-for-calling 20

        {params.gatk_cmd} IndexFeatureFile -I {output.gvcf}
        """


rule rnaseq_combine_gvcfs_by_tissue:
    input:
        ref = RNASEQ_VARIANT_REF,
        dict = RNASEQ_VARIANT_REF_DICT,
        fai = RNASEQ_VARIANT_REF + ".fai",
        gvcfs = lambda wildcards: expand(
            RNASEQ_VARIANT_DIR + "/{tissue}/03-gvcf/{locust}.rnaseq.g.vcf.gz",
            tissue = wildcards.tissue,
            locust = LOCUSTS[wildcards.tissue]
        )
    output:
        gvcf = RNASEQ_VARIANT_DIR + "/{tissue}/04-joint/" + SPECIES + ".{tissue}.combined.g.vcf.gz",
        tbi = RNASEQ_VARIANT_DIR + "/{tissue}/04-joint/" + SPECIES + ".{tissue}.combined.g.vcf.gz.tbi"
    threads: 4
    params:
        gatk_module = GATK_MODULE,
        gatk_conda_module = GATK_CONDA_MODULE,
        gatk_conda_env = GATK_CONDA_ENV,
        gatk_cmd = GATK_CMD,
        java_mem = "32g"
    shell:
        r"""
        module purge
        if [ -n "{params.gatk_module}" ]; then module load {params.gatk_module}; fi
        if [ -n "{params.gatk_conda_module}" ]; then module load {params.gatk_conda_module}; fi
        if [ -n "{params.gatk_conda_env}" ]; then source activate {params.gatk_conda_env}; fi

        mkdir -p $(dirname {output.gvcf})

        variant_args=""
        for gvcf in {input.gvcfs}; do
          variant_args="${{variant_args}} --variant ${{gvcf}}"
        done

        {params.gatk_cmd} --java-options "-Xmx{params.java_mem}" CombineGVCFs \
          -R {input.ref} \
          ${{variant_args}} \
          -O {output.gvcf}

        {params.gatk_cmd} IndexFeatureFile -I {output.gvcf}
        """


rule rnaseq_genotype_gvcfs_by_tissue:
    input:
        ref = RNASEQ_VARIANT_REF,
        dict = RNASEQ_VARIANT_REF_DICT,
        fai = RNASEQ_VARIANT_REF + ".fai",
        gvcf = RNASEQ_VARIANT_DIR + "/{tissue}/04-joint/" + SPECIES + ".{tissue}.combined.g.vcf.gz",
        tbi = RNASEQ_VARIANT_DIR + "/{tissue}/04-joint/" + SPECIES + ".{tissue}.combined.g.vcf.gz.tbi"
    output:
        vcf = RNASEQ_VARIANT_DIR + "/{tissue}/05-raw/" + SPECIES + ".{tissue}.rnaseq.raw.vcf.gz",
        tbi = RNASEQ_VARIANT_DIR + "/{tissue}/05-raw/" + SPECIES + ".{tissue}.rnaseq.raw.vcf.gz.tbi"
    threads: 4
    params:
        gatk_module = GATK_MODULE,
        gatk_conda_module = GATK_CONDA_MODULE,
        gatk_conda_env = GATK_CONDA_ENV,
        gatk_cmd = GATK_CMD,
        java_mem = "32g"
    shell:
        r"""
        module purge
        if [ -n "{params.gatk_module}" ]; then module load {params.gatk_module}; fi
        if [ -n "{params.gatk_conda_module}" ]; then module load {params.gatk_conda_module}; fi
        if [ -n "{params.gatk_conda_env}" ]; then source activate {params.gatk_conda_env}; fi

        mkdir -p $(dirname {output.vcf})

        {params.gatk_cmd} --java-options "-Xmx{params.java_mem}" GenotypeGVCFs \
          -R {input.ref} \
          -V {input.gvcf} \
          -O {output.vcf}

        {params.gatk_cmd} IndexFeatureFile -I {output.vcf}
        """


rule rnaseq_filter_pass_snps_by_tissue:
    input:
        ref = RNASEQ_VARIANT_REF,
        dict = RNASEQ_VARIANT_REF_DICT,
        fai = RNASEQ_VARIANT_REF + ".fai",
        vcf = RNASEQ_VARIANT_DIR + "/{tissue}/05-raw/" + SPECIES + ".{tissue}.rnaseq.raw.vcf.gz",
        tbi = RNASEQ_VARIANT_DIR + "/{tissue}/05-raw/" + SPECIES + ".{tissue}.rnaseq.raw.vcf.gz.tbi"
    output:
        snps = RNASEQ_VARIANT_DIR + "/{tissue}/06-filtered/" + SPECIES + ".{tissue}.rnaseq.raw.snps.vcf.gz",
        snps_tbi = RNASEQ_VARIANT_DIR + "/{tissue}/06-filtered/" + SPECIES + ".{tissue}.rnaseq.raw.snps.vcf.gz.tbi",
        filtered = RNASEQ_VARIANT_DIR + "/{tissue}/06-filtered/" + SPECIES + ".{tissue}.rnaseq.filtered.snps.vcf.gz",
        filtered_tbi = RNASEQ_VARIANT_DIR + "/{tissue}/06-filtered/" + SPECIES + ".{tissue}.rnaseq.filtered.snps.vcf.gz.tbi",
        pass_vcf = RNASEQ_VARIANT_DIR + "/{tissue}/06-filtered/" + SPECIES + ".{tissue}.rnaseq.pass.snps.vcf.gz",
        pass_tbi = RNASEQ_VARIANT_DIR + "/{tissue}/06-filtered/" + SPECIES + ".{tissue}.rnaseq.pass.snps.vcf.gz.tbi"
    threads: 4
    params:
        gatk_module = GATK_MODULE,
        gatk_conda_module = GATK_CONDA_MODULE,
        gatk_conda_env = GATK_CONDA_ENV,
        gatk_cmd = GATK_CMD,
        java_mem = "24g"
    shell:
        r"""
        module purge
        if [ -n "{params.gatk_module}" ]; then module load {params.gatk_module}; fi
        if [ -n "{params.gatk_conda_module}" ]; then module load {params.gatk_conda_module}; fi
        if [ -n "{params.gatk_conda_env}" ]; then source activate {params.gatk_conda_env}; fi

        mkdir -p $(dirname {output.pass_vcf})

        # Keep SNPs first. RNA-seq indels are more error-prone around splice junctions.
        {params.gatk_cmd} --java-options "-Xmx{params.java_mem}" SelectVariants \
          -R {input.ref} \
          -V {input.vcf} \
          --select-type-to-include SNP \
          -O {output.snps}
        {params.gatk_cmd} IndexFeatureFile -I {output.snps}

        # No known truth set is available for VQSR/BQSR, so use transparent hard filters.
        {params.gatk_cmd} --java-options "-Xmx{params.java_mem}" VariantFiltration \
          -R {input.ref} \
          -V {output.snps} \
          -O {output.filtered} \
          --filter-name "RNAseq_QD2" --filter-expression "QD < 2.0" \
          --filter-name "RNAseq_FS30" --filter-expression "FS > 30.0" \
          --filter-name "RNAseq_SOR3" --filter-expression "SOR > 3.0"
        {params.gatk_cmd} IndexFeatureFile -I {output.filtered}

        # Final file for downstream genetic variation analyses: SNPs that pass filters.
        {params.gatk_cmd} --java-options "-Xmx{params.java_mem}" SelectVariants \
          -R {input.ref} \
          -V {output.filtered} \
          --exclude-filtered true \
          -O {output.pass_vcf}
        {params.gatk_cmd} IndexFeatureFile -I {output.pass_vcf}
        """
