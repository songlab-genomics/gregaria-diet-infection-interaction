### =================================================================
### STAR (index + alignment)
### =================================================================

rule STAR_index:
    input:
        ref_genome = REF_GENOME,
        annotation = REF_GTF
    output:
        index_dir = directory(STAR_INDEX_DIR)
    params:
        sjdb_overhang = STAR_SJDB_OVERHANG,
        align_intron_max = STAR_ALIGN_INTRON_MAX
    threads: 20
    shell:
        r"""
        module purge
        module load star-2.7.10b-gcc-12.1.0

        mkdir -p {output.index_dir}

        STAR --runMode genomeGenerate \
          --runThreadN {threads} \
          --genomeDir {output.index_dir} \
          --genomeFastaFiles {input.ref_genome} \
          --sjdbGTFfile {input.annotation} \
          --alignIntronMax {params.align_intron_max} \
          --sjdbOverhang {params.sjdb_overhang}
        """


rule STAR_align:
    input:
        index = STAR_INDEX_DIR,
        annotation = REF_GTF,
        trimmed_read1 = WORKDIR + "/01-{tissue}-trimmed-fastp/{locust}_1.trimmed.fastq.gz",
        trimmed_read2 = WORKDIR + "/01-{tissue}-trimmed-fastp/{locust}_2.trimmed.fastq.gz"
    output:
        bam = WORKDIR + "/02-{tissue}-star/{locust}_Aligned.sortedByCoord.out.bam",
        bai = WORKDIR + "/02-{tissue}-star/{locust}_Aligned.sortedByCoord.out.bam.csi",
        readtable = WORKDIR + "/02-{tissue}-star/{locust}_ReadsPerGene.out.tab",
        counts = WORKDIR + "/03-{tissue}-DESeq2/{locust}_counts.txt"
    params:
        prefix = WORKDIR + "/02-{tissue}-star/{locust}_",
        mapq_unique = STAR_OUTSAM_MAPQ_UNIQUE,
        align_intron_max = STAR_ALIGN_INTRON_MAX,
        align_mates_gap_max = STAR_ALIGN_MATES_GAP_MAX
    threads: 16
    shell:
        r"""
        module purge
        module load star-2.7.10b-gcc-12.1.0

        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {output.counts})

        STAR --runThreadN {threads} \
          --genomeDir {input.index} \
          --genomeLoad NoSharedMemory \
          --limitBAMsortRAM 32000000000 \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMmapqUnique {params.mapq_unique} \
          --outSAMattrRGline ID:{wildcards.locust} SM:{wildcards.locust} LB:Stranded_Total_RNA_RiboZero PL:Illumina PU:NovaSeq6000 \
          --quantMode TranscriptomeSAM GeneCounts \
          --twopassMode Basic \
          --sjdbGTFfile {input.annotation} \
          --outSAMattributes NH HI AS NM MD \
          --alignIntronMax {params.align_intron_max} \
          --alignMatesGapMax {params.align_mates_gap_max} \
          --outSAMunmapped Within \
          --readFilesCommand zcat \
          --readFilesIn {input.trimmed_read1} {input.trimmed_read2} \
          --outFileNamePrefix {params.prefix}

        module purge
        module load samtools-1.21-gcc-12.1.0

        samtools index -c {output.bam}

        # STAR ReadsPerGene: keep gene + unstranded (col 2), first-strand (col 3), second-strand (col 4)
        # You were cutting 1,4 previously; keeping identical behavior:
        cut -f1,4 {output.readtable} | grep -v "_" > {output.counts}
        """
