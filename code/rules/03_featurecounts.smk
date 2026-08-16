### =================================================================
### FEATURECOUNTS
### =================================================================

rule featurecounts:
    input:
        bamfile = WORKDIR + "/02-{tissue}-star/{locust}_Aligned.sortedByCoord.out.bam",
        annotation = REF_GTF
    output:
        counts_txt = WORKDIR + "/04-{tissue}-featurecounts/{locust}_counts.txt",
        final_counts = WORKDIR + "/03-{tissue}-DESeq2/{locust}_featurecounts.txt"
    threads: 12
    shell:
        r"""
        module purge
        ml mamba
        source activate subread

        mkdir -p $(dirname {output.counts_txt})
        mkdir -p $(dirname {output.final_counts})

        featureCounts -p \
          --countReadPairs \
          -B \
          -C \
          -s {FEATURECOUNTS_STRAND} \
          -t transcript,exon \
          -g gene_id \
          --extraAttributes gene_name \
          --primary \
          -Q 10 \
          -a {input.annotation} \
          -R BAM {input.bamfile} \
          -T {threads} \
          -o {output.counts_txt}

        cut -f 1,8 {output.counts_txt} | tail -n +3 > {output.final_counts}
        """
