rule addRG:
    input:
        bam = rules.align.output.bam,
        bamIndex = rules.align.output.bamIndex,
    output:
        bam = temp(path.join(analysis.addRG_dir, "bam", "{SAMPLE}.bam")),
        bamIndex = temp(path.join(analysis.addRG_dir, "bam", "{SAMPLE}.bam.bai")),
        samstats = path.join(analysis.addRG_dir, "samstats", "{SAMPLE}.tsv"),
    params:
        RGID = lambda wildcard: analysis.RGID[wildcard.SAMPLE],
        RGSM = lambda wildcard: analysis.RGSM[wildcard.SAMPLE],
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 4000,
        time = "00-02:00:00",
    shell:
        """
        gatk \
            AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.bam} \
            -SORT_ORDER coordinate \
            -RGID {params.RGID} \
            -RGPU null \
            -RGSM {params.RGSM} \
            -RGPL ILLUMINA \
            -RGLB null \
            -CREATE_INDEX False

        samtools index {output.bam}
        samtools stats -d {output.bam} > {output.samstats}
        """