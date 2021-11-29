## When MarkDuplicates (Picard) is run on coordinate sorted BAM files, unmapped mates of mapped records and supplementary/secondary alignments are excluded from the duplication test
## For variant analysis with GATK this is not a problem because HaplotypeCaller filters unmapped reads and secondary alignments before analysing
rule markDuplicates:
    input:
        unpack(analysis.markDuplicates_inputs),
    output:
        bam = temp(path.join(analysis.markDuplicates_dir, "bam", "{SAMPLE}.bam")),
        bamIndex = temp(path.join(analysis.markDuplicates_dir, "bam", "{SAMPLE}.bai")),
        metrics = path.join(analysis.markDuplicates_dir, "metrics", "{SAMPLE}.tsv"),
        samstats = path.join(analysis.markDuplicates_dir, "samstats", "{SAMPLE}.tsv"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 32000,
        time = "00-04:00:00",
    shell:
        """
        gatk \
            MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam}  \
            --DUPLICATE_SCORING_STRATEGY RANDOM \
            --CREATE_INDEX true \
            --METRICS_FILE {output.metrics}

        samtools stats {output.bam} > {output.samstats}
        """