rule count_pe:
    input:
        bam = expand(rules.align_pe.output.bam, SAMPLE = samples),
        gtf = rules.refs_downloadGtf.output,
    output:
        counts = os.path.join("results", counts_dir, "counts.out")
    params:
        Q = 1,
        s = 0,
        minOverlap = 35,
        fracOverlap = 0.9
    conda:
        "../envs/count.yml"
    resources:
        cpu = 4,
        ntasks = 1,
        mem_mb = 32000,
        time = "00-02:00:00",
    shell:
        """
        featureCounts \
            -Q {params.Q} \
            -s {params.s} \
            --minOverlap {params.minOverlap} \
            --fracOverlap {params.fracOverlap} \
            -T {resources.cpu} \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bam}
        """