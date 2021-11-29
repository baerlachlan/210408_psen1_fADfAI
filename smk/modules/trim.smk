rule trim:
    input:
        unpack(analysis.trim_inputs),
    output:
        R1 = temp(path.join(analysis.trim_dir, "fastq", "{SAMPLE}" + analysis.pair_tags[0] + analysis.fastq_ext)),
        R2 = temp(path.join(analysis.trim_dir, "fastq", "{SAMPLE}" + analysis.pair_tags[1] + analysis.fastq_ext)),
        html = path.join(analysis.trim_dir, "log", "{SAMPLE}.html"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 2000,
        time = "00-02:00:00",
    shell:
        """
        fastp \
            -i {input.R1}  \
            -I {input.R2}  \
            -o {output.R1} \
            -O {output.R2} \
            --qualified_quality_phred 20 \
            --length_required 35 \
            --trim_poly_g \
            --thread 1 \
            --html {output.html} \
            --json /dev/null \
        """