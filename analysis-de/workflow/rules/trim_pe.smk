rule trim_pe:
    input:
        unpack(trim_inputs_pe)
    output:
        R1 = temp(os.path.join(
            "results", trim_dir, "fastq",
            "{SAMPLE}" + config["paired_end"]["tags"][0] + config["fastq_ext"]
        )),
        R2 = temp(os.path.join(
            "results", trim_dir, "fastq",
            "{SAMPLE}" + config["paired_end"]["tags"][1] + config["fastq_ext"]
        )),
        html = os.path.join("results", trim_dir, "log", "{SAMPLE}.html")
    conda:
        "../envs/trim.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-02:00:00",
    shell:
        """
        fastp \
            -i {input.R1} \
            -I {input.R2} \
            -o {output.R1} \
            -O {output.R2} \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 35 \
            --trim_poly_g \
            --thread 1 \
            --html {output.html} \
            --json /dev/null \
        """