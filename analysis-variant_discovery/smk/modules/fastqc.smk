rule fastqc_raw:
    input:
        path.join(analysis.raw_dir, "fastq", "{SAMPLE}" + analysis.fastq_ext),
    output:
        multiext(path.join(analysis.raw_dir, "FastQC", "{SAMPLE}_fastqc"), ".zip", ".html"),
    params:
        outDir = path.join(analysis.raw_dir, "FastQC"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-01:00:00",
    shell:
        "fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"

rule fastqc_trim:
    input:
        path.join(analysis.trim_dir, "fastq", "{SAMPLE}" + analysis.fastq_ext),
    output:
        multiext(path.join(analysis.trim_dir, "FastQC", "{SAMPLE}_fastqc"), ".zip", ".html"),
    params:
        outDir = path.join(analysis.trim_dir, "FastQC"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 2000,
        time = "00-01:00:00",
    shell:
        "fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"

rule fastqc_align:
    input:
        rules.align.output.bam,
    output:
        multiext(path.join(analysis.align_dir, "FastQC", "{SAMPLE}_fastqc"), ".zip", ".html"),
    params:
        outDir = path.join(analysis.align_dir, "FastQC"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-01:00:00",
    shell:
        "fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"