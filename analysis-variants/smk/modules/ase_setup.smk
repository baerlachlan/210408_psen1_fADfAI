rule ase_setup:
    input:
        bams = expand(os.path.join(analysis.bqsr_dir, "bam", "{SAMPLE}.bam"), SAMPLE=analysis.samples),
    output:
        snvDir = temp(directory(expand(os.path.join(analysis.wasp_dir, "1_snvs", "{SAMPLE}"), SAMPLE=analysis.samples))),
        intervals = temp(expand(os.path.join(analysis.aseRC_dir, "intervals", "{SAMPLE}.intervals"), SAMPLE=analysis.samples)),
    conda:
        "../envs/ase.yml"
    params:
        proj_root = os.getcwd(),
        variants_dir = analysis.variants_dir,
        wasp_dir = analysis.wasp_dir,
        aseRC_dir = analysis.aseRC_dir,
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-02:00:00",
    script:
        "../scripts/ase_setup.R"

