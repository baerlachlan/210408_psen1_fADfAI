rule intervals:
    output:
        os.path.join(analysis.refs_dir, "exons.intervals"),
    params:
        species = analysis.species,
        ensembl_release = analysis.ensembl_release,
    conda:
        "../envs/r.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-00:30:00",
    script:
        "../scripts/intervals.R"