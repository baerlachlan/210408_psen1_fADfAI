rule refs_downloadFa:
    output:
        temp(analysis.refFa_path),
    params:
        ref_url = analysis.ref_url,
        refGz_path = analysis.refGz_path
    threads:
        1
    shell:
        """
        # curl -L {params.ref_url} | gzip -d > {output}
        rsync -avP {params.ref_url} {params.refGz_path}
        gzip -d {params.refGz_path}
        """

rule refs_downloadGtf:
    output:
        temp(analysis.gtf_path),
    params:
        gtf_url = analysis.gtf_url,
    threads:
        1
    shell:
        """
        # curl -L {params.gtf_url} > {output}
        rsync -avP {params.gtf_url} {output}
        """

rule refs_downloadKnownVariants:
    output:
        temp(analysis.knownVariants_path),
    params:
        knownVariants_url = analysis.knownVariants_url,
    threads:
        1
    shell:
        """
        # curl -L {params.knownVariants_url} > {output}
        rsync -avP {params.knownVariants_url} {output}
        """

rule refs_refDict:
    input:
        rules.refs_downloadFa.output,
    output:
        temp(analysis.refFa_path.rstrip("fa") + "dict"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-00:10:00",
    shell:
        "gatk CreateSequenceDictionary -R {input}"

rule refs_refIndex:
    input:
        rules.refs_downloadFa.output,
    output:
        temp(analysis.refFa_path + ".fai"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-00:10:00",
    shell:
        "samtools faidx {input}"

rule refs_starIndex:
    input:
        ref_fa = rules.refs_downloadFa.output,
        gtf = rules.refs_downloadGtf.output,
    output:
        temp(directory(os.path.join(analysis.refs_dir, "star"))),
    params:
        overhang = analysis.read_length - 1,
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 16,
        ntasks = 1,
        mem_mb = 32000,
        time = "00-00:30:00",
    shell:
        """
        zcat {input.gtf} > temp.gtf

        STAR \
            --runThreadN {resources.cpu} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.ref_fa} \
            --sjdbGTFfile temp.gtf \
            --sjdbOverhang {params.overhang}

        rm temp.gtf
        """

rule refs_knownVariantsIndex:
    input:
        rules.refs_downloadKnownVariants.output,
    output:
        analysis.knownVariants_path + ".tbi",
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-00:10:00",
    shell:
        "tabix -p vcf {input}"