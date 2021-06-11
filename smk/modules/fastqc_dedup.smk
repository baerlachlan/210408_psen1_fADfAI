rule fastqc_dedup:
	input:
		"04_dedupUmis/bam/{SAMPLE}.bam"
	output:
		"04_dedupUmis/FastQC/{SAMPLE}_fastqc.zip",
		"04_dedupUmis/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "04_dedupUmis/FastQC/"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-02:00:00"
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"