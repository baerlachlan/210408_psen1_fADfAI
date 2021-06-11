rule fastqc_trim:
	input:
		"02_trimmedData/fastq/{SAMPLE}.fastq.gz"
	output:
		"02_trimmedData/FastQC/{SAMPLE}_fastqc.zip",
		"02_trimmedData/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "02_trimmedData/FastQC/"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-01:00:00"
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"