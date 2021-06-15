rule fastqc_align:
	input:
		"03_align/bam/{SAMPLE}Aligned.sortedByCoord.out.bam"
	output:
		"03_align/FastQC/{SAMPLE}Aligned.sortedByCoord.out_fastqc.zip",
		"03_align/FastQC/{SAMPLE}Aligned.sortedByCoord.out_fastqc.html"
	params:
		outDir = "03_align/FastQC/"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-01:00:00"
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"