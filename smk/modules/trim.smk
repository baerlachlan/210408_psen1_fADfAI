rule trim:
	input:
		R1 = "01_addUmi/fastq/{SAMPLE}_R1.fastq.gz",
		R2 = "01_addUmi/fastq/{SAMPLE}_R2.fastq.gz"
	output:
		R1 = temp("02_trim/fastq/{SAMPLE}_R1.fastq.gz"),
		R2 = temp("02_trim/fastq/{SAMPLE}_R2.fastq.gz"),
		html = "02_trim/log/{SAMPLE}.html"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 2000,
		time = "00-02:00:00"
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

rule fastqc_trim:
	input:
		"02_trim/fastq/{SAMPLE}.fastq.gz"
	output:
		"02_trim/FastQC/{SAMPLE}_fastqc.zip",
		"02_trim/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "02_trim/FastQC/"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-01:00:00"
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"