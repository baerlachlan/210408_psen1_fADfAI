rule trim:
	input:
		R1 = "01_addUmiHeader/fastq/{SAMPLE}_R1.fastq.gz",
		R2 = "01_addUmiHeader/fastq/{SAMPLE}_R2.fastq.gz"
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
			--unqualified_percent_limit 100 \
			--length_required 35 \
			--complexity_threshold 100 \
			--cut_front \
			--cut_tail \
			--trim_poly_g \
			--thread 1 \
			--html {output.html} \
			--json /dev/null \
		"""