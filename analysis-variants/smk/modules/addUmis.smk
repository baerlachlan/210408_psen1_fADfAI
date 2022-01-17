rule addUmis:
	input:
		R1 = os.path.join(analysis.raw_dir, "fastq", "{SAMPLE}" + analysis.pair_tags[0] + analysis.fastq_ext),
		R2 = os.path.join(analysis.raw_dir, "fastq", "{SAMPLE}" + analysis.pair_tags[1] + analysis.fastq_ext),
		UMI = os.path.join(analysis.raw_dir, "fastq", "{SAMPLE}" + analysis.umi_tag + analysis.fastq_ext),
	output:
		R1 = temp(os.path.join(analysis.addUmis_dir, "fastq", "{SAMPLE}" + analysis.pair_tags[0] + analysis.fastq_ext)),
		R2 = temp(os.path.join(analysis.addUmis_dir, "fastq", "{SAMPLE}" + analysis.pair_tags[1] + analysis.fastq_ext)),
		UMI1 = temp(os.path.join(analysis.addUmis_dir, "fastq", "{SAMPLE}_I1" + analysis.fastq_ext)),
		UMI2 = temp(os.path.join(analysis.addUmis_dir, "fastq", "{SAMPLE}_I2" + analysis.fastq_ext)),
		html = os.path.join(analysis.addUmis_dir, "log", "{SAMPLE}.html"),
	conda:
		"../envs/gatk.yml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 4000,
		time = "00-01:00:00",
	shell:
		"""
		fastp \
        	-i {input.R1}  \
        	-I {input.UMI}  \
        	-o {output.R1} \
        	-O {output.UMI1} \
        	--html {output.html} \
			--json /dev/null \
        	--umi --umi_loc=read2 --umi_len=8 \
        	-G -Q -A -L -w 1 -u 100 -n 8 -Y 100

   		fastp \
        	-i {input.R2}  \
        	-I {input.UMI}  \
        	-o {output.R2} \
        	-O {output.UMI2} \
			--html /dev/null \
			--json /dev/null \
        	--umi --umi_loc=read2 --umi_len=8 \
        	-G -Q -A -L -w 1 -u 100 -n 8 -Y 100
		"""