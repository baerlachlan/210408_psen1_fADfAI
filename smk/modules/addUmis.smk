rule addUmis:
	input:
		R1 = "00_rawData/fastq/{SAMPLE}_R1.fastq.gz",
		R2 = "00_rawData/fastq/{SAMPLE}_R2.fastq.gz",
		UMI = "00_rawData/fastq/{SAMPLE}_I1.fastq.gz"
	output:
		R1 = temp("01_addUmi/fastq/{SAMPLE}_R1.fastq.gz"),
		R2 = temp("01_addUmi/fastq/{SAMPLE}_R2.fastq.gz"),
		UMI1 = temp("01_addUmi/fastq/{SAMPLE}_I1.fastq.gz"),
		UMI2 = temp("01_addUmi/fastq/{SAMPLE}_I2.fastq.gz"),
		html = "01_addUmi/log/{SAMPLE}.html"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 4000,
		time = "00-01:00:00"
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