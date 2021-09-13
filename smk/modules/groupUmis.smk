## This step adds BAM tags that define UMI groups which are essentially duplicates
## The resulting BAM files are read for Picard MarkDuplicates
## Picard MarkDuplicates is preferred over Umi-tools as it allows for random selection of the representative read avoiding mapping bias
rule groupUmis:
	input:
		bam = "04_addRG/bam/{SAMPLE}_{LANE}.bam",
		bamIndex = "04_addRG/bam/{SAMPLE}_{LANE}.bai"
	output:
		bam = "05_groupUmis/bam/{SAMPLE}_{LANE}.bam",
		bamIndex = "05_groupUmis/bam/{SAMPLE}_{LANE}.bam.bai"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 8000,
		time = "00-05:00:00"
	shell:
		"""
		umi_tools group \
			-I {input.bam} \
			-S {output.bam} \
			--temp-dir=. \
			--output-bam \
			--method=unique \
			--extract-umi-method=read_id \
			--umi-separator=":" \
			--paired

		samtools index {output.bam}
		"""