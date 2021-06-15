rule align:
	input:
		R1 = "02_trim/fastq/{SAMPLE}_R1.fastq.gz",
		R2 = "02_trim/fastq/{SAMPLE}_R2.fastq.gz",
		starIndex = "refs/star/"
	output:
		bam = temp("03_align/bam/{SAMPLE}Aligned.sortedByCoord.out.bam"),
		bamIndex = temp("03_align/bam/{SAMPLE}Aligned.sortedByCoord.out.bam.bai")
	params:
		overhang = READ_LEN-1,
		bname = "03_align/bam/{SAMPLE}"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 16,
		ntasks = 1,
		mem_mb = 32000,
		time = "00-05:00:00"
	shell:
		"""
		STAR \
			--genomeDir {input.starIndex}\
			--runThreadN {resources.cpu} \
			--readFilesIn {input.R1} {input.R2} \
			--readFilesCommand "gunzip -c" \
			--sjdbOverhang {params.overhang} \
			--outSAMtype BAM SortedByCoordinate \
			--twopassMode Basic \
			--outFileNamePrefix {params.bname}


		mkdir -p 03_align/log
		mv {params.bname}*out 03_align/log
		mv {params.bname}*tab 03_align/log

		samtools index {output.bam}
		"""