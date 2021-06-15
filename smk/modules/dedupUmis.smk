rule dedupUmis:
	input:
		bam = "03_align/bam/{SAMPLE}Aligned.sortedByCoord.out.bam",
		bamIndex = "03_align/bam/{SAMPLE}Aligned.sortedByCoord.out.bam.bai"
	output:
		bam = temp("04_dedupUmis/bam/{SAMPLE}.bam"),
		bamIndex = temp("04_dedupUmis/bam/{SAMPLE}.bam.bai"),
		log = "04_dedupUmis/log/{SAMPLE}.log"
	params:
		statsPrefix = "04_dedupUmis/log/{SAMPLE}"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 32000,
		time = "00-06:00:00"
	shell:
		"""
		umi_tools dedup \
    		-I {input.bam} \
    		-L {output.log} \
    		-S {output.bam} \
    		--umi-separator=":" \
    		--temp-dir=. \
    		--paired \
    		--output-stats={params.statsPrefix}

    	samtools index {output.bam}
		"""