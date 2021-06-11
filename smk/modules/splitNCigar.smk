rule splitNCigar:
	input:
		bam = "04_dedupUmis/bam/{SAMPLE}.bam",
		bamIndex = "04_dedupUmis/bam/{SAMPLE}.bam.bai",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		bam = temp("05_splitNCigar/bam/{SAMPLE}.bam"),
		bamIndex = temp("05_splitNCigar/bam/{SAMPLE}.bai")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 8,
		ntasks = 1,
		mem_mb = 32000,
		time = "00-04:00:00"
	shell:
		"""
		gatk \
            SplitNCigarReads \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output.bam}
		"""