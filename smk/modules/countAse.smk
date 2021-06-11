rule aseReadCounts:
	input:
		bam = "08_recalBases/bam/{SAMPLE}.bam",
		bamIndex = "08_recalBases/bam/{SAMPLE}.bai",
		vcf = "09_callSnvs/selected/{SAMPLE}.vcf.gz",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa"
	output:
		tsv = "11_countAse/{SAMPLE}.tsv"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "00-04:00:00"
	shell:
		"""
		gatk \
			ASEReadCounter \
			-I {input.bam} \
			-V {input.vcf} \
			-R {input.refFa} \
			-O {output.tsv} \
			--min-mapping-quality 10 \
			--min-base-quality 20
		"""