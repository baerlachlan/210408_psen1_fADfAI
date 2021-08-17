rule featureCounts:
	input:
		bam = expand("03_align/bam/{SAMPLE}.bam", SAMPLE = SAMPLES),
		gtf = REFDIR + "Danio_rerio.GRCz11.101.chr.gtf.gz"
	output:
		counts = "03_align/featureCounts/counts.out",
		summary = "03_align/featureCounts/counts.out.summary",
		genes = "03_align/featureCounts/genes.out"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 2,
		ntasks = 2,
		mem_mb = 4000,
		time = "00-02:00:00"
	shell:
		"""
		featureCounts \
			-Q 10 \
			-s 0 \
			-T 4 \
			-p \
			-a {input.gtf} \
			-o {output.counts} {input.bam}

		## Storing the output in a single file
		cut -f1,7- {output.counts} | \
		sed 1d > {output.genes}
		"""