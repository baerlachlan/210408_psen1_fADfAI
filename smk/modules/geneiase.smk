rule geneiase:
	input:
		"11_geniase/input/{SAMPLE}.static.tsv"
	output:
		"11_geniase/output/{SAMPLE}.static.pval.tsv"
	conda:
		"../envs/geneiase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "01-00:00:00"
	shell:
		"""
		../packages/geneiase-1.0.1/bin/geneiase \
			-t static \
			-i {input} \
			-o {output}
		"""