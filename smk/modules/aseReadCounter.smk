rule aseRC:
	input:
		bam = "11_wasp/5_merge/{SAMPLE}.keep.merge.sort.bam",
		bamIndex = "11_wasp/5_merge/{SAMPLE}.keep.merge.sort.bam.bai",
		vcf = "10_callSnvs/4_selected/{SAMPLE}.vcf.gz",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		intervals = "10_callSnvs/intervals/exons.intervals"
	output:
		tsv = "12_aseReadCounter/wasp/{SAMPLE}.tsv"
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
			-L {input.intervals} \
			-O {output.tsv} \
			--min-mapping-quality 10 \
			--min-base-quality 20
		"""

rule aseRC_nowasp:
	input:
		bam = "09_recalBases/bam/{SAMPLE}.bam",
		bamIndex = "09_recalBases/bam/{SAMPLE}.bai",
		vcf = "10_callSnvs/4_selected/{SAMPLE}.vcf.gz",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		intervals = "10_callSnvs/intervals/exons.intervals"
	output:
		tsv = "12_aseReadCounter/nowasp/{SAMPLE}.tsv"
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
			-L {input.intervals} \
			-O {output.tsv} \
			--min-mapping-quality 10 \
			--min-base-quality 20
		"""