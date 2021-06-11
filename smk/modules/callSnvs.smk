rule callSnvs:
	input:
		bam = "08_recalBases/bam/{SAMPLE}.bam",
		bamIndex = "08_recalBases/bam/{SAMPLE}.bai",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		vcf = temp("09_callSnvs/called/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("09_callSnvs/called/{SAMPLE}.vcf.gz.tbi")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 8000,
		time = "03-00:00:00"
	shell:
		"""
		gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
			HaplotypeCaller \
			-R {input.refFa} \
			-I {input.bam} \
			-O {output.vcf} \
			-dont-use-soft-clipped-bases \
			--standard-min-confidence-threshold-for-calling 20
		"""

rule filterSnvs:
	input:
		vcf = "09_callSnvs/called/{SAMPLE}.vcf.gz",
		vcfIndex = "09_callSnvs/called/{SAMPLE}.vcf.gz.tbi",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		vcf = temp("09_callSnvs/filtered/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("09_callSnvs/filtered/{SAMPLE}.vcf.gz.tbi")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 1000,
		time = "00-00:30:00"
	shell:
		"""
		gatk \
		    VariantFiltration \
			--R {input.refFa} \
			--V {input.vcf} \
			--window 35 \
			--cluster 3 \
			--filter-name "FS" \
			--filter "FS > 30.0" \
			--filter-name "QD" \
			--filter "QD < 2.0" \
			-O {output.vcf}
		"""

rule selectSnvs:
	input:
		vcf = "09_callSnvs/filtered/{SAMPLE}.vcf.gz",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa"
	output:
		vcf = "09_callSnvs/selected/{SAMPLE}.vcf.gz",
		vcfIndex = "09_callSnvs/selected/{SAMPLE}.vcf.gz.tbi"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-00:30:00"
	shell:
		"""
		gatk \
			SelectVariants \
			-R {input.refFa} \
			-V {input.vcf} \
			--select-type-to-include SNP \
			-O {output.vcf}
		"""