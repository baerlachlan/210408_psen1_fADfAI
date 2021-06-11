## Define a known set of variants based on the data, opposed to using a pre-existing database
## First we detect all potenital SNVs
rule detectSnvs:
	input:
		bam = "06_addRG/bam/{SAMPLE}.bam",
		bamIndex = "06_addRG/bam/{SAMPLE}.bai",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		vcf = temp("07_knownSnvs/detected/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("07_knownSnvs/detected/{SAMPLE}.vcf.gz.tbi")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
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

## Here we filter the detected SNVs with stringent criteria
## These SNVs will be used for measuring base recalibration
rule knownSnvs:
	input:
		vcf = "07_knownSnvs/detected/{SAMPLE}.vcf.gz",
		vcfIndex = "07_knownSnvs/detected/{SAMPLE}.vcf.gz.tbi",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		vcf = temp("07_knownSnvs/filtered/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("07_knownSnvs/filtered/{SAMPLE}.vcf.gz.tbi")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 1000,
		time = "00-02:00:00"
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