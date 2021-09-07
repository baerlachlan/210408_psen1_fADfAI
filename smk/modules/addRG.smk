## gatk HaplotypeCaller requires read group (RG) tags
## An explanation of this is found at https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
rule addRG:
	input:
		bam = "06_splitNCigar/bam/{SAMPLE}.bam",
		bamIndex = "06_splitNCigar/bam/{SAMPLE}.bai"
	output:
		bam = temp("07_addRG/bam/{SAMPLE}.bam"),
		bamIndex = temp("07_addRG/bam/{SAMPLE}.bai"),
		samstats = "07_addRG/samstats/{SAMPLE}.tsv"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 4000,
		time = "00-02:00:00"
	shell:
		"""
		gatk \
			AddOrReplaceReadGroups \
    		-I {input.bam} \
   			-O {output.bam} \
    		-SORT_ORDER coordinate \
    		-RGID default \
    		-RGLB default \
    		-RGSM {wildcards.SAMPLE} \
			-RGPU default \
    		-RGPL ILLUMINA \
    		-CREATE_INDEX True

		samtools stats -d {output.bam} > {output.samstats}
		"""