## gatk HaplotypeCaller requires read group (RG) tags
## An explanation of this is found at https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
rule addRG:
	input:
		bam = "03_align/bam/{SAMPLE}_{LANE}.bam",
		bamIndex = "03_align/bam/{SAMPLE}_{LANE}.bam.bai"
	output:
		bam = temp("04_addRG/bam/{SAMPLE}_{LANE}.bam"),
		bamIndex = temp("04_addRG/bam/{SAMPLE}_{LANE}.bai"),
		samstats = "04_addRG/samstats/{SAMPLE}_{LANE}.tsv"
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
    		-RGID {wildcards.SAMPLE}_{wildcards.LANE} \
			-RGPU null \
    		-RGSM {wildcards.SAMPLE} \
    		-RGPL ILLUMINA \
    		-RGLB null \
    		-CREATE_INDEX True

		samtools stats -d {output.bam} > {output.samstats}
		"""