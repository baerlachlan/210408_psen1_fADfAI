## When MarkDuplicates (Picard) is run on coordinate sorted BAM files, unmapped mates of mapped records and supplementary/secondary alignments are excluded from the duplication test
## For variant analysis with GATK this is not a problem because HaplotypeCaller filters unmapped reads and secondary alignments before analysing
##
## Here we also merge lanes in the same step
rule markDuplicates:
	input:
		bam_lane1 = "05_groupUmis/bam/{SAMPLE}_L001.bam",
		bamIndex_lane1 = "05_groupUmis/bam/{SAMPLE}_L001.bam.bai",
		bam_lane2 = "05_groupUmis/bam/{SAMPLE}_L002.bam",
		bamIndex_lane2 = "05_groupUmis/bam/{SAMPLE}_L002.bam.bai"
	output:
		bam = temp("06_markDuplicates/bam/{SAMPLE}.bam"),
		bamIndex = temp("06_markDuplicates/bam/{SAMPLE}.bai"),
		metrics = "06_markDuplicates/metrics/{SAMPLE}.tsv",
		samstats = "06_markDuplicates/samstats/{SAMPLE}.tsv"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 32000,
		time = "00-04:00:00"
	shell:
		"""
		gatk \
			MarkDuplicates \
			--INPUT {input.bam_lane1} \
			--INPUT {input.bam_lane2} \
			--OUTPUT {output.bam}  \
			--BARCODE_TAG BX \
			--DUPLICATE_SCORING_STRATEGY RANDOM \
			--CREATE_INDEX true \
			--METRICS_FILE {output.metrics}

		samtools stats {output.bam} > {output.samstats}
		"""