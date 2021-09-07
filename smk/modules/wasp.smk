rule findIntersecting:
	input:
		bam = "09_recalBases/bam/{SAMPLE}.bam",
		snpDir = "11_wasp/1_snvs/{SAMPLE}"
	output:
		fq1 = temp("11_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.remap.fq1.gz"),
		fq2 = temp("11_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.remap.fq2.gz"),
		single = temp("11_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.remap.single.fq.gz"),
		to_remap = temp("11_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.to.remap.bam"),
		keep_intersect = temp("11_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.keep.bam")
	params:
		outDir = temp(directory("11_wasp/2_findIntersecting/{SAMPLE}"))
	conda:
		"../envs/wasp.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-08:00:00"
	shell:
		"""
		python ../packages/WASP/mapping/find_intersecting_snps.py \
    	    --is_paired_end \
    	    --is_sorted \
			--output_dir {params.outDir} \
    	    --snp_dir {input.snpDir} \
    	    {input.bam}
		"""

rule remap:
	input:
		fq1 = "11_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.remap.fq1.gz",
		fq2 = "11_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.remap.fq2.gz",
		starIndex = "refs/star/"
	output:
		remapped_unsorted = temp("11_wasp/3_remap/{SAMPLE}Aligned.out.bam"),
		remapped_sorted = temp("11_wasp/3_remap/{SAMPLE}sorted.out.bam"),
		index = temp("11_wasp/3_remap/{SAMPLE}sorted.out.bam.bai"),
		STARgenome = temp(directory("11_wasp/3_remap/{SAMPLE}_STARgenome")),
		STARpass1 = temp(directory("11_wasp/3_remap/{SAMPLE}_STARpass1"))
	params:
		overhang = READ_LEN-1,
		bname = "11_wasp/3_remap/{SAMPLE}"
	conda:
		"../envs/wasp.yaml"
	resources:
		cpu = 16,
		ntasks = 1,
		mem_mb = 32000,
		time = "00-04:00:00"
	shell:
		"""
		STAR \
			--genomeDir {input.starIndex}\
			--runThreadN {resources.cpu} \
			--readFilesIn {input.fq1} {input.fq2} \
			--readFilesCommand "gunzip -c" \
			--sjdbOverhang {params.overhang} \
			--outSAMtype BAM Unsorted \
			--twopassMode Basic \
			--outFileNamePrefix {params.bname}

		mkdir -p 11_wasp/3_remap/log
		mv {params.bname}*out 11_wasp/3_remap/log
		mv {params.bname}*tab 11_wasp/3_remap/log

		samtools sort -o {output.remapped_sorted} {output.remapped_unsorted}
		samtools index {output.remapped_sorted}
		"""

rule filterRemapped:
	input:
		to_remap = "11_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.to.remap.bam",
		remapped_unsorted = "11_wasp/3_remap/{SAMPLE}Aligned.out.bam",
		remapped_sorted = "11_wasp/3_remap/{SAMPLE}sorted.out.bam"
	output:
		keep_filter = temp("11_wasp/4_filterRemapped/{SAMPLE}.keep.bam")
	conda:
		"../envs/wasp.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 24000,
		time = "00-02:00:00"
	shell:
		"""
		python ../packages/WASP/mapping/filter_remapped_reads.py \
			{input.to_remap} \
			{input.remapped_sorted} \
			{output.keep_filter}
		"""

rule merge:
	input:
		keep_filter = "11_wasp/4_filterRemapped/{SAMPLE}.keep.bam",
		keep_intersect = "11_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.keep.bam"
	output:
		keep_merged = temp("11_wasp/5_merge/{SAMPLE}.keep.merge.bam"),
		keep_sorted = temp("11_wasp/5_merge/{SAMPLE}.keep.merge.sort.bam"),
		keep_sortedIndex = temp("11_wasp/5_merge/{SAMPLE}.keep.merge.sort.bam.bai")
	conda:
		"../envs/wasp.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "00-02:00:00"
	shell:
		"""
		samtools merge {output.keep_merged} \
			{input.keep_filter} \
			{input.keep_intersect}
		samtools sort -o {output.keep_sorted} \
			{output.keep_merged}
		samtools index {output.keep_sorted}
		"""
