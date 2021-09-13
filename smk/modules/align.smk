rule align:
	input:
		R1 = "02_trim/fastq/{SAMPLE}_{LANE}_R1.fastq.gz",
		R2 = "02_trim/fastq/{SAMPLE}_{LANE}_R2.fastq.gz",
		starIndex = "refs/star/"
	output:
		bamRenamed = temp("03_align/bam/{SAMPLE}_{LANE}.bam"),
		bamIndex = temp("03_align/bam/{SAMPLE}_{LANE}.bam.bai"),
		STARgenome = temp(directory("03_align/bam/{SAMPLE}_{LANE}_STARgenome")),
		STARpass1 = temp(directory("03_align/bam/{SAMPLE}_{LANE}_STARpass1"))
	params:
		overhang = READ_LEN-1,
		bname = "03_align/bam/{SAMPLE}_{LANE}",
		bamBeforeRename = "03_align/bam/{SAMPLE}_{LANE}Aligned.sortedByCoord.out.bam"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 16,
		ntasks = 1,
		mem_mb = 32000,
		time = "00-02:00:00"
	shell:
		"""
		STAR \
			--genomeDir {input.starIndex}\
			--runThreadN {resources.cpu} \
			--readFilesIn {input.R1} {input.R2} \
			--readFilesCommand "gunzip -c" \
			--sjdbOverhang {params.overhang} \
			--outSAMtype BAM SortedByCoordinate \
			--twopassMode Basic \
			--outFileNamePrefix {params.bname}

		mv {params.bamBeforeRename} {output.bamRenamed}

		mkdir -p 03_align/log
		mv {params.bname}*out 03_align/log
		mv {params.bname}*tab 03_align/log

		samtools index {output.bamRenamed}
		"""

rule fastqc_align:
	input:
		"03_align/bam/{SAMPLE}.bam"
	output:
		"03_align/FastQC/{SAMPLE}_fastqc.zip",
		"03_align/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "03_align/FastQC/"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-01:00:00"
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"