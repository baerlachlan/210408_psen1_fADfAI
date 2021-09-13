rule recalFirstPass:
	input:
		bam = "07_splitNCigar/bam/{SAMPLE}.bam",
		bamIndex = "07_splitNCigar/bam/{SAMPLE}.bai",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict",
		snvs = "08_dbsnp/4_selected/{SAMPLE}_snvs.vcf.gz",
		indels = "08_dbsnp/4_selected/{SAMPLE}_indels.vcf.gz"
	output:
		temp("09_recalBases/recal/{SAMPLE}.firstPass.table")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 8000,
		time = "00-04:00:00"
	shell:
		"""
		gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output} \
            --known-sites {input.snvs} \
            --known-sites {input.indels}
		"""

rule applyRecal:
	input:
		bam = "07_splitNCigar/bam/{SAMPLE}.bam",
		bamIndex = "07_splitNCigar/bam/{SAMPLE}.bai",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict",
		recal = "09_recalBases/recal/{SAMPLE}.firstPass.table"
	output:
		bam = temp("09_recalBases/bam/{SAMPLE}.bam"),
		bamIndex = temp("09_recalBases/bam/{SAMPLE}.bai"),
		metrics = "09_recalBases/metrics/{SAMPLE}.tsv"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 8000,
		time = "00-04:00:00"
	shell:
		"""
		gatk \
            --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output.bam} \
            --bqsr-recal-file {input.recal}

		samtools stats -d {output.bam} > {output.metrics}
		"""

rule recalSecondPass:
	input:
		bam = "09_recalBases/bam/{SAMPLE}.bam",
		bamIndex = "09_recalBases/bam/{SAMPLE}.bai",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict",
		snvs = "08_dbsnp/4_selected/{SAMPLE}_snvs.vcf.gz",
		indels = "08_dbsnp/4_selected/{SAMPLE}_indels.vcf.gz"
	output:
		temp("09_recalBases/recal/{SAMPLE}.secondPass.table")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 8000,
		time = "00-04:00:00"
	shell:
		"""
		gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output} \
            --known-sites {input.snvs} \
            --known-sites {input.indels}
		"""

rule analyzeCovariates:
	input:
		firstPass = "09_recalBases/recal/{SAMPLE}.firstPass.table",
		secondPass = "09_recalBases/recal/{SAMPLE}.secondPass.table"
	output:
		csv = "09_recalBases/recal/{SAMPLE}.analyzeCovariates.csv"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "00-01:00:00"
	shell:
		"""
		gatk AnalyzeCovariates \
			-before {input.firstPass} \
     		-after {input.secondPass} \
     		-csv {output.csv}
		"""