rule findIntersecting:
	input:
		bam = "08_recalBases/bam/{SAMPLE}.bam"
		snpDir = directory("10_wasp/snvs/{SAMPLE}")
	output:
		outDir = temp(directory("10_wasp/findIntersecting"))
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 2,
		ntasks = 2,
		mem_mb = 8000,
		time = "00-10:00:00"
	shell:
		"""
		python ../packages/WASP/mapping/find_intersecting_snps.py \
    	    --is_paired_end \
    	    --is_sorted \
			--output_dir {output.outDir} \
    	    --snp_dir {input.snpDir} \
    	    {input.bam}
		"""

rule remap:

rule filterRemapped:

rule merge: