#!/bin/bash -l
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-core=1
#SBATCH --time=72:00:00
#SBATCH --mem=8GB
#SBATCH -o /hpcfs/users/a1647910/210408_psen1_fADfAI_snv/slurm/%x_%j.out
#SBATCH -e /hpcfs/users/a1647910/210408_psen1_fADfAI_snv/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

## HPC modules
module load arch/haswell
module load arch/skylake

## Modules
module load Anaconda3/2020.07

## Activate conda env
conda activate SNV

## Reference
REFFA=/hpcfs/users/a1647910/refs/ensembl-release-101/danio_rerio/Danio_rerio.GRCz11.dna.primary_assembly.fa

##Directories
PROJROOT=/hpcfs/users/a1647910/210408_psen1_fADfAI_snv
ADDRG=${PROJROOT}/06_addRG
NORECAL=${PROJROOT}/07_callVariants_noRecal

## Create output directory
mkdir -p ${NORECAL}/vcf

## Define input
INBAM=$1
echo -e "input bam is ${INBAM}"

## Define output
OUTVCF=${NORECAL}/vcf/$(basename ${INBAM%.bam}.vcf.gz)
echo -e "output vcf will be ${OUTVCF}"

## Run HaplotypeCaller
gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
	HaplotypeCaller \
	-R ${REFFA} \
	-I ${INBAM} \
	-O ${OUTVCF} \
	-dont-use-soft-clipped-bases \
	--standard-min-confidence-threshold-for-calling 20

conda deactivate