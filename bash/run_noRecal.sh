#!/bin/bash

## Directories
PROJROOT=/hpcfs/users/a1647910/210408_psen1_fADfAI_snv
ADDRG=${PROJROOT}/06_addRG

## Samples
BAMS=$(ls ${ADDRG}/bam/*.bam)
echo -e "Found files:\n${BAMS}"

## Send jobs with sbatch
for BAM in ${BAMS}
do
    sbatch ${PROJROOT}/bash/callVariants_noRecal.sh ${BAM}
done