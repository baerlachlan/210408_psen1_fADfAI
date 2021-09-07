#!/bin/bash

## This script follows the recommended gatk workflow for RNA-seq short variant discovery (SNPs + Indels):
## - https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-

## Minor adjustments have been made as follows:
## - The workflow assumes raw data is in the unmapped BAM format
##   - As the raw data for this dataset is in FASTQ format, the initial SamToFastq method was skipped
##   - Later in the workflow, the MergeBamAlignment step was also skipped
##     - Unmapped BAM files contain metadata and this step allows for merging of that metadata into the aligned BAMS
##     - This was not needed as raw data was FASTQ which does not contain metadata
## - A bootstrapping approach was taken to define a set of known variants
##   - Our zebrafish are most likely quite different to the reference models due to years of breeding
##   - It is better to therefore generate a set of known variants from our data than to use pre-defined ones in a database (dbSNP , Ensembl)
##   - The reasoning and process for this is located at https://gatk.broadinstitute.org/hc/en-us/articles/360035890531?id=44 (Section 3)

SAMPLES = [
	"1_KB_A10", "2_KB_A11", "3_KB_A1", "4_KB_A6",
	"5_KB_A7", "6_KB_B11", "7_KB_B12", "8_KB_B3",
	"9_KB_B5", "10_KB_B6", "11_KB_B7", "12_KB_B8",
	"13_KB_B9", "14_KB_C10", "15_KB_C11", "16_KB_C12",
	"17_KB_C1", "18_KB_C2", "19_KB_C3", "20_KB_C4",
	"21_KB_C5", "22_KB_C6", "23_KB_C7", "24_KB_C8"
]
PAIR_ID = ["R1", "R2"]
REF_EXT = ["dict", "fa.fai"]
FQC_EXT = ["zip", "html"]
VCF_EXT = ["vcf.gz", "vcf.gz.tbi"]

## Set variables that may change between datasets
REFDIR = "/hpcfs/users/a1647910/refs/ensembl-release-101/danio_rerio/"
READ_LEN = 100

rule all:
	input:
		expand("00_rawData/FastQC/{SAMPLE}_{PAIR}_fastqc.{EXT}", SAMPLE = SAMPLES, PAIR = PAIR_ID, EXT = FQC_EXT),
		expand("02_trim/FastQC/{SAMPLE}_{PAIR}_fastqc.{EXT}", SAMPLE = SAMPLES, PAIR = PAIR_ID, EXT = FQC_EXT),
		expand("03_align/FastQC/{SAMPLE}_fastqc.{EXT}", SAMPLE = SAMPLES, PAIR = PAIR_ID, EXT = FQC_EXT),
		"03_align/featureCounts/genes.out",
		expand("08_dbsnp/4_selected/{SAMPLE}_snvs.vcf.gz", SAMPLE = SAMPLES),
		# expand("09_recalBases/recal/{SAMPLE}.analyzeCovariates.csv", SAMPLE = SAMPLES),
		# expand("10_callSnvs/4_selected/{SAMPLE}.vcf.gz", SAMPLE = SAMPLES),
		# expand("12_aseReadCounter/{DIR}/{SAMPLE}.tsv", DIR = ["wasp", "nowasp"], SAMPLE = SAMPLES),
		# expand("13_geneiase/2_ase/{SAMPLE}.static.pval.tsv", SAMPLE = SAMPLES)

include: "smk/modules/refs.smk"
include: "smk/modules/fastqc_raw.smk"
include: "smk/modules/addUmis.smk"
include: "smk/modules/trim.smk"
include: "smk/modules/align.smk"
include: "smk/modules/featureCounts.smk"
include: "smk/modules/groupUmis.smk"
include: "smk/modules/markDuplicates.smk"
include: "smk/modules/splitNCigar.smk"
include: "smk/modules/addRG.smk"
include: "smk/modules/dbsnp.smk"
include: "smk/modules/recalBases.smk"
include: "smk/modules/callSnvs.smk"
# include: "smk/modules/wasp.smk"
# include: "smk/modules/aseReadCounter.smk"
# include: "smk/modules/geneiase.smk"
