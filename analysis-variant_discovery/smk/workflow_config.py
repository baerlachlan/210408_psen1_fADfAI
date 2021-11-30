from os import path
import pandas as pd
from snakemake.io import expand

class RNASeqVariantDiscovery:

    def __init__(
            self, species, assembly, ensembl_release, read_length, fastq_ext,
            pair_tags, merge_samples = False, merge_tags = "",
            add_umi_header=False, umi_tag="", umi_dedup=False,
            bootstrap_known_variants=False, samples_tsv="smk/data/samples.tsv"):

        self.species = species
        self.assembly = assembly
        self.ensembl_release = ensembl_release
        self.read_length = read_length
        self.fastq_ext = fastq_ext
        self.pair_tags = pair_tags
        self.merge_samples = merge_samples
        self.merge_tags = merge_tags
        self.add_umi_header = add_umi_header
        self.umi_tag = umi_tag
        self.umi_dedup = umi_dedup
        self.bootstrap_known_variants = bootstrap_known_variants
        self.samples_tsv = samples_tsv

        ## Run methods
        self.dirs()
        self.read_samples_tsv()
        self.references()
        self.targets()

    ####
    ## Build directory structure
    ####

    def dirs(self):
        ## Define dir indexes
        raw_ind = 0
        addUmis_ind = 1
        trim_ind = 2
        align_ind = 3
        addRG_ind = 4
        groupUmis_ind = 5
        mergeSamples_ind = 6
        markDuplicates_ind = 7
        splitNCigar_ind = 8
        knownVariants_ind = 9
        bqsr_ind = 10
        variants_ind = 11
        ## Index decreases by 1 for any steps following optional steps that are not run
        if not self.add_umi_header:
            trim_ind -= 1
            align_ind -= 1
            addRG_ind -= 1
            mergeSamples_ind -= 1
            markDuplicates_ind -= 1
            splitNCigar_ind -= 1
            knownVariants_ind -= 1
            bqsr_ind -= 1
            variants_ind -= 1
        if not self.umi_dedup:
            mergeSamples_ind -= 1
            markDuplicates_ind -= 1
            splitNCigar_ind -= 1
            knownVariants_ind -= 1
            bqsr_ind -= 1
            variants_ind -= 1
        if not self.merge_samples:
            markDuplicates_ind -= 1
            splitNCigar_ind -= 1
            knownVariants_ind -= 1
            bqsr_ind -= 1
            variants_ind -= 1
        if not self.bootstrap_known_variants:
            bqsr_ind -= 1
            variants_ind -= 1
        ## Build dir names
        self.refs_dir = "refs"
        self.raw_dir = str(raw_ind).zfill(2) + "_rawData"
        self.trim_dir = str(trim_ind).zfill(2) + "_trim"
        self.align_dir = str(align_ind).zfill(2) + "_align"
        self.addRG_dir = str(addRG_ind).zfill(2) + "_addRG"
        self.markDuplicates_dir = str(markDuplicates_ind).zfill(2) + "_markDuplicates"
        self.splitNCigar_dir = str(splitNCigar_ind).zfill(2) + "_splitNCigar"
        self.bqsr_dir = str(bqsr_ind).zfill(2) + "_bqsr"
        self.variants_dir = str(variants_ind).zfill(2) + "_variants"
        if self.add_umi_header:
            self.addUmis_dir = str(addUmis_ind).zfill(2) + "_addUmis"
        if self.umi_dedup:
            self.groupUmis_dir = str(groupUmis_ind).zfill(2) + "_groupUmis"
        if self.merge_samples:
            self.mergeSamples_dir = str(mergeSamples_ind).zfill(2) + "_mergeSamples"
        if self.bootstrap_known_variants:
            self.knownVariants_dir = str(knownVariants_ind).zfill(2) + "_knownVariants"

    ## For debugging
    def print_dirs(self):
        print(self.raw_dir)
        if self.add_umi_header:
            print(self.addUmis_dir)
        print(self.trim_dir)
        print(self.align_dir)
        print(self.addRG_dir)
        if self.umi_dedup:
            print(self.groupUmis_dir)
        if self.merge_samples:
            print(self.mergeSamples_dir)
        print(self.markDuplicates_dir)
        print(self.splitNCigar_dir)
        if self.bootstrap_known_variants:
            print(self.knownVariants_dir)
        print(self.bqsr_dir)
        print(self.variants_dir)

    ####
    ## Samples
    ####

    def read_samples_tsv(self):
        # self.samples = os.listdir(path.join(self.raw_dir, "fastq"))
        # self.samples = [
        #     sample.replace(self.fastq_ext, "") for sample in self.samples
        # ]
        # for id in self.pair_tags:
        #     self.samples = [sample.replace(id, "") for sample in self.samples]
        # self.samples = [sample.replace(self.umi_tag, "") for sample in self.samples]
        # self.samples = list(set(self.samples))  ## Remove duplicates after removal of tags
        samples_df = pd.read_csv(self.samples_tsv, sep="\t")
        self.samples = list(set(samples_df["sample"]))
        self.RGID = dict(zip(samples_df["basename"], samples_df["read_group"]))
        self.RGSM = dict(zip(samples_df["basename"], samples_df["sample"]))

    ####
    ## Reference filenames, paths and urls for downloading
    ####

    def references(self):
        self.ref_fa = path.join(
            ".".join([
                self.species.capitalize(),
                self.assembly,
                "dna.primary_assembly.fa.gz"
            ])
        )
        self.ref_path = path.join(self.refs_dir, self.ref_fa.rstrip(".gz"))  ## File will be unzipped during download
        self.ref_url = path.join(
            "http://ftp.ensembl.org/pub",
            "release-" + str(self.ensembl_release),
            "fasta",
            self.species,
            "dna",
            self.ref_fa
        )
        self.gtf = path.join(
            ".".join([
                self.species.capitalize(),
                self.assembly,
                str(self.ensembl_release),
                "chr.gtf.gz"
            ])
        )
        self.gtf_path = path.join(self.refs_dir, self.gtf)
        self.gtf_url = path.join(
            "http://ftp.ensembl.org/pub",
            "release-" + str(self.ensembl_release),
            "gtf",
            self.species,
            self.gtf
        )
        self.knownVariants = path.join(".".join([self.species, "vcf.gz"]))
        self.knownVariants_path = path.join(self.refs_dir, self.knownVariants)
        self.knownVariants_url = path.join(
            "http://ftp.ensembl.org/pub",
            "release-" + str(self.ensembl_release),
            "variation/vcf",
            self.species,
            self.knownVariants
        )

    ####
    ## Input functions for rules succeeding optional rules
    ####

    def trim_inputs(self, wildcards):
        return {
            "R1": path.join(
                self.addUmis_dir if self.add_umi_header else self.raw_dir,
                "fastq",
                wildcards.SAMPLE + self.pair_tags[0] + self.fastq_ext
            ),
            "R2": path.join(
                self.addUmis_dir if self.add_umi_header else self.raw_dir,
                "fastq",
                wildcards.SAMPLE + self.pair_tags[1] + self.fastq_ext
            ),
        }

    def mergeSamples_inputs(self, wildcards):
        return {
            "bam" : expand(path.join(
                self.groupUmis_dir if self.umi_dedup else self.addRG_dir,
                "bam",
                wildcards.SAMPLE + "{MERGETAG}" + ".bam"
            ), MERGETAG = self.merge_tags),
            "bamIndex" : expand(path.join(
                self.groupUmis_dir if self.umi_dedup else self.addRG_dir,
                "bam",
                wildcards.SAMPLE + "{MERGETAG}" + ".bam.bai"
            ), MERGETAG = self.merge_tags)
        }

    def markDuplicates_inputs(self, wildcards):
        if self.merge_samples:
            input_dir = self.mergeSamples_dir
        elif self.umi_dedup:
            input_dir = self.groupUmis_dir
        else:
            input_dir = self.addRG_dir
        return {
            "bam" : path.join(
                input_dir,
                "bam",
                wildcards.SAMPLE + ".bam"
            ),
            "bamIndex" : path.join(
                input_dir,
                "bam",
                wildcards.SAMPLE + ".bam.bai"
            ),
        }

    def knownVariants_files(self, wildcards):
        if self.bootstrap_known_variants:
            return {
                "knownVariants" : path.join(
                    self.knownVariants_dir,
                    "6_select",
                    "known_variants.vcf.gz"
                ),
                "knownVariantsIndex" : path.join(
                    self.knownVariants_dir,
                    "6_select",
                    "known_variants.vcf.gz.tbi"
                )
            }
        else:
            return {
                "knownVariants" : self.knownVariants_path,
                "knownVariantsIndex" : self.knownVariants_path + ".tbi",
                }

    ####
    ## All file endpoints
    ####

    def targets(self):

        ## FastQC reports
        self.fqc = []
        fqc_raw = expand(
            path.join(self.raw_dir, "FastQC/{SAMPLE}{MERGETAG}{PAIRTAG}_fastqc.{EXT}"),
            SAMPLE=self.samples,
            MERGETAG=self.merge_tags,
            PAIRTAG=self.pair_tags,
            EXT=["html", "zip"]
        )
        self.fqc.extend(fqc_raw)
        fqc_trim = expand(
            path.join(self.trim_dir, "FastQC/{SAMPLE}{MERGETAG}{PAIRTAG}_fastqc.{EXT}"),
            SAMPLE=self.samples,
            MERGETAG=self.merge_tags,
            PAIRTAG=self.pair_tags,
            EXT=["html", "zip"]
        )
        self.fqc.extend(fqc_trim)
        fqc_align = expand(
            path.join(self.align_dir, "FastQC/{SAMPLE}{MERGETAG}_fastqc.{EXT}"),
            SAMPLE=self.samples,
            MERGETAG=self.merge_tags,
            EXT=["html", "zip"]
        )
        self.fqc.extend(fqc_align)

        ## BQSR summary for QC purposes
        self.bqsr = []
        # bqsr_firstPass = expand(
        #     path.join(self.bqsr_dir, "recal/{SAMPLE}.firstPass.table"),
        #     SAMPLE=self.samples
        # )
        # self.bqsr.extend(bqsr_firstPass)
        # bqsr_secondPass = expand(
        #     path.join(self.bqsr_dir, "recal/{SAMPLE}.secondPass.table"),
        #     SAMPLE=self.samples
        # )
        # self.bqsr.extend(bqsr_secondPass)
        bqsr_analyzeCovariates = expand(
            path.join(self.bqsr_dir, "recal/{SAMPLE}.analyzeCovariates.csv"),
            SAMPLE=self.samples
        )
        self.bqsr.extend(bqsr_analyzeCovariates)

        ## End of analysis variants
        self.variants = [path.join(self.variants_dir, "6_select/all_samples.vcf.gz")]

        ## All outputs
        self.outputs = []
        self.outputs.extend(self.fqc)
        self.outputs.extend(self.bqsr)
        self.outputs.extend(self.variants)