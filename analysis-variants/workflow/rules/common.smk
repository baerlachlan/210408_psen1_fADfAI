####
## Manage config
####

if not config["umi"]["activate"]:
    config["umi"]["add_header"]["activate"] = False
if not config["umi"]["add_header"]["activate"]:
    config["umi"]["add_header"]["tag"] = []
if not config["merge_samples"]["activate"]:
    config["merge_samples"]["tags"] = []

####
## Directory structure
####

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
wasp_ind = 12
aseRC_ind = 13
## Index decreases by 1 for any steps following optional steps that are not run
if not config["umi"]["add_header"]["activate"]:
    trim_ind -= 1
    align_ind -= 1
    addRG_ind -= 1
    mergeSamples_ind -= 1
    groupUmis_ind -= 1
    markDuplicates_ind -= 1
    splitNCigar_ind -= 1
    knownVariants_ind -= 1
    bqsr_ind -= 1
    variants_ind -= 1
    wasp_ind -= 1
    aseRC_ind -= 1
if not config["umi"]["activate"]:
    mergeSamples_ind -= 1
    markDuplicates_ind -= 1
    splitNCigar_ind -= 1
    knownVariants_ind -= 1
    bqsr_ind -= 1
    variants_ind -= 1
    wasp_ind -= 1
    aseRC_ind -= 1
if not config["merge_samples"]["activate"]:
    markDuplicates_ind -= 1
    splitNCigar_ind -= 1
    knownVariants_ind -= 1
    bqsr_ind -= 1
    variants_ind -= 1
    wasp_ind -= 1
    aseRC_ind -= 1
if not config["bootstrap_known_variants"]["activate"]:
    bqsr_ind -= 1
    variants_ind -= 1
    wasp_ind -= 1
    aseRC_ind -= 1

## Build dir names
raw_dir = str(raw_ind).zfill(2) + "_rawData"
trim_dir = str(trim_ind).zfill(2) + "_trim"
align_dir = str(align_ind).zfill(2) + "_align"
addRG_dir = str(addRG_ind).zfill(2) + "_addRG"
markDuplicates_dir = str(markDuplicates_ind).zfill(2) + "_markDuplicates"
splitNCigar_dir = str(splitNCigar_ind).zfill(2) + "_splitNCigar"
bqsr_dir = str(bqsr_ind).zfill(2) + "_bqsr"
variants_dir = str(variants_ind).zfill(2) + "_variants"
if config["umi"]["add_header"]["activate"]:
    addUmis_dir = str(addUmis_ind).zfill(2) + "_addUmis"
if config["umi"]["activate"]:
    groupUmis_dir = str(groupUmis_ind).zfill(2) + "_groupUmis"
if config["merge_samples"]["activate"]:
    mergeSamples_dir = str(mergeSamples_ind).zfill(2) + "_mergeSamples"
if config["bootstrap_known_variants"]["activate"]:
    knownVariants_dir = str(knownVariants_ind).zfill(2) + "_knownVariants"
if config["ase_counts"]["activate"]:
    wasp_dir = str(wasp_ind).zfill(2) + "_wasp"
    aseRC_dir = str(aseRC_ind).zfill(2) + "_aseRC"

####
## Sample names
####

if config["samples_tsv"]:
    samples_df = pd.read_csv(config["samples_tsv"], sep="\t")
    samples = list(dict.fromkeys(samples_df["sample"])) ## Remove duplicates
else:
    samples = os.listdir(os.path.join("results", raw_dir, "fastq"))
    samples = [
        sample.replace(config["fastq_ext"], "") for sample in samples
    ]
    for tag in config["pair_tags"]:
        samples = [sample.replace(tag, "") for sample in samples]
    if config["merge_samples"]["activate"]:
        for tag in config["merge_samples"]["tags"]:
            samples = [sample.replace(tag, "") for sample in samples]
    if config["umi"]["add_header"]["activate"]:
        samples = [sample.replace(config["umi"]["add_header"]["tag"], "") for sample in samples]
    samples = list(dict.fromkeys(samples))  ## Remove duplicates

####
## Reference filenames, paths and urls for downloading
####

refGz = os.path.join(
    ".".join([
        config["ref"]["species"].capitalize(),
        config["ref"]["assembly"],
        "dna.primary_assembly.fa.gz"
    ])
)
refGz_path = os.path.join("resources", refGz)
refFa = refGz.rstrip(".gz")
refFa_path = os.path.join("resources", refFa)
ref_url = os.path.join(
    "rsync://ftp.ensembl.org/ensembl/pub",
    "release-" + str(config["ref"]["ensembl_release"]),
    "fasta",
    config["ref"]["species"],
    "dna",
    refGz
)
gtf = os.path.join(
    ".".join([
        config["ref"]["species"].capitalize(),
        config["ref"]["assembly"],
        str(config["ref"]["ensembl_release"]),
        "chr.gtf.gz"
    ])
)
gtf_path = os.path.join("resources", gtf)
gtf_url = os.path.join(
    "rsync://ftp.ensembl.org/ensembl/pub",
    "release-" + str(config["ref"]["ensembl_release"]),
    "gtf",
    config["ref"]["species"],
    gtf
)
knownVariants = os.path.join(".".join([config["ref"]["species"], "vcf.gz"]))
knownVariants_path = os.path.join("resources", knownVariants)
knownVariants_url = os.path.join(
    "rsync://ftp.ensembl.org/ensembl/pub",
    "release-" + str(config["ref"]["ensembl_release"]),
    "variation/vcf",
    config["ref"]["species"],
    knownVariants
)

####
## Input functions for rules affected by optional rules
####

def trim_inputs(wildcards):
    return {
        "R1": os.path.join(
            "results",
            addUmis_dir if config["umi"]["add_header"]["activate"] else raw_dir,
            "fastq",
            wildcards.SAMPLE + wildcards.MERGETAG + config["pair_tags"][0] + config["fastq_ext"]
        ),
        "R2": os.path.join(
            "results",
            addUmis_dir if config["umi"]["add_header"]["activate"] else raw_dir,
            "fastq",
            wildcards.SAMPLE + wildcards.MERGETAG + config["pair_tags"][1] + config["fastq_ext"]
        ),
    }

def mergeSamples_inputs(wildcards):
    return {
        "bam" : expand(os.path.join(
            "results",
            groupUmis_dir if config["umi"]["activate"] else addRG_dir,
            "bam",
            wildcards.SAMPLE + "{MERGETAG}" + ".bam"
        ), MERGETAG = config["merge_samples"]["tags"]),
        "bamIndex" : expand(os.path.join(
            "results",
            groupUmis_dir if config["umi"]["activate"] else addRG_dir,
            "bam",
            wildcards.SAMPLE + "{MERGETAG}" + ".bam.bai"
        ), MERGETAG = config["merge_samples"]["tags"])
    }

def markDuplicates_inputs(wildcards):
    if config["merge_samples"]["activate"]:
        input_dir = mergeSamples_dir
    elif config["umi"]["activate"]:
        input_dir = groupUmis_dir
    else:
        input_dir = addRG_dir
    return {
        "bam" : os.path.join(
            "results",
            input_dir,
            "bam",
            wildcards.SAMPLE + ".bam"
        ),
        "bamIndex" : os.path.join(
            "results",
            input_dir,
            "bam",
            wildcards.SAMPLE + ".bam.bai"
        ),
    }

def knownVariants_files(wildcards):
    if config["bootstrap_known_variants"]["activate"]:
        return {
            "knownVariants" : os.path.join(
                "results",
                knownVariants_dir,
                "6_select",
                "known_variants.vcf.gz"
            ),
            "knownVariantsIndex" : os.path.join(
                "results",
                knownVariants_dir,
                "6_select",
                "known_variants.vcf.gz.tbi"
            )
        }
    else:
        return {
            "knownVariants" : knownVariants_path,
            "knownVariantsIndex" : knownVariants_path + ".tbi",
            }

####
## All file endpoints
####

def workflow_outputs():

    outputs = []

    ## FastQC reports
    fqc_raw = expand(
        os.path.join("results", raw_dir, "FastQC/{SAMPLE}{MERGETAG}{PAIRTAG}_fastqc.{EXT}"),
        SAMPLE=samples,
        MERGETAG=config["merge_samples"]["tags"],
        PAIRTAG=config["pair_tags"],
        EXT=["html", "zip"]
    )
    outputs.extend(fqc_raw)
    fqc_trim = expand(
        os.path.join("results", trim_dir, "FastQC/{SAMPLE}{MERGETAG}{PAIRTAG}_fastqc.{EXT}"),
        SAMPLE=samples,
        MERGETAG=config["merge_samples"]["tags"],
        PAIRTAG=config["pair_tags"],
        EXT=["html", "zip"]
    )
    outputs.extend(fqc_trim)
    fqc_align = expand(
        os.path.join("results", align_dir, "FastQC/{SAMPLE}{MERGETAG}_fastqc.{EXT}"),
        SAMPLE=samples,
        MERGETAG=config["merge_samples"]["tags"],
        EXT=["html", "zip"]
    )
    outputs.extend(fqc_align)

    ## BQSR summary for QC purposes
    bqsr_analyzeCovariates = expand(
        os.path.join("results", bqsr_dir, "recal/{SAMPLE}.analyzeCovariates.csv"),
        SAMPLE=samples
    )
    outputs.extend(bqsr_analyzeCovariates)

    ## Variants
    variants = [os.path.join("results", variants_dir, "6_select/all_samples.vcf.gz")]
    outputs.extend(variants)

    ## ASE read counts
    if config["ase_counts"]["activate"]:
        aseRC = expand(
            os.path.join("results", aseRC_dir, "{DIR}", "{SAMPLE}.tsv"),
            DIR=["wasp", "no_wasp"],
            SAMPLE=samples
        )
        outputs.extend(aseRC)

    return(outputs)