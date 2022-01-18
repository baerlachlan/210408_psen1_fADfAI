####
## Manage config
####

if not config["paired_end"]["activate"]:
    config["paired_end"]["tags"] = None
if not config["merge_lanes"]["activate"]:
    config["merge_lanes"]["tags"] = None

####
## Directory structure
####

## Define dir indexes
raw_ind = 0
merge_ind = 1
trim_ind = 2
align_ind = 3
counts_ind = 4
## Index decreases by 1 for any steps following optional steps that are not run
if not config["merge_lanes"]["activate"]:
    trim_ind -= 1
    align_ind -= 1
    counts_ind -= 1

## Build dir names
raw_dir = str(raw_ind).zfill(2) + "_rawData"
merge_dir = str(merge_ind).zfill(2) + "_merge"
trim_dir = str(trim_ind).zfill(2) + "_trim"
align_dir = str(align_ind).zfill(2) + "_align"
counts_dir = str(counts_ind).zfill(2) + "_counts"

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
    if config["paired_end"]["activate"]:
        for tag in config["paired_end"]["tags"]:
            samples = [sample.replace(tag, "") for sample in samples]
    if config["merge_lanes"]["activate"]:
        for tag in config["merge_lanes"]["tags"]:
            samples = [sample.replace(tag, "") for sample in samples]
    samples = list(dict.fromkeys(samples))  ## Remove duplicates

####
## Reference files
####

species = config["species"]
assembly = config["assembly"]
ensembl_release = config["ensembl_release"]

refGz = os.path.join(
    ".".join([
        species.capitalize(),
        assembly,
        "dna.primary_assembly.fa.gz"
    ])
)
refGz_path = os.path.join("resources", refGz)
refFa = refGz.rstrip(".gz")
refFa_path = os.path.join("resources", refFa)
ref_url = os.path.join(
    "rsync://ftp.ensembl.org/ensembl/pub",
    "release-" + str(ensembl_release),
    "fasta",
    species,
    "dna",
    refGz
)
gtf = os.path.join(
    ".".join([
        species.capitalize(),
        assembly,
        str(ensembl_release),
        "chr.gtf.gz"
    ])
)
gtf_path = os.path.join("resources", gtf)
gtf_url = os.path.join(
    "rsync://ftp.ensembl.org/ensembl/pub",
    "release-" + str(ensembl_release),
    "gtf",
    species,
    gtf
)

####
## Input functions
####

def trim_inputs_pe(wildcards):
    # if config["merge_lanes"]["activate"]:
    return {
        "R1": os.path.join(
            "results",
            merge_dir if config["merge_lanes"]["activate"] else raw_dir,
            "fastq",
            wildcards.SAMPLE + config["paired_end"]["tags"][0] + config["fastq_ext"]
        ),
        "R2": os.path.join(
            "results",
            merge_dir if config["merge_lanes"]["activate"] else raw_dir,
            "fastq",
            wildcards.SAMPLE + config["paired_end"]["tags"][1] + config["fastq_ext"]
        )
    }

####
## All file endpoints
####

def workflow_outputs():

    outputs = []

    ## FastQC reports
    fqc_raw = expand(
        os.path.join("results", raw_dir, "FastQC/{SAMPLE}{LANETAG}{PAIRTAG}_fastqc.{EXT}"),
        SAMPLE=samples,
        LANETAG=config["merge_lanes"]["tags"],
        PAIRTAG=config["paired_end"]["tags"],
        EXT=["html", "zip"],
    )
    outputs.extend(fqc_raw)
    # fqc_trim = expand(
    #     os.path.join(trim_dir, "FastQC/{SAMPLE}{MERGETAG}{PAIRTAG}_fastqc.{EXT}"),
    #     SAMPLE=samples,
    #     MERGETAG=lane_tags,
    #     PAIRTAG=pair_tags,
    #     EXT=["html", "zip"]
    # )
    # outputs.extend(fqc_trim)
    # fqc_align = expand(
    #     os.path.join(align_dir, "FastQC/{SAMPLE}{MERGETAG}_fastqc.{EXT}"),
    #     SAMPLE=samples,
    #     MERGETAG=lane_tags,
    #     EXT=["html", "zip"]
    # )
    # outputs.extend(fqc_align)

    ## Tmp
    # tmp = expand(os.path.join("results", trim_dir, "fastq", "{SAMPLE}{PAIRTAG}" + config["fastq_ext"]), SAMPLE=samples, PAIRTAG=config["paired_end"]["tags"])
    # tmp = expand(os.path.join("results", align_dir, "bam", "{SAMPLE}.bam"), SAMPLE=samples)
    tmp = os.path.join("results", counts_dir, "counts.out")
    outputs.extend([tmp])

    return(outputs)