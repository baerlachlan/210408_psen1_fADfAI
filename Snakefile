from os import path

from smk.workflow_config import RNASeqVariantDiscovery


analysis = RNASeqVariantDiscovery(
    species="danio_rerio",
    assembly="GRCz11",
    ensembl_release=101,
    read_length=100,
    fastq_ext=".fastq.gz",
    pair_tags=["_R1", "_R2"],
    merge_samples=True,
    merge_tags=["_L001", "_L002"],
    add_umi_header=True,
    umi_tag="_I1",
    umi_dedup=True,
    bootstrap_known_variants=False,
)

rule all:
    input:
        analysis.outputs

include: "smk/modules/refs.smk"
include: "smk/modules/intervals.smk"
if analysis.add_umi_header:
    include: "smk/modules/addUmis.smk"
include: "smk/modules/trim.smk"
include: "smk/modules/align.smk"
include: "smk/modules/fastqc.smk"
include: "smk/modules/addRG.smk"
if analysis.umi_dedup:
    include: "smk/modules/groupUmis.smk"
if analysis.merge_samples:
    include: "smk/modules/mergeSamples.smk"
include: "smk/modules/markDuplicates.smk"
include: "smk/modules/splitNCigar.smk"
if analysis.bootstrap_known_variants:
    include: "smk/modules/knownVariants.smk"
include: "smk/modules/bqsr.smk"
include: "smk/modules/variants.smk"

localrules: refs_downloadFa, refs_downloadGtf, refs_downloadKnownVariants  ## Requires internet access but not much compute so run on head node
