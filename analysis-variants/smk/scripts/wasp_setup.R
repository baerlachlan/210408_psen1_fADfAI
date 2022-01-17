library(magrittr)
library(stringr)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(VariantAnnotation)

## Allow user to run interactively or as part of snakemake workflow
if (exists("snakemake")) {
  proj_root <- snakemake@params[["proj_root"]]
  variants_dir <- snakemake@params[["variants_dir"]]
  wasp_dir <- snakemake@params[["wasp_dir"]]
  ## Otherwise if running via snakemake, variables will be setup automatically
} else {
  proj_root <- "/hpcfs/users/a1647910/210408_psen1_fADfAI/analysis-variant_discovery"
  variants_dir <- "10_variants"
  wasp_dir <- "11_wasp"
}

vcf_file <- file.path(
  proj_root, variants_dir, "6_select", "all_samples.vcf.gz"
)
stopifnot(file.exists(vcf_file))

svp <- ScanVcfParam(info = "", geno = c("GT", "GQ"))
vcf <- suppressWarnings({
  readVcf(vcf_file, param = svp)
  # readVcf(vcf_file)
})
samples <- samples(header(vcf))

gr <- rowRanges(vcf)[,c("REF", "ALT")]
gr$REF <- as.character(gr$REF)
gr$ALT <- CharacterList(gr$ALT) %>%
  unstrsplit(sep = ",")
gr$alleles_all <- paste(gr$REF, gr$ALT, sep = ",")

hetSnvs <- lapply(samples, function(sample){
  tbl <- tibble(
    chromosome = as.vector(seqnames(gr)),
    position = start(gr),
    alleles_all = gr$alleles_all,
    GT = geno(vcf)$GT[,sample] %>%
      str_replace(., "\\|", "\\/"),
    GQ = geno(vcf)$GQ[,sample],
    sample = sample
  ) %>%
    dplyr::filter(!is.na(GQ)) %>%  # Remove un-called
    separate(
      col = "GT", into = c("allele_1", "allele_2"), sep = "/",
      remove = FALSE, convert = TRUE, fill = "right"
    ) %>%
    dplyr::filter(allele_1 != allele_2)  # Remove homozygous
  tbl$allele_1 <- map2(tbl$alleles_all, tbl$allele_1, function(x, y){
    unlist(str_split(x, ","))[y + 1]  # GT is 0 based so need to add 1
  }) %>%
    unlist()
  tbl$allele_2 <- map2(tbl$alleles_all, tbl$allele_2, function(x, y){
    unlist(str_split(x, ","))[y + 1]  # GT is 0 based so need to add 1
  }) %>%
    unlist()
  return(tbl)
})

lapply(hetSnvs, function(x){
  hetSnvs_byChr <- x %>%
    dplyr::select(chromosome, position, allele_1, allele_2, sample) %>%
    split(.[,'chromosome'])
  lapply(hetSnvs_byChr, function(x){
    chr <- unique(x$chromosome)
    sample <- unique(x$sample)
    path <- file.path(
      proj_root,
      wasp_dir,
      "1_snvs",
      sample,
      paste0(chr, ".snvs.txt.gz")
    )
    if (!dir.exists(dirname(path))) {
      dir.create(dirname(path), recursive = TRUE)
    }
    tibble(
      position = x$position,
      allele_1 = x$allele_1,
      allele_2 = x$allele_2
    ) %>%
      write_delim(file = path, delim = " ", col_names = FALSE)
  })
})
