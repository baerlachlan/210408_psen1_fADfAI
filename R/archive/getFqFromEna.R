library(tidyverse)
library(here)

## This script gets URL's for downloading publicly available FASTQ files
## when given a study's or multiple studies' accession, and saves them as 
## a text file ready to be downloaded with `wget -i`

## Manually enter studies here
enaStudy <- c("PRJNA690794")

## Check for output directory
dirPath <- here("files")
makeDir <- !dir.exists(dirPath)
if (makeDir) {
    dir.create(dirPath)
}

## Get FASTQ download URLs
fastqDL <- lapply(enaStudy, function(x){
    read_tsv(
        paste0(
            "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
            x,
            "&result=read_run&fields=fastq_ftp&format=tsv"
        )
    ) %>% 
        .$fastq_ftp %>% 
        str_split(";") %>% 
        unlist() %>%
        vapply(function(x){paste0("ftp://", x)}, character(1), USE.NAMES = FALSE)
})

## Export URLs to their own files ready for use with bash script
if (length(enaStudy) == 1) {
    lapply(fastqDL, function(x){
        write_lines(
            x,
            here(paste0(
                "files/",
                enaStudy,
                ".txt"
            ))
        )
    })
} else {
    counter <- 0
    lapply(fastqDL, function(x){
        counter <<- counter + 1
        write_lines(
            x,
            here(paste0(
                "files/",
                counter,
                "_",
                enaStudy[counter],
                ".txt"
            ))
        )
    })
}
