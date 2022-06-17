library(tidyverse)

ofg.dir <- '/scratch/users/k2142172/outputs/ofg/'
ofg.files <- list.files(path = ofg.dir, pattern = 'ofg_chr.*.vcf.gz.batch.vep.vcf.processed.tsv',
full.names = T)
ofg <- lapply(ofg.files, read_table)
ofg.names <- str_extract(ofg.files, 'chr[0-9]+')
names(ofg) <- ofg.names

save.image(paste0(ofg.dir, 'processed_chr.RData'))
