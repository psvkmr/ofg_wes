library(tidyverse)

ofg.dir <- '/scratch/users/k2142172/outputs/ofg/'
#var.file.names <- list.files(path = ofg.dir, pattern = '*filtered_variants.tsv', full.names = T)
#ofg <- lapply(var.file.names, read_table)
load(paste0(ofg.dir, 'processed_chr.RData'))

# genes to exclude
exclude.genes <- read_table('/scratch/users/k2142172/resources/flagged_genes_to_exclude.txt', col_names = F)

# should be done already
filterData <- function(df){
    df <- df %>% select(-c(contains('EGAN'), 'AFR_AF', 'AMR_AF', 'EAS_AF', 'SAS_AF',
        'AA_AF', 'gnomAD_AFR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_SAS_AF')) %>%
        filter(IMPACT %in% c('HIGH', 'MODERATE'))
    df <- as.data.frame(apply(df, 2, function(x) gsub("^$", NA, x)))
    df$CHROM <- gsub('chr', '', df$CHROM)
    df <- mutate(df, across(.cols = c('CHROM', 'POS', 'AF',	'EUR_AF',	'EA_AF',	'gnomAD_AF',	'gnomAD_AMR_AF',	'gnomAD_NFE_AF',	'gnomAD_OTH_AF',	'MAX_AF',
                                      'MaxEntScan_alt',	'MaxEntScan_diff',	'MaxEntScan_ref',	'CADD_PHRED',	'CADD_RAW',	'Conservation',
                                      'SIFT_score',	'PolyPhen_score',	'VCF_AF',	'AQ',	'AC',	'AN'), as.numeric))
    df <- filter(df, AN > 400)
    df <- filter(df, (is.na(Existing_variation)) | (is.na(MAX_AF) & is.na(AF)) | (is.na(CADD_PHRED) & is.na(SIFT_score) & is.na(PolyPhen_score))
        | (MAX_AF < 0.001) | (CADD_PHRED > 30) | (SIFT_score < 0.05) | (PolyPhen_score > 0.15))
    return(df)
}

ofg.filt <- lapply(ofg, filterData)
all.ofg.filt <- Reduce(rbind, ofg.filt) %>% arrange(CHROM, POS)

for (i in 1:length(ofg.filt)){
  write.table(ofg.filt[[i]], paste0(ofg.dir, 'ofg_', names(ofg.filt)[[i]], '_filtered_variants.tsv'), sep = '\t', row.names = F, quote = F)
}
write.table(all.ofg.filt, paste0(ofg.dir, 'all_ofg_filtered_variants.tsv'), sep = '\t', row.names = F, quote = F)

# more stringent
stringentFilt <- function(df){
  df <- filter(df, ((MAX_AF < 0.0001) | (is.na(Existing_variation) & is.na(MAX_AF))))
  df <- df[!df$SYMBOL %in% exclude.genes$X1, ]
  return(df)
}

ofg.stringent <- lapply(ofg.filt, stringentFilt)
all.ofg.stringent <- Reduce(rbind, ofg.stringent) %>% arrange(CHROM, POS)

for (i in 1:length(ofg.stringent)){
  write.table(ofg.stringent[[i]], paste0(ofg.dir, 'ofg_', names(ofg.stringent)[[i]], '_stringent_variants.tsv'), sep = '\t', row.names = F, quote = F)
}
write.table(all.ofg.stringent, paste0(ofg.dir, 'all_ofg_stringent_variants.tsv'), sep = '\t', row.names = F, quote = F)
