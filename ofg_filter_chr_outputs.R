library(tidyverse)

terms <- read.csv('/scratch/users/k2142172/resources/ensembl_consequence_rankings.csv')
ofg.dir <- '/scratch/users/k2142172/outputs/ofg/'

load(paste0(ofg.dir, 'processed_chr.RData'))

filterData <- function(df){
    df <- df %>% select(-c(contains('EGAN'), 'AFR_AF', 'AMR_AF', 'EAS_AF', 'SAS_AF',
        'AA_AF', 'gnomAD_AFR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_SAS_AF')) %>%
        filter(IMPACT %in% c('HIGH', 'MODERATE'))
    df <- as.data.frame(apply(df, 2, function(x) gsub("\"\"", NA, x)))
    df <- mutate(df, 'SIFT_score' = as.numeric(str_extract(df$SIFT, '0\\.*[0-9]*')),
         'PolyPhen_score' = as.numeric(str_extract(df$PolyPhen, '0\\.*[0-9]*')),
         'CADD_PHRED' = as.numeric(CADD_PHRED))
    df <- filter(df, (is.na(ID)) | (is.na(MAX_AF) & is.na(AF)) | (is.na(CADD_PHRED) & is.na(SIFT_score) & is.na(PolyPhen_score))
        | (MAX_AF < 0.001) | (CADD_PHRED > 30) | (SIFT_score < 0.05) | (PolyPhen_score > 0.15))
    return(df)
}

ofg.filt <- lapply(ofg, filterData)

for (i in 1:length(ofg.filt)){
  write.table(ofg.filt[[i]], paste0(ofg.dir, 'ofg_', names(ofg.filt)[[i]], '_filtered_variants.tsv'), sep = '\t', row.names = F, quote = F)
}
