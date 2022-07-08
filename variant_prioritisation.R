# srun with 20G
library(tidyverse)
library(UpSetR)

#ofg.dir <- '/scratch/users/k2142172/outputs/ofg/'
ofg.dir <- 'C:/Users/Prasanth/Documents/ofg_wes/'

#===============================================================================
# Data prep
#===============================================================================

setwd(ofg.dir)
# load all unfiltered data with genotypes
#load(paste0(ofg.dir, 'processed_chr.RData'))
# load premade csv for filtered data
ofg <- read.table('all_ofg_filtered_variants.tsv', header = T)
ofg <- select(ofg, -c(contains('EGAN')))

# genes published as supposedly false red flags
exclude.genes <- read_table('flagged_genes_to_exclude.txt', col_names = F)

# remove initially unnecessary cols from filtered dataframes
key.cols <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'Allele', 'Consequence', 'IMPACT', 'SYMBOL',
'BIOTYPE', 'Existing_variation', 'MAX_AF', 'CLIN_SIG', 'CADD_PHRED', 'Conservation', 'SIFT_score',
'PolyPhen_score', 'VCF_AF', 'AQ', 'AC', 'AN')

ofg.keycols <- select(ofg, all_of(key.cols))
ofg.keycols <- ofg.keycols[(ofg.keycols$MAX_AF < 0.01) | (is.na(ofg.keycols$MAX_AF)), ]

#===============================================================================
# Filter within VCF data
#===============================================================================


# top genes by summed Allele counts in VCF
top.filt.genes <- ofg.keycols %>% group_by(SYMBOL) %>% summarise(count = sum(AC)) %>% 
  filter(count > 10) %>% arrange(desc(count))

# top variant hits filtering by various metrics in VEP annotations
top.hits <- list()
# impact is redundant, they are all moderate high
#top.hits$IMPACT <- ofg[ofg$IMPACT %in% c('MODERATE', 'HIGH'), 'ID']
top.hits$novel <- ofg.keycols[is.na(ofg.keycols$Existing_variation), 'ID']
#top.hits$indels <- ofg.keycols[ofg.keycols$VARIANT_CLASS %in% c('deletion', 'indel', 'insertion'), 'ID']
#top.hits$gene_pheno <- ofg.keycols[ofg.keycols$GENE_PHENO == 1, 'ID']
#top.hits$rare <- ofg.keycols[(ofg.keycols$MAX_AF < 0.01) | (is.na(ofg.keycols$MAX_AF)), 'ID']
top.hits$high_cadd <- ofg.keycols[ofg.keycols$CADD_PHRED > 20, 'ID']
top.hits$gene <- ofg.keycols[ofg.keycols$SYMBOL %in% top.filt.genes$SYMBOL, 'ID']
top.hits$polyphen <- ofg.keycols[ofg.keycols$PolyPhen_score > 0.8, 'ID']
top.hits$sift <- ofg.keycols[ofg.keycols$SIFT_score <= 0.05, 'ID']
top.hits$gerp <- ofg.keycols[ofg.keycols$Conservation >= 2, 'ID']

# upset plot to show interactions
# upset(fromList(top.hits), nsets = 10, nintersects = NA, scale.intersections = 'log2', 
#       scale.sets = 'log2')

var.hits <- Reduce(intersect, list(top.hits$gene, top.hits$high_cadd, 
                                   top.hits$polyphen, top.hits$sift, top.hits$gerp))
ofg.hits <- ofg.keycols[ofg.keycols$ID %in% var.hits, ]
ofg.hits %>% group_by(SYMBOL) %>% summarise(count = sum(AC)) %>% arrange(desc(count))


#===============================================================================
# Filter by specific targets
#===============================================================================


# ofg associated genes from Phil
ofg.genes <- read.table('ofg_associated_genes.txt')

# check overlapping ids hit gene NOD2 in filtered data
nod2 <- filter(ofg.keycols, SYMBOL == 'NOD2', (VCF_AF > MAX_AF) | is.na(MAX_AF))
nod2p <- filter(ofg.keycols, SYMBOL == 'NOD2')
#write.table(nod2p, 'nod2_prioritised_variants.tsv', quote = F, sep = '\t', row.names = F)

# ofg filtered variants found in associated genes
ofg.var.in.genes <- filter(ofg, SYMBOL %in% ofg.genes$V1)
#write.table(ofg.var.in.genes, 'prioritised_variants_in_known_ofg_genes.tsv', sep = '\t', quote = F, row.names = F)
# ofg associated genes not found in our filtered variant prioritisation data frame
ofg.genes.not.in.filt <- ofg.genes$V1[!ofg.genes$V1 %in% names(table(ofg.var.in.genes$SYMBOL))]


#===============================================================================
# Compare to UK biobank exomes
#===============================================================================


# write data out with IDs in ukbb format
as.ukbb.ids <- ofg.hits
as.ukbb.ids$ID <- gsub('chr', '', as.ukbb.ids$ID)
as.ukbb.ids$ID <- gsub('_', ':', as.ukbb.ids$ID)
as.ukbb.ids <- split(as.ukbb.ids, as.ukbb.ids$CHROM)
as.ukbb.ids <- lapply(as.ukbb.ids, function(x) unlist(strsplit(x$ID, ';')))
# for (i in seq(1, 22)){
#   write.table(as.ukbb.ids[[i]], paste0('top_hits_ukbb/vars_as_ukbb_ids_', i, '.txt'), quote = F, row.names = F, col.names = F)
# }

# load vcf data file for vars of interest, takes a while since so many columns
vcf <- vroom::vroom('top_hits_ukbb/vars_from_ukbb_22.vcf', delim = '\t', comment = '##')
# remove info fields, convert genotypes to integers to make matrix
vcf.mat <- vcf[, -c(1:9)]
#vcf.mat <- vcf.mat[, 1:20000]

vcfToMatGsubs <- function(x){
  xa <- gsub(':.*$', '', x)
  
  if(xa == '0/0'){
    xb <- gsub('0/0', 0, xa)
  } else if(xa == '0/1'){
    xb <- gsub('0/1', 1, xa)
  } else if(xa == '1/0'){
    xb <- gsub('1/0', 1, xa)
  } else if(xa == '1/1'){
    xb <- gsub('1/1', 2, xa)
  } else if(xa == './.'){
    xb <- gsub('\\./\\.', NA, xa)
  } else{
    print(xa)
    stop('unrecognised genotype')
  }
  
  return(xb)
}


Sys.time()
vcf.mat <- apply(vcf.mat, c(1,2), vcfToMatGsubs)
Sys.time()
print('done')

vcf.mat <- apply(vcf.mat, c(1,2), as.integer)
vcf.mat <- as.matrix(vcf.mat)
# divide total counts for alleles by total allele number to get frequency
vcf.rowsums <- rowSums(vcf.mat, na.rm = T) / (ncol(vcf.mat)*2)
vcf.df <- data.frame(vcf[, 1:8], vcf.rowsums)
vcf.df$ID <- gsub(':', '_', vcf.df$ID)
vcf.df$ID <- gsub('^', 'chr', vcf.df$ID)
ofg.hits.ukbb <- left_join(ofg.hits, vcf.df, by = 'ID')
ofg.hits.ukbb.nona <- ofg.hits.ukbb[!is.na(ofg.hits.ukbb$vcf.rowsums), ]
ofg.hits.ukbb.nona$diff <- ofg.hits.ukbb.nona$VCF_AF / ofg.hits.ukbb.nona$vcf.rowsums
ofg.hits.ukbb.nona[which((ofg.hits.ukbb.nona$VCF_AF / ofg.hits.ukbb.nona$vcf.rowsums) > 5), ] %>% View()

#===============================================================================
# GWAS comparisons
#===============================================================================


# load gwas catalgoue data for crohns
gwascat <- read_csv('crohns_gwas_catalogue_EFO_0000384.csv')

# fix formatting
gw <- separate(gwascat, Variant, c('Existing_variation', 'Allele'), '-<') %>%
  mutate(Allele = gsub('b>', '', Allele)) %>%
  mutate(Allele = gsub('<.*$', '', Allele))

# find overlapping rsIDs
prts <- filter(ng, Existing_variation %in% gw$Existing_variation)
# gives 1 variant, an insertion in NOD2
# in gw catalogue, this variant with rsID is identified in 3 crohns papers
