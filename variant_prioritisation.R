# srun with 20G
library(tidyverse)

ofg.dir <- '/scratch/users/k2142172/outputs/ofg/'

# load all unfiltered data with genotypes
load(paste0(ofg.dir, 'processed_chr.RData'))
# load premade csv for filtered data
all.ofg.filt <- read.table('all_ofg_filtered_variants.tsv', header = T)
# stringent filtering of 0.0001 maf and exclude common gene hits
all.ofg.stringent <- read.table('all_ofg_stringent_variants.tsv', header = T)

# ofg associated genes
ofg.genes <- read.table('ofg_associated_genes.txt')
# load gwas catalgoue data for crohns
gwascat <- read_csv('crohns_gwas_catalogue_EFO_0000384.csv')

# remove initially unnecessary cols from filtered dataframes
ng <- select(all.ofg.filt, c(CHROM, POS, REF, ALT, Allele, Consequence, IMPACT, SYMBOL,
BIOTYPE, Existing_variation, MAX_AF, CLIN_SIG, CADD_PHRED, Conservation, SIFT_score,
PolyPhen_score, VCF_AF, AQ, AC, AN))

ngf <- select(all.ofg.stringent, c(CHROM, POS, REF, ALT, Allele, Consequence, IMPACT, SYMBOL,
BIOTYPE, Existing_variation, MAX_AF, CLIN_SIG, CADD_PHRED, Conservation, SIFT_score,
PolyPhen_score, VCF_AF, AQ, AC, AN))

# fix formatting
gw <- separate(gwascat, Variant, c('Existing_variation', 'Allele'), '-<') %>%
mutate(Allele = gsub('b>', '', Allele)) %>%
mutate(Allele = gsub('<.*$', '', Allele))

# find overlapping rsIDs
prts <- filter(ng, Existing_variation %in% gw$Existing_variation)
# gives 1 variant, an insertion in NOD2
# in gw catalogue, this variant with rsID is identified in 3 crohns papers

# check overlapping ids hit gene NOD2 in filtered data
nod2 <- filter(ng, SYMBOL == 'NOD2', (VCF_AF > MAX_AF) | is.na(MAX_AF))
nod2p <- filter(all.ofg.filt, SYMBOL == 'NOD2')
#write.table(nod2p, 'nod2_prioritised_variants.tsv', quote = F, sep = '\t', row.names = F)

# ofg filtered variants found in associated genes
ofg.var.in.genes <- filter(all.ofg.filt, SYMBOL %in% ofg.genes$V1)
#write.table(ofg.var.in.genes, 'prioritised_variants_in_known_ofg_genes.tsv', sep = '\t', quote = F, row.names = F)
# ofg associated genes not found in our filtered variant prioritisation data frame
ofg.genes.not.in.filt <- ofg.genes$V1[!ofg.genes$V1 %in% names(table(ofg.var.in.genes$SYMBOL))]

# prioritise genes by variant counts
top.filt.genes <- table(ngf$SYMBOL) %>% as.data.frame() %>% arrange(desc(Freq))
gw.genes.in.ngf <- filter(ngf, SYMBOL %in% gw$Mapped_gene)
gw.ngf.genes.frq <- filter(top.filt.genes, Var1 %in% gw.genes.in.ngf$SYMBOL)
