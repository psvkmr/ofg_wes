library(tidyverse)
library(UpSetR)

wes.resources.dir <- 'C:/Users/Prasanth/Documents/resources/wes_exclude/'

moesm4 <- read_table('12920_2014_64_MOESM4_ESM.txt')
s7 <- read_table('Table_S7_gene_exclusion_list_final.txt', skip = 1, col_names = F)

all.exclude.genes <- list(flags = moesm4$FLAGS[!is.na(moesm4$FLAGS)], 
                          omim = moesm4$OMIM[!is.na(moesm4$OMIM)],
                          hgmd = moesm4$HGMD[!is.na(moesm4$HGMD)],
                          wes = moesm4$WES[!is.na(moesm4$WES)],
                          s7 = s7$X1)

exclude.gene.intersections <- upset(fromList(exclude.gene.intersections), order.by = 'freq')

uniq.names <- unique(unlist(all.exclude.genes))

write.table(uniq.names, paste0(wes.resources.dir, 'flagged_genes_to_exclude.txt'), row.names = F, col.names = F, quote = F)
