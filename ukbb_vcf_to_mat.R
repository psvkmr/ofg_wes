# load vcf data file for vars of interest, takes a while since so many columns
vcfs <- list()
Sys.time()
for (i in seq(1, 22)){
  vcf <- read_delim(paste0('/scratch/users/k2142172/outputs/ofg/top_hits_ukbb/vars_from_ukbb_', i, '.vcf'), delim = '\t', comment = '##')
  vcfs[[i]] <- vcf
}
Sys.time()

# remove info fields, convert genotypes to integers to make matrix
vcf.mats <- lapply(vcfs, function(x) x[, -c(1:9)])
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
vcf.mats2 <- lapply(vcf.mats, function(x) apply(x, c(1,2), vcfToMatGsubs))
Sys.time()
vcf.mats2 <- lapply(vcf.mats2, function(x) apply(x, c(1,2), as.integer))
vcf.mats2 <- lapply(vcf.mats2, as.matrix)
save(vcf.mats2, file = '/scratch/users/k2142172/outputs/ofg/top_hits_ukbb/ofg_ukbb_vcf_mat.RData')
print('done')
