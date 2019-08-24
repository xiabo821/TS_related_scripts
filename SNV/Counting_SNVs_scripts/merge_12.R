#!/usr/bin/env Rscript
# Yun Yan (yy1533@nyu.edu)

args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(magrittr); library(gtools)})

workdir <- args[1]
fout_count_unstrand <- args[2]
fout_count_stranded <- args[3]
# workdir <- '.'
# setwd('/Users/yunyan/Downloads/icm_count')
infpath_l <- list.files(
  file.path(workdir, 'count/vcf12'),
  full.names = T, recursive = F) %>% sort()

to_basepair <- function(x = c('A', 'T', 'G', 'C', 'N')){
  x <- match.arg(x)
  switch (x,
    A = 'T', T = 'A', G = 'C', C = 'G', N = 'N'
  )
}

#define %not_in% funtion:
'%not_in%' <- function(x,y)!('%in%'(x,y))

# generate all 4x3=12 possible snps
generate_all_12_snps <- function(){
  snp_expected <- gtools::permutations(n=4, r=2, v=c('A', 'C', 'G','T')) #12

  snp_expected <- apply(snp_expected, 1, function(r) paste0(r, collapse = '')) %>%
    unique() %>% sort()
  return(snp_expected)
}
snp_12 <- generate_all_12_snps()
str(snp_12)

#fix SNV type depends on the strand info. 
fix_snv_type_with_strand <- function(x, strand_int=c(-1,1)){
  # fix_snv_type_with_strand('AG', -1)
  # AG,-1 ==> TC
  if (strand_int == 1) return(x)
  x <- strsplit(x=x, split = '') %>% unlist() # A,G
  res <- c(to_basepair(x[1]), to_basepair(x[2])) # T, C
  return(paste0(res, collapse = ''))
}
  

#generate an all zero output variable with column names being snp_12
res_unstrand <- NULL
cnt <- read.delim(file = infpath_l[[1]], header = F, stringsAsFactors = F)
colnames(cnt) <- c('chr_name', 'start', 'end',
                   'gene_id', 'gene_name', 'strand_int', 'cnt')
res_unstrand <- as_tibble(cnt[, 1:6])
for (snp in snp_12) res_unstrand[[snp]] <- 0
# dim(res_unstrand)




########### generate the unstranded SNV counts ###########
cat('\n\n##############################################\n')
cat('\n[Start the unstranded counting.]\n\n')

validbase = c('A', 'C', 'G','T') 
snp_avail <- c()
for (i in seq_len(length(infpath_l))){
  infpath <- infpath_l[[i]]
  snp_raw <- basename(infpath)
  # exclude NN, xN, Nx
  snp_raw_l <- strsplit(snp_raw, split='') %>% unlist()
  if (snp_raw_l[1] %not_in% validbase | snp_raw_l[2] %not_in% validbase){
    cat('excluding ', infpath, '.\n')
    next()
  }
  cat('filling ', infpath, '.\n')
  cnt <- read.delim(file = infpath, header = F, stringsAsFactors = F)

  snp_avail <- c(snp_avail, snp_raw)
  colnames(cnt) <- c('chr_name', 'start', 'end',
                     'gene_id', 'gene_name', 'strand_int', 'cnt')
  stopifnot(all.equal(res_unstrand$gene_id, cnt$gene_id))
  res_unstrand[[snp_raw]] <- cnt$cnt
  cat('[Done with ',snp_raw,']','\n')
}

cat(length(snp_avail), 'available snps:\n')
str(snp_avail)

#prepare for negative strand processing.
snp_12_strand_dict <- sapply(snp_avail, function(x){
  fix_snv_type_with_strand(x=x, strand = -1)})


cat('\n[Done with the unstranded counting.]\n')
cat('R variable name of the unstranded counting results is [res_unstrand].\n')


########### generate the stranded SNV counts ###########
cat('\n\n##############################################\n')
cat('\n[Start the stranded counting.]\n\n')

res_stranded <- res_unstrand
g_is_pos <- res_stranded$strand_int == 1; stopifnot(!any(is.na(g_is_pos)))
res_stranded_pos <- res_stranded[g_is_pos, ]
res_stranded_neg <- res_stranded[!g_is_pos, ]
res_stranded_neg_ginfo <- res_stranded_neg[, 1:6]     #minus-strand gene info
res_stranded_neg_val   <- res_stranded_neg[, snp_avail]  #minus-strand gene SNV counts
# head(colnames(res_stranded_neg_val))
res_stranded_neg_val <- res_stranded_neg_val[, as.character(snp_12_strand_dict)]
# head(colnames(res_stranded_neg_val))
colnames(res_stranded_neg_val) <- snp_avail
res_stranded_neg <- cbind(res_stranded_neg_ginfo, res_stranded_neg_val)
res_stranded <- rbind(res_stranded_pos, res_stranded_neg)
cat('\n[Done with the stranded counting.]\n')
cat('R variable name of the stranded counting results is [res_stranded].\n')



########### Aggregate by genes; Quality check of SNV counts ###########
cat('\n\n##############################################\n')
cat('\n[Start aggregation by genes and Quality check.]\n\n')

res_unstrand_f <- cbind(res_unstrand$gene_id, res_unstrand[, snp_avail])
colnames(res_unstrand_f) <- c('gene_id',snp_avail)
res_unstrand_final <- aggregate(. ~ gene_id, res_unstrand_f, sum)
res_stranded_f <- cbind(res_stranded$gene_id, res_stranded[, snp_avail])
colnames(res_stranded_f) <- c('gene_id',snp_avail)
res_stranded_final <- aggregate(. ~ gene_id, res_stranded_f, sum)


idx <- match(res_unstrand_final$gene_id, res_stranded_final$gene_id)
res_stranded_final <- res_stranded_final[idx, ]
stopifnot(all.equal(res_unstrand_final$gene_id, res_stranded_final$gene_id))


cat('\n[Done with aggregation and quolity checking.]\n')
cat('Export results... \n')
write_tsv(x=res_unstrand_final, path = fout_count_unstrand)
write_tsv(x=res_stranded_final, path = fout_count_stranded)
cat('[DONE!!!]\n')
