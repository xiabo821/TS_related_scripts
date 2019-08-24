#!/usr/bin/env Rscript
# Yun Yan (yy1533@nyu.edu)

args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(magrittr); library(gtools)})

workdir <- args[1]
fout_count_unstrand <- args[2]
fout_count_stranded <- args[3]
# workdir <- '.'
# setwd('/Users/yunyan/Downloads/icm-count')
infpath_l <- list.files(
  file.path(workdir, 'count/vcf192'),
  full.names = T, recursive = F) %>% sort()

to_basepair <- function(x = c('A', 'T', 'G', 'C', 'N')){
  x <- match.arg(x)
  switch (x,
    A = 'T', T = 'A', G = 'C', C = 'G', N = 'N'
  )
}

#define %not_in% funtion:
'%not_in%' <- function(x,y)!('%in%'(x,y))


# generate all 4x12x5=192 possible snp contexts
generate_all_192_snp_contexts <- function(){
  snp4_expected <- NULL
  snp_expected <- gtools::permutations(n=4, r=2, v=c('A', 'C', 'G','T')) #12
  snp_expected_m <- gtools::permutations(n=4, r=2, v=c('A', 'C', 'G','T'), repeats.allowed = T) #25
  for (r in seq_len(nrow(snp_expected))) {
    snp4_snp <- do.call(
      rbind,
      replicate(nrow(snp_expected_m), snp_expected[r, , drop=F], simplify = F))
    snp4_expected_r <- cbind(
      snp_expected_m[, 1, drop=F],
      snp4_snp[, 1, drop=F],
      snp_expected_m[, 2, drop=F],
      snp4_snp[, 2, drop=F]
    )
    snp4_expected <- rbind(snp4_expected, snp4_expected_r)
  }
  snp4_expected <- apply(snp4_expected, 1, function(r) paste0(r, collapse = '')) %>%
    unique() %>% sort()
  return(snp4_expected)
}

snp_192_contexts <- generate_all_192_snp_contexts()
str(snp_192_contexts)

#Get single base SNV type from the tri-base reference SNV type. 
#get_snv_type_from_context <- function(x){
#  # get_snv_type_from_context('AGTC') ==> 'GC'
##  res <- strsplit(x=x, split = '') %>% unlist()
#  paste0(res[2], res[4])
#}

#fix SNV type depends on the strand info. 
fix_snv_type_with_strand <- function(x, strand_int=c(-1,1)){
  # AACG,-1 ==> _T_C ==> GTTC
  if (strand_int == 1) return(x)
  x <- strsplit(x=x, split = '') %>% unlist() # A,A,C,G
  res <- c(to_basepair(x[3]), to_basepair(x[2]), to_basepair(x[1]),
           to_basepair(x[4])) # G, T, T, C
  return(paste0(res, collapse = ''))
}

  
#generate an all zero output variable with column names being snp_192_contexts
res_unstrand <- NULL
cnt <- read.delim(file = infpath_l[[1]], header = F, stringsAsFactors = F)
colnames(cnt) <- c('chr_name', 'start', 'end',
                   'gene_id', 'gene_name', 'strand_int', 'cnt')
res_unstrand <- as_tibble(cnt[, 1:6])
for (snp_c in snp_192_contexts) res_unstrand[[snp_c]] <- 0
# dim(res_unstrand)




########### generate the unstranded SNV counts ###########
cat('\n\n##############################################\n')
cat('\n[Start the unstranded counting.]\n\n')

validbase = c('A', 'C', 'G','T') 

snp_contexts_avail <- c()
for (i in seq_len(length(infpath_l))){
  infpath <- infpath_l[[i]]
  snp_c_raw <- basename(infpath)
  # exclude NXxY, xXNY, orNXNY (X>Y SNV), "N" here refers to the unidentified genome reference base.
  snp_c_raw_l <- strsplit(snp_c_raw, split='') %>% unlist()
  if (snp_c_raw_l[1] %not_in% validbase | snp_c_raw_l[2] %not_in% validbase | snp_c_raw_l[3] %not_in% validbase | snp_c_raw_l[4] %not_in% validbase){
    cat('excluding ', infpath, '.\n')
    next()
  }
  cat('filling ', infpath, '.\n')
  cnt <- read.delim(file = infpath, header = F, stringsAsFactors = F)

  snp_contexts_avail <- c(snp_contexts_avail, snp_c_raw)
  colnames(cnt) <- c('chr_name', 'start', 'end',
                     'gene_id', 'gene_name', 'strand_int', 'cnt')
  stopifnot(all.equal(res_unstrand$gene_id, cnt$gene_id))
  res_unstrand[[snp_c_raw]] <- cnt$cnt
  cat('[Done with ',snp_c_raw,']','\n')
}

cat(length(snp_contexts_avail), 'available snp contexts:\n')
str(snp_contexts_avail)
#prepare for negative strand processing.
snp_192_contexts_strand_dict <- sapply(snp_contexts_avail, function(x){
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
res_stranded_neg_ginfo <- res_stranded_neg[, 1:6]				#minus-strand gene info
res_stranded_neg_val   <- res_stranded_neg[, snp_contexts_avail]  #minus-strand gene SNV counts
# head(colnames(res_stranded_neg_val))
res_stranded_neg_val <- res_stranded_neg_val[, as.character(snp_192_contexts_strand_dict)]
# head(colnames(res_stranded_neg_val))
colnames(res_stranded_neg_val) <- snp_contexts_avail
res_stranded_neg <- cbind(res_stranded_neg_ginfo, res_stranded_neg_val)
res_stranded <- rbind(res_stranded_pos, res_stranded_neg)
cat('\n[Done with the stranded counting.]\n')
cat('R variable name of the stranded counting results is [res_stranded].\n')




########### Aggregate by genes; Quality check of SNV counts ###########
cat('\n\n##############################################\n')
cat('\n[Start aggregation by genes and Quality check.]\n\n')

res_unstrand_f <- cbind(res_unstrand$gene_id, res_unstrand[, snp_contexts_avail])
colnames(res_unstrand_f) <- c('gene_id',snp_contexts_avail)
res_unstrand_final <- aggregate(. ~ gene_id, res_unstrand_f, sum)
res_stranded_f <- cbind(res_stranded$gene_id, res_stranded[, snp_contexts_avail])
colnames(res_stranded_f) <- c('gene_id',snp_contexts_avail)
res_stranded_final <- aggregate(. ~ gene_id, res_stranded_f, sum)

idx <- match(res_unstrand_final$gene_id, res_stranded_final$gene_id)
res_stranded_final <- res_stranded_final[idx, ]
stopifnot(all.equal(res_unstrand_final$gene_id, res_stranded_final$gene_id))

cat('\n[Done with aggregation and quolity checking.]\n')
cat('Export results... \n')
write_tsv(x=res_unstrand_final, path = fout_count_unstrand)
write_tsv(x=res_stranded_final, path = fout_count_stranded)
cat('[DONE!!!]\n')
