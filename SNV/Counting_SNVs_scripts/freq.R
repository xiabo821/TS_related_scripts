#!/usr/bin/env Rscript
# Yun Yan (yy1533@nyu.edu)
#
# Goal: count 1) single base freq, 2) tri-nucleotide freq
#


# Modify these 4 lines to specify input and outputs
species_name <- 'human'   #either 'human', or 'mouse'
infpath <- 'YOUR_PATH/human_ensembl90.tsv'
# absolute path. No header, 6 cols, 1-based; Strand should be 1 or -1;  e.g.:
# Chromosome/scaffold name    Gene start (bp) Gene end (bp)   Gene stable ID  Gene name   Strand
# 11	78215297	78418348	ENSG00000033327	GAB2	-1
fout_s <- 'YOUR_PATH/human_ensembl90.genebody.singleFreq.tsv'
# Ourput file of single base frequencies (stranded)
fout_s_unstranded <- 'YOUR_PATH/human_ensembl90.genebody.singleFreq.unstranded.tsv'
# Ourput file of single base frequencies (unstranded)
fout_tri <- 'YOUR_PATH/human_ensembl90.genebody.triFreq.tsv'
# Ourput file of tri-base reference frequencies (stranded)





##### set input genome and required packages #####
suppressPackageStartupMessages({
  library(dplyr); library(magrittr)})

if (species_name == 'human') {
# if using human species
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
suppressPackageStartupMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))
genome_obj <- Hsapiens   
} else if (species_name == 'mouse') {
# if using mouse species
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
genome_obj <- Mmusculus  
} else {
	stop('Error. Species name is not recognisable!')
}


##### Count begin #####

# read in tsv
ginfo <- read.delim(file = infpath, header = F, stringsAsFactors = F)
colnames(ginfo) <- c('chromosome_name', 'start_position', 'end_position',
                     'gene_id', 'gene_name', 'strand_int')
if (is.numeric(ginfo$strand_int)) {
  ginfo$strand_char <- ifelse(ginfo$strand_int > 0, '+', '-')
} else if (('+' %in% ginfo$strand_int) & ('-' %in% ginfo$strand_int)) {
  ginfo$strand_char <- ginfo$strand_int
}

### Speicies-specific information ###
if (species_name == 'mouse') {
#Add 'chr' for consistency with chr_names in BSgenome.Mmusculus.UCSC.mm10 version 3.8, suitable for R3.5
ginfo$chromosome_name = paste0('chr', ginfo$chromosome_name)
#Change 'chrMT' to 'chrM' for consistency with chr_names in BSgenome.Mmusculus.UCSC.mm10 version 3.8, suitable for R3.5
ginfo$chromosome_name <- replace(ginfo$chromosome_name, ginfo$chromosome_name == 'chrMT', 'chrM')
#remove any chr except the chr1-19,X,Y,M
valid_gene <- ginfo$chromosome_name %in% c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                             'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY','chrM')
ginfo <- ginfo[valid_gene,]
#############################
}  else if (species_name == 'human') {
#remove any chr except the chr1-22,X,Y,MT
valid_gene <- ginfo$chromosome_name %in% c('1','2','3','4','5','6','7','8','9','10',
                             '11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT')
ginfo <- ginfo[valid_gene,]
}

# convert to a GenomicRange object
g_gr <- GRanges(
  seqnames = ginfo$chromosome_name,
  ranges = IRanges(ginfo$start_position,
                   ginfo$end_position),
  strand = Rle(strand(ginfo$strand_char)))
mcols(g_gr)$gene_id <- ginfo$gene_id
mcols(g_gr)$gene_name <- ginfo$gene_name
names(g_gr) <- ginfo$gene_id



# Helper functions
alphabetfreq_in_grange <- function(genome_obj, gene_range, gid_abnormal){
  gid_normal <- names(gene_range)
  gene_dna <- getSeq(genome_obj, gene_range) # consider strand info by default
  gene_dna_freq <- alphabetFrequency(gene_dna, as.prob = F)[, c('A', 'C', 'G', 'T')]
  
  out <- data.frame(cbind(
    gene_id = gid_normal,
    width = width(gene_dna),
    gene_dna_freq), stringsAsFactors = F)
  if (!missing(gid_abnormal)){
    tmp <- data.frame(
      gene_id = gid_abnormal,
      width = 0,
      `A` = 0,
      `C` = 0,
      `G` = 0,
      `T` = 0)
    tmp <- tmp[, colnames(out)]
    out <- rbind(out, tmp) %>% dplyr::arrange(., gene_id)
  }
  return(out)
}

trinucleotidefreq_in_grange <- function(genome_obj, gene_range, gid_abnormal){
  gid_normal <- names(gene_range)
  start(gene_range) <- start(gene_range) -1 #extruding 1bp
  end(gene_range) <- end(gene_range) + 1 #extruding 1bp
  gene_dna <- getSeq(genome_obj, gene_range) # consider strand info by default
  
  dna <- getSeq(genome_obj, gene_range)
  dna_cnt <- trinucleotideFrequency(dna, step = 1)
  dna_cnt <- as.data.frame(dna_cnt)
  
  out <- data.frame(cbind(gene_id = gid_normal,
                          width = width(gene_dna)-2, #return extruding 2bps
                          dna_cnt),
                    stringsAsFactors = F)
  
  if (!missing(gid_abnormal)){
    tmp <- data.frame(
      gene_id = gid_abnormal,
      width = 0)
    tmp <- cbind(tmp, as.data.frame(
      matrix(data=0, nrow=length(gid_abnormal), ncol=ncol(out)-2)))
    colnames(tmp) <- colnames(out)
    out <- rbind(out, tmp) %>% dplyr::arrange(., gene_id)
  }
  return(out)
}


#COUNT AND EXPORT
####################################################################
###### Count single base freq
ginfo_base_freq <- alphabetfreq_in_grange(
  genome_obj = genome_obj, gene_range = g_gr)
ginfo_base_freq_unstranded <- ginfo_base_freq;
#group by gene and sum
ginfo_base_freq2 <- cbind(ginfo_base_freq[1],
                          apply(as.matrix(ginfo_base_freq[,-1]), 2, as.numeric))
ginfo_base_freq_aggre <- aggregate(. ~ gene_id, ginfo_base_freq2, sum)


###### Prepare the gene-unstranded counting of base frequency.
# Swap the A-T G-C frequencies of genes on the negative strand.
# ginfo_base_freq_unstranded   ---input file
# ginfo$strand_int   ----strand info
g_is_pos <- ginfo$strand_int == 1; stopifnot(!any(is.na(g_is_pos)))
res_stranded_neg <- ginfo_base_freq_unstranded[!g_is_pos, ]
colnames(res_stranded_neg) <- colnames(ginfo_base_freq_unstranded) 
ginfo_base_freq_unstranded$A[!g_is_pos] <- res_stranded_neg$T
ginfo_base_freq_unstranded$T[!g_is_pos] <- res_stranded_neg$A
ginfo_base_freq_unstranded$G[!g_is_pos] <- res_stranded_neg$C
ginfo_base_freq_unstranded$C[!g_is_pos] <- res_stranded_neg$G
#group by gene and sum
ginfo_base_freq_unstranded2 <- cbind(ginfo_base_freq_unstranded[1],
                          apply(as.matrix(ginfo_base_freq_unstranded[,-1]), 2, as.numeric))
ginfo_base_freq_unstranded_aggre <- aggregate(. ~ gene_id, ginfo_base_freq_unstranded2, sum)





###### Count tri-bases freq
ginfo_tribase_freq <- trinucleotidefreq_in_grange(
  genome_obj = genome_obj, gene_range = g_gr)
#group by gene and sum
ginfo_tribase_freq2 <- cbind(ginfo_tribase_freq[1],
                         apply(as.matrix(ginfo_tribase_freq[,-1]), 2, as.numeric))
ginfo_tribase_freq_aggre <- aggregate(. ~ gene_id, ginfo_tribase_freq2, sum)






# Export
write.table(ginfo_base_freq_aggre,
            file=fout_s,
            quote = F, sep = '\t', row.names = F, col.names = T)
write.table(ginfo_base_freq_unstranded_aggre,
            file=fout_s_unstranded,
            quote = F, sep = '\t', row.names = F, col.names = T)
write.table(ginfo_tribase_freq_aggre,
            file=fout_tri,
            quote = F, sep = '\t', row.names = F, col.names = T)
