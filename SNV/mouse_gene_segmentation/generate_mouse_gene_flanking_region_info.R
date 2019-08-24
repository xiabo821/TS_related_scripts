#!/usr/bin/env Rscript
#
#Goal: generate the upstream and downstream bed-like files, using the genebody info as input.
#
#rm(list=ls())


###### set the input filename with path #####
#Note that in all gene_infor files here, the start and end DOES NOT consider strand of genes, just from small to large position.

infpath <- 'YOUR_PATH/mouse_ensembl90.txt'
# absolute path. No header, 6 cols, 1-based; Strand should be 1 or -1;  e.g.:
# Chromosome/scaffold name    Gene start (bp) Gene end (bp)   Gene stable ID  Gene name   Strand
# 11	78215297	78418348	ENSG00000033327	GAB2	-1
###### set parameters #####
flanking_length <- 5000
###### set the outout filenames with path #####
output_upstream <- 'YOUR_PATH/mouse_ensembl90.upstream.tsv'
# Ourput file of gene upstream
output_downstream <- 'YOUR_PATH/mouse_ensembl90.downstream.tsv'
# Ourput file of gene upstream


####################################################################################
# read in input file
ginfo <- read.delim(file = infpath, header = F, stringsAsFactors = F)
colnames(ginfo) <- c('chromosome_name', 'start_position', 'end_position',
                     'gene_id', 'gene_name', 'strand_int')
# setup output file
upstream <- ginfo
downstream <- ginfo
# find indeces of gene strand
minus_strand <- ginfo$strand_int == -1
plus_strand  <- ginfo$strand_int == 1

#Calculate the start and end position of upstream and downstream sequences, respectively.
#Consider the strand. 
#Note that here the start and end DOES NOT consider strand of genes, just from small to large position.
upstream$start_position[minus_strand] <- ginfo$end_position[minus_strand] + 1
upstream$end_position[minus_strand] <- ginfo$end_position[minus_strand] + flanking_length
upstream$start_position[plus_strand] <- ginfo$start_position[plus_strand] - flanking_length
upstream$end_position[plus_strand] <- ginfo$start_position[plus_strand] - 1

downstream$start_position[minus_strand] <- ginfo$start_position[minus_strand] - flanking_length
downstream$end_position[minus_strand] <- ginfo$start_position[minus_strand] - 1
downstream$start_position[plus_strand] <- ginfo$end_position[plus_strand] + 1 
downstream$end_position[plus_strand] <- ginfo$end_position[plus_strand] + flanking_length

#REMOVE gene records which contains minus numbers in start_position and/or end_position 
upstr_drop_list   <- (upstream$start_position < 0 | upstream$end_position < 0)
downstr_drop_list <- (downstream$start_position < 0 | downstream$end_position < 0)
upstream <- upstream[!upstr_drop_list,]
downstream <- downstream[!downstr_drop_list,]




########## Export ##########
# Export results may contain non-printable characters if running from windows machine.
# Check by 'cat -A file_name' in bash/linux command for ^M if containing non-printable characters.
# If yes, remove by 'dos2unix file_name' in bash/linux command.
# Export
write.table(upstream,
            file=output_upstream,
            quote = F, sep = '\t', row.names = F, col.names = F)
# Export
write.table(downstream,
            file=output_downstream,
            quote = F, sep = '\t', row.names = F, col.names = F)