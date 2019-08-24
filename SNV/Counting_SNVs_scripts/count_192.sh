#!/bin/bash
#SBATCH --job-name=SNV_count # Job name 
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --mail-user=bx347@nyulangone.org 
#SBATCH --ntasks=1 # Run on a single CPU 
#SBATCH --mem=50gb # Job memory request 
#SBATCH --time=10:00:00 # Time limit hrs:min:sec 
#SBATCH --output=SNV_count_%j.log # Standard output and error log
#SBATCH -p cpu_short

# Yun Yan (yy1533@nyu.edu)

set -e

module load bedops/2.4.35
module load r/3.5.1


run_prefix=''
#your prefered run name prefix.

workdir='YOUR_PATH/gene_body' 
#absolute path

infile_gene='YOUR_PATH/human_ensembl90.tsv' 
# absolute path. No header, 6 cols; 1-based; Strand should be 1 or -1; e.g.:
# Chromosome/scaffold name    Gene start (bp) Gene end (bp)   Gene stable ID  Gene name   Strand
# 11	78215297	78418348	ENSG00000033327	GAB2	-1

infile_vcf='YOUR_PATH/human_processed_SNV_ALL.txt' 
#absolute path to the clean_SNV files; No header, 9 cols, 1-based. e.g.: 
# 10    68570385    rs191485190 A   CAT G   MAF=0.00878594  AA=A    reassAAfixSNV





########## Split the SNV type and sort ##########
echo "Part-0 1-based convert to 0-based"; 
date
mkdir -p ${workdir}/database
cd ${workdir}/database
awk '{$2=$2-1; print}' OFS='\t' $infile_gene | sort-bed --unique - > ${run_prefix}_gr.bed
awk '{print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' $infile_vcf > vcf.bed
cd ${workdir}; date; echo

echo "Part-1.1 Split 0-based VCF to 192 files"; date
mkdir -p ${workdir}/database/vcf192
cd ${workdir}/database/vcf192
awk ' {name=$6$7".bed"; print $0 >> name} ' OFS='\t' ${workdir}/database/vcf.bed
cd ${workdir}

echo "Part-1.2 For each subseted VCF file, sort."
for f in ${workdir}/database/vcf192/*.bed; do
    sort-bed --unique $f > $f.sorted
    rm $f
done
cd ${workdir}; date; echo






########## counting begin ##########
echo "Part-2 For each subseted VCF file, count upon genes."; date
mkdir -p ${workdir}/count/vcf192
cd ${workdir}/count/vcf192
for f in ${workdir}/database/vcf192/*.bed.sorted; do
    fb=$(basename $f) #AAAG.bed.sorted
    type_name=${fb%%.*}; echo 'counting' ${type_name}
    bedmap --echo --count --delim "\t" \
        ${workdir}/database/${run_prefix}_gr.bed $f > ${type_name}
done
cd ${workdir}; date; echo

echo "Part-3 Merge all count files and fix strand"; date
cd ${workdir}
Rscript merge_192.R \
    $workdir \
    ${run_prefix}_result.192.unstrand.tsv \
    ${run_prefix}_result.192.stranded.tsv
echo "Finished"; 
date