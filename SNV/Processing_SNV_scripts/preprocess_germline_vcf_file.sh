#!/bin/bash 
#SBATCH --job-name=VCF_file_processing # Job name 
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --mail-user=bx347@nyulangone.org 
#SBATCH --ntasks=1 # Run on a single CPU 
#SBATCH --mem=50gb # Job memory request 
#SBATCH --time=10:00:00 # Time limit hrs:min:sec 
#SBATCH --output=VCF_file_processing_%j.log # Standard output and error log
#SBATCH -p cpu_short


input_vcf_file="YOUR_PATH/homo_sapiens"    
#but remove the .vcf of the filename
#e.g. YOUR_PATH/mus_musculus

input_fasta_file="YOUR_PATH/hs_ref_GRCh38.p12_chr"  
#Consider directory. e.g. YOUR_PATH/mm_ref_GRCm38.p4_chr

species="human"  
#e.g. human, mouse 

workdir='YOUR_PATH/human_germline_SNV' 
#absolute path of your working directory, also the results will be kepted there. 
#e.g. YOUR_PATH/mouse_germline_SNV






################################################################################


# your specific codes of bash script for pre-processing the human/mouse germline SNV VCF files:

echo -e "####################################################################################################"
echo -e "Starting time: "
date
cd ${workdir}
# check species.
if [ $species = "human" ]
then
	declare -a chromsome=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" "MT")
	echo -e "Working species: ${species} \n"
elif [ $species = "mouse" ]
then
	declare -a chromsome=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "X" "Y" "MT")
	echo -e "Working species: ${species} \n"
else
	echo -e "ATTENTION!!! Species name not recognized!!! Inputed: ${species} \n"
fi

echo -e "\nInput genome files: $input_fasta_file \n"
echo -e "Input VCF file line length:   "
wc -l ${input_vcf_file}.vcf 
echo -e "\n\n"

#select the 'dbSNP' database & 'SNV' variants
mkdir result_vcf_files
awk '/dbSNP/ && /SNV/' ${input_vcf_file}.vcf > ./result_vcf_files/dbSNP_SNV.vcf
echo -e "dbSNP_derived ${species} SNV recordlength in total:   "
wc -l ./result_vcf_files/dbSNP_SNV.vcf
echo -e "\n\n"



#arr=("a 1" "b" "c")
#for i in "${arr[@]}"; do echo $i; done

echo -e "Chr\tInput_vcf_num\toutput_SNV_num\tReassign_ref_acd_AA(only)\tSwap_VCF_ref&alt(only)\treassign_AA&Correcting_VCF\n" > ${species}_SNV_stat.txt
mkdir processed_SNV_by_chr

for i in "${chromsome[@]}"
do
	echo -e "##########################################"
	echo -e "Starting chr$i at:"
	date
	#awk with '$1 == "XXX"' is perfect for selecting the lines who first elements matches the criteria XXX.
	#Note: awk requires its own varianble reading into 'var', so that to read in the program.
	awk '$1 == var' var="${i}" ./result_vcf_files/dbSNP_SNV.vcf > input.vcf  
	cp ${input_fasta_file}$i.fa input.fasta
	cp input.vcf ./result_vcf_files/chr${i}_dbSNP_SNV.vcf
	#perl processing of VCF file for each chromosome
	perl parse_SNPs.pl > ./processed_SNV_by_chr/${species}_processed_SNV_chr$i.txt

	#generating stats of variants from the processed VCF.
	#chr
	#total_vcf_num
	stat2=$(wc -l < input.vcf)  #count lines and only display the number
	#processed_SNV_num, could be larger than input SNV number because of multiple alternative alleles
	stat3=$(wc -l < ./processed_SNV_by_chr/${species}_processed_SNV_chr$i.txt)
	#fix_AA_only_num
	stat4=$(grep -c -w "reassAA" ./processed_SNV_by_chr/${species}_processed_SNV_chr$i.txt)
	#fix_SNV_only_num
	stat5=$(grep -c -w "fixSNV" ./processed_SNV_by_chr/${species}_processed_SNV_chr$i.txt)
	#fix_AA&SNV_num
	stat6=$(grep -c -w "reassAAfixSNV" ./processed_SNV_by_chr/${species}_processed_SNV_chr$i.txt)
	
	echo -e "chr$i\t$stat2\t$stat3\t$stat4\t$stat5\t$stat6" >> ${species}_SNV_stat.txt

	rm -f input.vcf
	rm -f input.fasta
	echo -e "\nProcessed $stat2 records from chr${i}_dbSNP_SNV.vcf file.\n"
	echo -e "\nGenerated $stat3 unique SNV records and saved in ./processed_SNV_by_chr/${species}_processed_SNV_chr$i.txt.\n"
	
	#Save the output as a combined file which 
	cat ./processed_SNV_by_chr/${species}_processed_SNV_chr$i.txt >> ${species}_processed_SNV_ALL.txt 

	echo -e "Finished chr$i at: "
	date
	echo -e "\n\n"

done
 



echo -e "Finished processing $species VCF file into single SNVs.\n"
echo -e "Finished at : "
date
echo -e "\n"
echo -e "####################################################################################################"
echo -e "\n\n\n\n\n"