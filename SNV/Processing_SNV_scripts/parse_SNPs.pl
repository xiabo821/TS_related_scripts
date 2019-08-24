#NOTE: Perl start form 0, wihle the VCF file and the genome starts from 1.

$MAF_cutoff = 0.05;

#open(IN,"Homo_sapiens.GRCh38.dna.chromosome.Y.fa") or die "Error_chr.fasta";
open(IN,"input.fasta") or die "Error_input.fasta";
<IN>; $i = 0;
while(<IN>) {
    $chrom{$i} = $_;
    $i++;}
#calculate the length of each chr DNA sequence row
$length = length($chrom{1})-1;  
open (IN, "input.vcf") or die "Error_input.vcf";
$i = 0;


while (<IN>) {
    @line = split /\t/, $_;
	#Only start with the SNV record which has a position equal or larger than previous ones. 
	#This is to filter out potential contamination from the SNV records of other chromasomes.
	#if $location > $line[1], meaning the coming SNV has a position smaller that the previous one, 
	#which indicates SNV of a new chr (VCF is sorted by chromosome position).
    if ($location <= $line[1]) {
	$i++;
	$location   = $line[1];
	$rsID       = $line[2];
	$ref_allele = $line[3];
	$alt_allele = $line[4];
	$info	    = $line[7];
	
	## 	Filter out the records with MAF > $MAF_cutoff
	@info_array = split /;/, $info;
	chomp(@info_array);
	@matches = grep /MAF/, @info_array;
	@MAF_array = split /=/, $matches[0];
	if ($MAF_array[1] > $MAF_cutoff){
		next;
	}
	
	
	## Calculate the chr_DNA rows and positions.
	$chrom_line     = int($location / $length);
	$chrom_position = ($location % $length);
	@DNA = split //, $chrom{$chrom_line}; #Capture the specific chr DNA sequence row
	#In most cases, just capture the sequences.  Perl start from 0.
	$before = $DNA[$chrom_position-2];
	$ref    = $DNA[$chrom_position-1];
	$after  = $DNA[$chrom_position]; 
	#Special: head and tail base
	if ($chrom_position == 1) {
		#Correct the 5' base capture
		@DNA_above = split //, $chrom{$chrom_line-1} ;
	    $before = $DNA_above[$length-1]; #Perl start from 0. Last base.
	}
	elsif ($chrom_position == 0) {
		#Correct the chr DNA sequence row, followed by correcting the 5'/ref/3' base capture
	    $chrom_position = $length; #indicate the last base of each DNA sequence row
		$chrom_line = $chrom_line - 1;   #the remainder == 0 indicates the last base of the above-row
		@DNA = split //, $chrom{$chrom_line} ; #Redefine the DNA variable to capture the right row of DNA sequence.
		
		$before = $DNA[$chrom_position-2]; #Perl start from 0.
		$ref    = $DNA[$chrom_position-1]; #Perl start from 0.
		
		@DNA_below = split //, $chrom{$chrom_line+1} ;
	    $after = $DNA_below[0]; #Perl start from 0.
	}
	
	
	## Fix genome (input.fasta) reference allele if the AA(Acenstral Allele) is different from the reference allele.
	$mark = "";
	@AA = grep /AA=/, @info_array;
	@AA_array = split /=/, $AA[0];
	$size  = @AA_array;
	if 	($size == 2) {
		if ($AA_array[1] !~ $ref){
			$mark = "reassAA";  # revAA indicates the inconsistent VCF reference to the AA(Acenstral Allele). 
			$ref  = $AA_array[1];
		}
	}


	## Check whether there is inconsistency in the .VCF file. Possibly have multiple_alternatives
	@multi_alt = split /,/, $alt_allele;
	$size  = @multi_alt; #Will be used in the end for printing out results.
	#check the VCF file ref_allele to the genome reference (.fasta)
	if ($ref_allele !~ $ref){
		$mark .= "fixSNV";  # Append fixSNV, indicating the inconsistent VCF reference to the genome sequence (AA reassigned). Then probably swap the ref and alt alleles. 
		if (grep($ref, @multi_alt)) {
			for(@multi_alt){s/$ref/$ref_allele/g}; #These bases are muturally exclusive. Then the $ref can be replaced by $ref_allele
			$ref_allele = $ref;	
			} 
		else {
			next;   #Could not identify reference allele, thus drop this record and move on to the next SNV.
			}
		}
	
	#export mutation types
	$tri_ref = "$before$ref$after";

		
	## Print out the variant records. Consider multiple variant alleles.
	if    ($size == 1){
		print "$line[0]\t$location\t$rsID\t$ref_allele\t$tri_ref\t$multi_alt[0]\t$matches[0]\t$AA[0]\t$mark\n";
	}
	elsif ($size == 2) {
		print "$line[0]\t$location\t$rsID\t$ref_allele\t$tri_ref\t$multi_alt[0]\t$matches[0]\t$AA[0]\t$mark\n";
		print "$line[0]\t$location\t$rsID\t$ref_allele\t$tri_ref\t$multi_alt[1]\t$matches[0]\t$AA[0]\t$mark\n";
	}
	elsif ($size == 3) {
		print "$line[0]\t$location\t$rsID\t$ref_allele\t$tri_ref\t$multi_alt[0]\t$matches[0]\t$AA[0]\t$mark\n";
		print "$line[0]\t$location\t$rsID\t$ref_allele\t$tri_ref\t$multi_alt[1]\t$matches[0]\t$AA[0]\t$mark\n";
		print "$line[0]\t$location\t$rsID\t$ref_allele\t$tri_ref\t$multi_alt[2]\t$matches[0]\t$AA[0]\t$mark\n";
	}
	
    }
	else {
		last; #if $location > $line[1], break the while loop and move on. See the begining for explanation.
		}  

}
