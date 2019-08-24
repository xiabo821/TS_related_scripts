# File structure
Put all the scripts and inputs inside a working directory, e.g., '/Users/yunyan/Downloads/icm-count/'. 

```
/Users/yunyan/Downloads/icm-count/
├── TET1_AG.tsv           # input.1
├── count_192.sh          # Run this script to count 192 snps occurence
├── count_12.sh           # Run this script to count 12 snps occurence
├── human_ensembl90.tsv   # input 2
├── merge_12.R            # helper function; don't run
└── merge_192.R           # helper function; don't run
```

- TET1_AG.tsv: example of SNP with contexts created by "Processing SNVs" codes .
- human_ensembl90.tsv: download from Ensembl 90. It is a 6-column file: Chromosome/scaffold name, Gene start (bp), Gene end (bp), Gene stable ID, Gene name, Strand.

# Dependencies

Tool: [bedops](https://bedops.readthedocs.io/en/latest/content/installation.html)

R packages: [gtools](https://cran.r-project.org/web/packages/gtools/index.html), [readr](https://readr.tidyverse.org/), dplyr, magrittr

# Run

In the 'count.sh' Bash script file, change the working directory and input file path accordingly.

```bash
workdir='/Users/yunyan/Downloads/icm-count' #absolute path

infile_gene='/Users/yunyan/Downloads/icm-count/human_ensembl90.tsv' #absolute path

infile_bo_vcf='/Users/yunyan/Downloads/icm-count/TET1_AG.tsv' #absolute path
```

Then go to the specified working directory and run the 'count.sh'.
```
cd /Users/yunyan/Downloads/icm-count
sh count_192.sh  # Run this script to count 192 snps occurence
sh count_12.sh   # Run this script to count 12 snps occurence
```

Finally, you will have the results as below.

```
./
├── TET1_AG.tsv
├── count
│   ├── vcf12
│   │   └── AG
│   └── vcf192
│       ├── AAAG
│       └── ...
├── count_12.sh
├── count_192.sh
├── database
│   ├── gr.bed
│   ├── vcf.bed
│   ├── vcf12
│   │   └── AG.bed.sorted
│   └── vcf192
│       ├── AAAG.bed.sorted
│       └── ...
├── human_ensembl90.tsv
├── merge_12.R
├── merge_192.R
├── result.12.stranded.tsv
├── result.12.unstrand.tsv
├── result.192.stranded.tsv
└── result.192.unstrand.tsv

```