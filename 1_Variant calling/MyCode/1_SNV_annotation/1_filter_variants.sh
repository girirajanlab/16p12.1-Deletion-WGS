#!/bin/bash

# Perform QC on MyCode VCFs
# Filters/updates:
# 1. Non-ref alleles
# 2. Split and normalize left
# 3. QUAL>=50
# 4. Read depth>=8
# 5. Alternative allele depth>0
# 6. Allele balance >=0.25 and <=0.75 or >=0.9
# 7. QUAL/(Alternative allele depth)>=1.5

# Input and output files
FASTA=/path/to/hg38/reference.fasta
INPUT_VCF=/path/to/input.vcf.gz # Input VCFs were provided by MyCode
OUTPUT_FILE=/path/to/output.vcf.gz

# Paths to tools
BCFPATH=/path/to/bcftools_v1.12

# Perform initital filtering
$BCFPATH/bcftools norm -m-both $INPUT_VCF | bcftools view -i 'GT="alt" & QUAL>=50 & FORMAT/DP>=8 & FORMAT/AD[:1]>0 & (FORMAT/AD[:1])/(FORMAT/DP)>=0.25 & QUAL/(FORMAT/AD[:1])>=1.5' | bcftools view -i '(FORMAT/AD[:1])/(FORMAT/DP)<=0.75 | (FORMAT/AD[:1])/(FORMAT/DP)>=0.9' |  bcftools norm -f $FASTA | bgzip > $OUTPUT_FILE
tabix -p vcf $OUTPUT_FILE

