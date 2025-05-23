#!/bin/bash

# Convert PLINK files into VCF for peddy

# Input and output files
INPUT_PLINK=/path/to/input/plink/file # Use the final output of script 3_vcf2plink.sh
OUTPUT_DIR=/path/to/output/directory

# Reference files
PEDDY_SITES=/path/to/peddy/GRCH37.sites # These sites are available on GitHub: https://github.com/brentp/peddy/blob/master/peddy/GRCH37.sites

# Tools
PLINK=/path/to/PLINK_v1.90
BCFPATH=/path/to/bcftools_v1.12

# peddy will need to be available in your environment
# Instructions for installing peddy are available from: https://github.com/brentp/peddy

# Convert peddy sites to a format usable by PLINK
sed 's/:/\t/g' $PEDDY_SITES | cut -f 1-2 > tmp1
sed 's/:/\t/g' $PEDDY_SITES | cut -f 2 > tmp2
paste tmp1 tmp2 $PEDDY_SITES > $OUTPUT_DIR/peddy_sites.list
rm tmp1
rm tmp2

# Subset PLINK to only peddy sites
$PLINK --bfile $INPUT_PLINK --extract range $OUTPUT_DIR/peddy_sites.list --recode vcf bgz --out $OUTPUT_DIR/PEDDY_SUBSET

# Sort and index VCF
$BCFPATH/bcftools sort $OUTPUT_DIR/PEDDY_SUBSET.vcf.gz | bgzip > $OUTPUT_DIR/PEDDY_SUBSET_SORTED.vcf.gz
tabix -p vcf $OUTPUT_DIR/PEDDY_SUBSET_SORTED.vcf.gz

# Peddy requires a read-depth field in the VCF (not present because VCFs were made from BED files)
# Populate the VCF with dummy DP values
gunzip -c $OUTPUT_DIR/PEDDY_SUBSET_SORTED.vcf.gz | sed 's/0\/0/20/g' | sed 's/1\/1/20/g' | sed 's/0\/1/20/g' | sed 's/GT/DP/g' | sed 's/Number=1/Number=R/g' | sed 's/Type=String/Type=Integer/g' | sed 's/Genotype/Read Depth/g' | $BCFPATH/bcftools sort | bgzip > $OUTPUT_DIR/DUMMY_DP.vcf.gz
tabix -p vcf $OUTPUT_DIR/DUMMY_DP.vcf.gz

$BCFPATH/bcftools annotate -a $OUTPUT_DIR/DUMMY_DP.vcf.gz -c FORMAT/DP $OUTPUT_DIR/PEDDY_SUBSET_SORTED.vcf.gz | bgzip > $OUTPUT_DIR/DP_ADDED.vcf.gz
tabix -p vcf $OUTPUT_DIR/DP_ADDED.vcf.gz

# Create a PED file for peddy
$PLINK --bfile $INPUT_PLINK --recode tab --out $OUTPUT_DIR/RECODE

# Run peddy
python -m peddy --plot --prefix $OUTPUT_DIR/peddy_output $OUTPUT_DIR/DP_ADDED.vcf.gz $OUTPUT_DIR/RECODE.ped

