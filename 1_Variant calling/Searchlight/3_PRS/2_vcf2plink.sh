#!/bin/bash

# After imputation, convert VCFs to PLINK and perform QC filtering

# Input and output files
IMPUTED_VCFS=/path/to/imputed/vcfs
OUTPUT_DIR=/path/to/output/directory

# Reference files
HG38_HG19_CHAIN=/path/to/hg38ToHg19.over.chain.gz # This file can be downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

# Tools
BCFPATH=/path/to/bcftools_v1.12
PLINK=/path/to/PLINK_v1.90
LIFTOVER=/path/to/liftOver
LOPLINK=/path/to/liftOverPlink.py # This scrpit can be downloaded from https://github.com/sritchie73/liftOverPlink/blob/master/liftOverPlink.py

# Merge VCFs across chromosomes
$BCFPATH/bcftools concat -o $OUTPUT_DIR/TOPMED_impute.vcf.gz -Oz $IMPUTED_VCFS/*dose.vcf.gz

# Convert VCF to PLINK format
$PLINK --vcf $OUTPUT_DIR/TOPMED_impute.vcf.gz --out $OUTPUT_DIR/TOPMED_impute

# Rerun QC steps
# SNP filters: Cohort MAF <0.05, Hardy-Weinberg equilibrium <1e-6 (under selection), geno<0.01 (SNPs missing in >1% of subjects)
$PLINK --bfile $OUTPUT_DIR/TOPMED_impute --maf 0.05 --hwe 1e-6 --geno 0.01 --make-bed --out $OUTPUT_DIR/SNP_QC
# Sample filters: mind<0.01 (Individuals with >1% missing genotypes)
$PLINK --bfile $OUTPUT_DIR/SNP_QC --mind 0.01 --write-snplist --make-bed --out $OUTPUT_DIR/SAMPLE_QC

# Remove samples with high heterozygosity
# Select SNPs in 200 SNP/50 SNP sliding windows that have LD r^2 scores of >0.25
$PLINK --bfile $OUTPUT_DIR/SAMPLE_QC --keep $OUTPUT_DIR/SAMPLE_QC.fam --extract $OUTPUT_DIR/SAMPLE_QC.snplist --indep-pairwise 200 50 0.25 --out $OUTPUT_DIR/SNP_WINDOWS
$PLINK --bfile $OUTPUT_DIR/SNP_WINDOWS --extract $OUTPUT_DIR/SNP_WINDOWS.prune.in --keep $OUTPUT_DIR/SAMPLE_QC.fam --het --out $OUTPUT_DIR/HETEROZYGOSITY

# Use R script heterozygosity_filter.R to identify samples with high heterozygosity
Rscript helper_scripts/heterozygosity_filter.R $OUTPUT_DIR
$PLINK --bfile $OUTPUT_DIR/SAMPLE_QC --extract $OUTPUT_DIR/SNP_WINDOWS.prune.in --keep $OUTPUT_DIR/HETEROZYGOSITY.valid.sample --check-sex --out $OUTPUT_DIR/HET_FILTER

# Remove duplicate SNPs, if any
cut -f 2 $OUTPUT_DIR/HET_FILTER.bim | sort | uniq -d > $OUTPUT_DIR/dup_snps.txt
$PLINK --bfile $OUTPUT_DIR/HET_FILTER --exclude $OUTPUT_DIR/dup_snps.txt --make-bed  --write-snplist --out $OUTPUT_DIR/DEDUP

# Recode PLINK bfiles to make PED/MAP format
$PLINK --bfile $OUTPUT_DIR/DEDUP --recode tab --out $OUTPUT_DIR/RECODE

# Liftover SNPs to hg19 from hg38
python2 $LOPLINK -m $OUTPUT_DIR/RECODE.map -o $OUTPUT_DIR/HG19 -p $OUTPUT_DIR/RECODE.ped -c $HG38_HG19_CHAIN -e $LIFTOVER


