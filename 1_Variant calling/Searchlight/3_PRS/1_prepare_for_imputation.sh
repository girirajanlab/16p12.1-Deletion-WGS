#!/bin/bash

# Input and output files
INPUT_PLINK=/path/to/input/array/plink/files
# Input array data can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
# Note that the downloaded array data is separated into three batches - arrays should be merged before PRS calculation
OUTPUT_DIR=/path/to/output/plink_files
HRC_CHECK_BIM_SCRIPT=/path/to/HRC-1000G-check-bim.pl # This script can be downloaded from http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
HRC_SITES=/path/to/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz # This file can be downloaded from ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

# Reference files
REF=/path/to/b37/human_g1k_v37.fasta

# Tools
PLINK=/path/to/PLINK_v1.90
VCF_COOKER=/path/to/vcfCooker
BCFPATH=/path/to/bcftools_v1.12

# Perform QC filtering
# SNP filters: Cohort MAF <0.05, Hardy-Weinberg equilibrium <1e-6 (under selection), geno<0.01 (SNPs missing in >1% of subjects)
$PLINK --bfile $INPUT_PLINK --maf 0.05 --hwe 1e-6 --geno 0.01 --make-bed --out $OUTPUT_DIR/SNP_QC
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

# Create a frequency file
$PLINK --freq --bfile $OUTPUT_DIR/DEDUP --out $OUTPUT_DIR/IMPUTE

# Do pre-imputation QC
perl $HRC_CHECK_BIM_SCRIPT -b $OUTPUT_DIR/DEDUP.bim -f $OUTPUT_DIR/IMPUTE.frq -r $HRC_SITES -h
sh Run-plink.sh # Run the resulting script

# Generate VCF files per chromosome
for i in {1..23}; do $VCF_COOKER vcfCooker --in-bfile $OUTPUT_DIR/DEDUP-updated-chr${i} --ref $REF --out $OUTPUT_DIR/DEDUP-updated-chr${i}.vcf --write-vcf; done
for i in {1..23}; do bgzip $OUTPUT_DIR/DEDUP-updated-chr${i}.vcf; tabix -p vcf $OUTPUT_DIR/DEDUP-updated-chr${i}.vcf.gz; done
for i in {1..23}; do $BCFPATH/bcftools sort $OUTPUT_DIR/DEDUP-updated-chr${i}.vcf.gz -Oz -o $OUTPUT_DIR/DEDUP-updated-chr${i}-sort.vcf.gz; done

# These VCFs are then be uploaded to the TOPMed Imputation Server (https://imputation.biodatacatalyst.nhlbi.nih.gov/#!pages/home)
# The following parameters are used:
# Reference Panel: TOPMed r2
# Input Files: (bgzipped VCFs for all chromosomes)
# Array Build: GRCh37/hg19
# rsq Filter: 0.3
# Phasing: Eagle v2.4 (phased output)
# Population: vs. TOPMed Panel
# Mode: Quality Control & Imputation
