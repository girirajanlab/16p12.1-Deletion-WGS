#!/bin/bash

# Input and output files
INPUT_PLINK_PATH=/path/to/input/UKB/BED/BIM/FAM/files
INPUT_FAM=/path/to/input/fam/file.fam
SAMPLES_16P=/path/to/16p12.1/deletion/samples.fam # A list of samples of interest in FAM format
BLACKLIST=/path/to/long/range/LD/blacklist/from/Price/2008/AJHG.txt # A BED file of regions long-range high LD regions in the genome
OUTPUT_DIR=/path/to/output/plink_files
HRC_CHECK_BIM_SCRIPT=/path/to/HRC-1000G-check-bim.pl # This script can be downloaded from http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
HRC_SITES=/path/to/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz # This file can be downloaded from ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

# Reference files
REF=/path/to/b37/human_g1k_v37.fasta

# Tools
PLINK=/path/to/PLINK_v1.90
VCF_COOKER=/path/to/vcfCooker
BCFPATH=/path/to/bcftools_v1.12

# Subset out desired 16p12.1 deletion samples
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y XY MT
do
	cp $INPUT_PLINK_PATH/ukb_snp_chr${chrom}_v2.bim $INPUT_PLINK_PATH/ukb_cal_chr${chrom}_v2.bim
	cp $INPUT_FAM $INPUT_PLINK_PATH/ukb_cal_chr${chrom}_v2.fam

	$PLINK --bfile $INPUT_PLINK_PATH/ukb_cal_chr${chrom}_v2 --keep $SAMPLES_16P  --make-bed --out $OUTPUT_DIR/SAMPLES_16P_chr${chrom}

	rm $INPUT_PLINK_PATH/ukb_cal_chr${chrom}_v2.bim
	rm $INPUT_PLINK_PATH/ukb_cal_chr${chrom}_v2.fam

	echo $OUTPUT_DIR/SAMPLES_16P_chr${chrom} >> $OUTPUT_DIR/plink_merge.list
done

# Merge PLINK per chromosome
$PLINK --merge-list $OUTPUT_DIR/plink_merge.list --out $OUTPUT_DIR/CHR_MERGED

# Perform QC filtering
# SNP filters: Cohort MAF <0.05, Hardy-Weinberg equilibrium <1e-5 (under selection), geno<0.05 (SNPs missing in >5% of subjects)
$PLINK --bfile $OUTPUT_DIR/SAMPLES_16P --maf 0.05 --hwe 1e-5 --geno 0.05 --chr 1-23 --make-bed --out $OUTPUT_DIR/SNP_QC
# Sample filters: mind<0.01 (Individuals with >1% missing genotypes)
$PLINK --bfile $OUTPUT_DIR/SNP_QC --mind 0.01 --me 0.05 0.05 --make-bed --out $OUTPUT_DIR/SAMPLE_QC

# Remove SNPs within long-range LD regions (from Price et al, AJHG 2002)
$PLINK --bfile $OUTPUT_DIR/SAMPLE_QC --exclude 'range' $BLACKLIST --write-snplist --make-bed --out $OUTPUT_DIR/LR_LD

# Remove samples with high heterozygosity
# Select SNPs in 200 SNP/50 SNP sliding windows that have LD r^2 scores of >0.25
$PLINK --bfile $OUTPUT_DIR/LR_LD --keep $OUTPUT_DIR/LR_LD.fam --extract $$OUTPUT_DIR/LR_LD.snplist --indep-pairwise 200 50 0.25 --out $OUTPUT_DIR/SNP_WINDOWS
$PLINK --bfile $OUTPUT_DIR/SNP_WINDOWS --extract $OUTPUT_DIR/SNP_WINDOWS.prune.in --keep $OUTPUT_DIR/LR_LD.fam --het --out $OUTPUT_DIR/HETEROZYGOSITY

# Use R script heterozygosity_filter.R to identify samples with high heterozygosity
Rscript helper_scripts/heterozygosity_filter.R $OUTPUT_DIR
$PLINK --bfile $OUTPUT_DIR/LR_LD --extract $OUTPUT_DIR/SNP_WINDOWS.prune.in --keep $OUTPUT_DIR/HETEROZYGOSITY.valid.sample --check-sex --out $OUTPUT_DIR/HET_FILTER

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
