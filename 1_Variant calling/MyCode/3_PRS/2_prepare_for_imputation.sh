#!/bin/bash

# Input and output files
INPUT_LGEN=/path/to/input/array_data # Use the output directory from script 1_signal2legn.py
BLACKLIST=/path/to/long/range/LD/blacklist/from/Price/2008/AJHG.txt # A BED file of regions long-range high LD regions in the genome
OUTPUT_DIR=/path/to/output/plink_files
HRC_CHECK_BIM_SCRIPT=/path/to/HRC-1000G-check-bim.pl # This script can be downloaded from http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
HRC_SITES=/path/to/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz # This file can be downloaded from ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

# Reference files
HG38_HG19_CHAIN=/path/to/hg38ToHg19.over.chain.gz # This file can be downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
REF=/path/to/b37/human_g1k_v37.fasta

# Tools
PLINK=/path/to/PLINK_v1.90
LIFTOVER=/path/to/liftOver
LOPLINK=/path/to/liftOverPlink.py # This scrpit can be downloaded from https://github.com/sritchie73/liftOverPlink/blob/master/liftOverPlink.py
VCF_COOKER=/path/to/vcfCooker
BCFPATH=/path/to/bcftools_v1.12

# Concatenate LGEN files by array
cat $INPUT_LGEN/OMNI/* > $OUTPUT_DIR/OMNI.lgen
cat $INPUT_LGEN/GSA/* > $OUTPUT_DIR/GSA.lgen

# Rename the FAM and MAP files from script 1_signal2lgen.py
cp $INPUT_LGEN/signal_to_lgen.fam $OUTPUT_DIR/OMNI.fam
cp $INPUT_LGEN/signal_to_lgen.fam $OUTPUT_DIR/GSA.fam

cp $INPUT_LGEN/singal_to_lgen_OMNI.map $OUTPUT_DIR/OMNI.map
cp $INPUT_LGEN/singal_to_lgen_GSA.map $OUTPUT_DIR/GSA.map

# Convert LGEN to BED/BIM/FAM
$PLINK --lfile $OUTPUT_DIR/OMNI --make-bed --out $OUTPUT_DIR/OMNI
$PLINK --lfile $OUTPUT_DIR/GSA --make-bed --out $OUTPUT_DIR/GSA

# Merge arrays
echo $OUTPUT_DIR/OMNI > $OUTPUT_DIR/merge.list
echo $OUTPUT_DIR/GSA >> $OUTPUT_DIR/merge.list
$PLINK --merge-list $OUTPUT_DIR/merge.list --out $OUTPUT_DIR/MERGED

# Perform QC filtering
# SNP filters: Cohort MAF <0.05, Hardy-Weinberg equilibrium <1e-6 (under selection), geno<0.01 (SNPs missing in >1% of subjects)
$PLINK --bfile $OUTPUT_DIR/MERGED --maf 0.05 --hwe 1e-6 --geno 0.01 --make-bed --out $OUTPUT_DIR/SNP_QC
# Sample filters: mind<0.01 (Individuals with >1% missing genotypes)
$PLINK --bfile $OUTPUT_DIR/SNP_QC --mind 0.01 --make-bed --out $OUTPUT_DIR/SAMPLE_QC

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

# Lift over PLINK to GRCh37 from GRCh38
$PLINK --bfile $OUTPUT_DIR/DEDUP --recode tab --out $OUTPUT_DIR/TO_LIFT
python2 $LOPLINK -m $OUTPUT_DIR/TO_LIFT.map -o $OUTPUT_DIR/HG19 -c $HG38_HG19_CHAIN -e $LIFTOVER

# Remove unlifted SNPs
cut -f 4 $OUTPUT_DIR/HG19.bed.unlifted | sed "/^#/d" > $OUTPUT_DIR/unlifted_SNPs.list
$PLINK --bfile $OUTPUT_DIR/DEDUP --recode --out $OUTPUT_DIR/EXCLUDE_UNLIFTED --exclude $OUTPUT_DIR/unlifted_SNPs.list
$PLINK --ped $OUTPUT_DIR/EXCLUDE_UNLIFTED.ped --map $OUTPUT_DIR/HG19.map --chr 1-23 --recode --make-bed --out $OUTPUT_DIR/LIFTED

# Create a frequency file
$PLINK --freq --bfile $OUTPUT_DIR/LIFTED --out $OUTPUT_DIR/IMPUTE

# Do pre-imputation QC
perl $HRC_CHECK_BIM_SCRIPT -b $OUTPUT_DIR/DEDUP.bim -f $OUTPUT_DIR/IMPUTE.frq -r $HRC_SITES -h
sh Run-plink.sh # Run the resulting script

# Generate VCF files per chromosome
for i in {1..22}; do $VCF_COOKER vcfCooker --in-bfile $OUTPUT_DIR/DEDUP-updated-chr${i} --ref $REF --out $OUTPUT_DIR/DEDUP-updated-chr${i}.vcf --write-vcf; done
for i in {1..22}; do bgzip $OUTPUT_DIR/DEDUP-updated-chr${i}.vcf; tabix -p vcf $OUTPUT_DIR/DEDUP-updated-chr${i}.vcf.gz; done
for i in {1..22}; do $BCFPATH/bcftools sort $OUTPUT_DIR/DEDUP-updated-chr${i}.vcf.gz -Oz -o $OUTPUT_DIR/DEDUP-updated-chr${i}-sort.vcf.gz; done

# These VCFs are then be uploaded to the TOPMed Imputation Server (https://imputation.biodatacatalyst.nhlbi.nih.gov/#!pages/home)
# The following parameters are used:
# Reference Panel: TOPMed r2
# Input Files: (bgzipped VCFs for all chromosomes)
# Array Build: GRCh37/hg19
# rsq Filter: 0.3
# Phasing: Eagle v2.4 (phased output)
# Population: vs. TOPMed Panel
# Mode: Quality Control & Imputation
