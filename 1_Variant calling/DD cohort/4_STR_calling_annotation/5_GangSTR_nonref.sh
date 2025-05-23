#!/bin/bash

# Call variants again on non-reference loci identified for a family
# Input and output files
INPUT_BAMS="/path/to/bams/for/family.bam,/as/a/comma/separated/list.bam" # Use the BAMs generated in 1_GATK/1_GATK_pipeilne.sh here
SEXES="comma,separated,list,of,sample,sexes" # Use 1/2 chrX ploidy encoding where 1 indicates male and 2 indicates female
AUTO_BED=/path/to/input/autosome/regions.bed # Use the autosome BED generated from script 4_filter_calls.sh here
CHRX_BED=/path/to/input/autosome/regions.bed # Use the chrX BED  generated from script 4_filter_calls.sh here
OUTPUT_DIR=/path/to/output/directory

# Reference files
REF=/path/to/hg19_reference

# Tools
GANGSTR=/path/to/GangSTR-2.5
BCFPATH=/path/to/bcftools_v1.3.1

# Call STRs for non-reference loci using GangSTR
# Autosomes
$GANGSTR/build/GangSTR \
	--bam $INPUT_BAMS \
	--ref $REF \
	--regions $AUTO_BED \
	--out $OUTPUT_DIR/autosome \
	--include-ggl 
# ChrX
/data5/software/GangSTR-2.5/build/GangSTR \
	--bam $INPUT_BAMS \
	--ref $REF \
	--regions $CHRX_BED  \
	--out $OUTPUT_DIR/chrX \
	--include-ggl \
	--chrom chrX \
	--samp-sex $SEXES

# Zip and index
bgzip $OUTPUT_DIR/autosome.vcf
tabix -p vcf $OUTPUT_DIR/autosome.vcf.gz

bgzip $OUTPUT_DIR/chrX.vcf
tabix -p vcf $OUTPUT_DIR/chrX.vcf.gz

# Remove sex chromosomes from autosome file
$BCFPATH/bcftools view -t ^chrX,chrY $OUTPUT_DIR/autosome.vcf.gz | bgzip > $OUTPUT_DIR/autosome_nosex.vcf.gz

# Concat autosomal and sex chromsome files
$BCFPATH/bcftools concat $OUTPUT_DIR/autosome_nosex.vcf.gz $OUTPUT_DIR/chrX.vcf.gz | bgzip > $OUTPUT_DIR/str_calls.vcf.gz
tabix -p vcf $OUTPUT_DIR/str_calls.vcf.gz

