#!/bin/bash

# Use GangSTR to call STR variants on a single family and dumpSTR to filter STR calls

# Input and output files
INPUT_BAMS="/path/to/bams/for/family.bam,/as/a/comma/separated/list.bam" # Use the BAMs generated in 1_GATK/1_GATK_pipeilne.sh here
SEXES="comma,separated,list,of,sample,sexes" # Use 1/2 chrX ploidy encoding where 1 indicates male and 2 indicates female
OUTPUT_DIR=/path/to/output/directory

# Reference files
REF=/path/to/hg19_reference
GANGSTR_REF=/path/to/hg19/gangstr/reference/hg19_ver13.1.bed  # This file can be downloaded here: https://s3.amazonaws.com/gangstr/hg19/genomewide/hg19_ver13_1.bed.gz

# Tools
GANGSTR=/path/to/GangSTR-2.4
DUMPSTR=/path/to/dumpSTR

# Call STRs with GangSTR
# Autosomes
$GANGSTR/build/GangSTR \
	--bam $INPUT_BAMS \
	--ref $REF \
	--regions $GANGSTR_REF \
	--out $OUTPUT_DIR/GangSTR \
	--include-ggl
# chrX
$GANGSTR/build/GangSTR \
	--bam $INPUT_BAMS \
	--ref $REF \
	--regions $GANGSTR_REF \
	--out $OUTPUT_DIR/GangSTR_chrX \
	--include-ggl \
	--chrom chrX \
	--ploidy $SEXES

# Compress and index
bgzip $OUTPUT_DIR/GangSTR.vcf
tabix -p vcf $OUTPUT_DIR/GangSTR.vcf.gz
bgzip $OUTPUT_DIR/GangSTR_chrX.vcf
tabix -p vcf $OUTPUT_DIR/GangSTR_chrX.vcf.gz

# Filter STR calls with dumpSTR
$DUMPSTR \
	--verbose \
	--vcf $OUTPUT_DIR/GangSTR.vcf.gz \
	--out $OUTPUT_DIR/dumpSTR \
	--min-call-DP 20 \
	--max-call-DP 1000 \
	--filter-spanbound-only \
	--filter-badCI \
	--readlen 150

$DUMPSTR \
	--verbose \
	--vcf $OUTPUT_DIR/GangSTR_chrX.vcf.gz \
	--out $OUTPUT_DIR/dumpSTR_chrX \
	--min-call-DP 10 \
	--max-call-DP 1000 \
	--filter-spanbound-only \
	--filter-badCI \
	--readlen 150

# Compress and index
bgzip $OUTPUT_DIR/dumpSTR.vcf
tabix -p vcf $OUTPUT_DIR/dumpSTR.vcf.gz
bgzip $OUTPUT_DIR/dumpSTR_chrX.vcf
tabix -p vcf $OUTPUT_DIR/dumpSTR_chrX.vcf.gz


