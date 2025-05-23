#!/bin/bash

# Merge STR calls across families using mergeSTR and apply population-level filters using dumpSTR

# Input and output files
AUTO_VCF="/path/to/vcfs/for/family,/as/a/comma/separated/list" # Use the autosomal VCFs for all families generated from script 1_run_GangSTR_dumpSTR.sh here
CHRX_VCF="/path/to/vcfs/for/family,/as/a/comma/separated/list" # Use the chrX VCFs for all families generated from script 1_run_GangSTR_dumpSTR.sh here
FEMALE_SAMPLES="/list/of/female/samples,/as/a/comma/separated/list"
OUTPUT_DIR=/path/to/output/directory

# Reference files
SUPDUP=/path/to/GRCh37/genomicSuperDups.bed
# Note: this file was obtained from UCSC Table Browser
# with the options: Mammal, HUman, GRCh37/hg19, Repeats, Segemntal Dups, genomicSuperDups, output format BED
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1047622799_NYRO2yV9XWpcQdLYnWOEDqBIAVUJ&clade=mammal&org=&db=hg19&hgta_group=rep&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=


# Tools
MERGESTR=/path/to/mergeSTR
DUMPSTR=/path/to/dumpSTR
BCFPATH=/path/to/bcftools_v1.3.1

# Merge STR calls across families
mergeSTR \
	--vcfs $AUTO_VCF \
	--out $OUTPUT_DIR/autosome

mergeSTR \
	--vcfs $CHRX_VCF \
	--out $OUTPUT_DIR/chrX

# For chrX calls, make a copy of the VCF with only female samples
# This file will be used for chrX Hardy-Weinburg equilibrium calculations
$BCFPATH/bcftools view -s $FEMALE_SAMPLES $OUTPUT_DIR/chrX > $OUTPUT_DIR/chrX_Female_only.vcf

# Filter calls with dumpSTR for
# 1. Sites with a call rate >0.8
# 2. Remove loci overlapping segmental duplication regions
# 3. Departure from Hardy Weinberg equilibrium >0.00001
$DUMPSTR \
	--vcf $OUTPUT_DIR/autosome \
	--out $OUTPUT_DIR/autosome_filter \
	--min-locus-hwep 0.00001 \
	--min-locus-callrate 0.8 \
	--filter-regions $SUPDUP  \
	--filter-regions-names SEGDUP

$DUMPSTR \
	--vcf $OUTPUT_DIR/chrX \
	--out $OUTPUT_DIR/chrX_filter \
	--min-locus-callrate 0.8 \
	--filter-regions $SUPDUP  \
	--filter-regions-names SEGDUP

# For chrX, Hardy Weinberg equilibrium p value is based only on female samples
$DUMPSTR \
	--vcf $OUTPUT_DIR/chrX_Female_only.vcf \
	--min-locus-hwep 0.00001  \
	--out $OUTPUT_DIR/chrX_filter_hw

