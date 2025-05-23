#!/bin/bash

# Call SVs with CNVnator

# Input and output files
INPUT_BAM=/path/to/input.bam # Use the sorted BAM files generated from bwa-mem in script 1_GATK_pipeline.sh
OUTPUT_VCF_PREFIX=/path/to/output_prefix
TMP_DIR=/path/to/tmp/directory
PREFIX=intermediate_file_prefix

# Reference files
REF_DIR=/path/to/hg19/individual/chromosome/reference/fastas
CHROMS=/path/to/chromosome_bedfile.bed.gz

# Paths to tools
CNVNATOR_DIR=/path/to/CNVnator_v0.4.1/
BCFPATH=/path/to/bcftools_v1.12

# Call SVs using CNVnator
$CNVNATOR_DIR/src/cnvnator -root $TMP_DIR/$PREFIX.root -tree $INPUT_BAM -chrom $(seq -f 'chr%g' 1 22) chrX chrY chrM
$CNVNATOR_DIR/src/cnvnator -root $TMP_DIR/$PREFIX.root -his 200 -chrom $(seq -f 'chr%g' 1 22) chrX chrY chrM -d $REF_DIR
$CNVNATOR_DIR/src/cnvnator -root $TMP_DIR/$PREFIX.root -stat 200
$CNVNATOR_DIR/src/cnvnator -root $TMP_DIR/$PREFIX.root -partition 200
$CNVNATOR_DIR/src/cnvnator -root $TMP_DIR/$PREFIX.root -call 200 > $TMP_DIR/$PREFIX.cnvnator.out

# Export CNV calls into VCF format
$CNVNATOR_DIR/src/cnvnator2VCF.pl -prefix $PREFIX -reference hg19 $TMP_DIR/$PREFIX.cnvnator.out $REF_DIR > $OUTPUT_VCF_PREFIX.vcf

# Remove any CNVs with a natorQ0 value >= 0.5
$BCFPATH/bcftools filter -i "INFO/natorQ0 < 0.5" $OUTPUT_VCF_PREFIX.vcf > $OUTPUT_VCF_PREFIX.filtered.vcf

rm -r $TMP_DIR
