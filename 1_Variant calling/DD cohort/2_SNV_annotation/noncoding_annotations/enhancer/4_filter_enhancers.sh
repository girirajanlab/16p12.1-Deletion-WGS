#!/bin/bash

# Input and output files
INPUT_VCF=/path/to/input.vcf.gz # Use the VCF generated from script 3_annotate_fetal_brain_enhancer.sh
OUTPUT_VCF=/path/to/output.vcf.gz

# Paths to tools
BCFPATH=/path/to/bcftools_v1.12
TABIXPATH=/path/to/tabix

$BCFPATH/bcftools view -i 'INFO/fb_enhancer~"putative_fetal_brain_enhancer"' $INPUT_VCF | bgzip > $OUTPUT_VCF
$TABIXPATH/tabix -p vcf $OUTPUT_VCF
