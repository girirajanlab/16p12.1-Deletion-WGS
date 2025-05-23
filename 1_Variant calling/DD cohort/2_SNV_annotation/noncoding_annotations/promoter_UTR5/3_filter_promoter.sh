#!/bin/bash

# Filter for 5' UTR and promoter variants

# Input and output files
INPUT_VCF=/path/to/input.vcf.gz # Use the VCF generated from script 3_annotate_fetal_brain_enhancer.sh
OUTPUT_VCF=/path/to/output.vcf.gz

# Paths to tools
BCFPATH=/path/to/bcftools_v1.12
TABIXPATH=/path/to/tabix

# Filter variants for those within 5'UTRs or upstream of genes
$BCFPATH/bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV19~"UTR5" | INFO/Func.wgEncodeGencodeBasicV19~"upstream"' $INPUT_VCF | bgzip > $OUTPUT_VCF
$TABIXPATH/tabix -p vcf $OUTPUT_VCF
