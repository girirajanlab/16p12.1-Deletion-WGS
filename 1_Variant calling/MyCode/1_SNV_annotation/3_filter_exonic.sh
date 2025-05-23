#!/bin/bash

# Filter variants for exonic or splicing variants

# Paths to tools
BCFPATH=/path/to/bcftools_v1.12
TABIXPATH=/path/to/tabix

# Input and output files
INPUT_VCF=/path/to/input.hg38_multianno.vcf # Use the output of script 2_annovar.sh here
OUTPUT_FILE=/path/to/output.vcf.gz

# Filter variants for exonic or splicing variants
$BCFPATH/bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV38~"exonic" | INFO/Func.wgEncodeGencodeBasicV38~"splicing"' $INPUT_VCF | bgzip > $OUTPUT_FILE
$TABIXPATH/tabix -p vcf $OUTPUT_FILE

