#!/bin/bash

# Annotate variants for those in fetal brain active enhancer regions

# Input and output files
INPUT_VCF=/path/to/input.vcf.gz # Use the VCF generated from script 2_annovar.sh
OUTPUT_VCF=/path/to/output.vcf.gz

# VCFanno needs a TOML file to describe how annotation should be performed
TOML=/path/to/fetal_brain_enhancer.toml

# Paths to tools
VCFANNO_PATH=/path/to/vcfanno
TABIXPATH=/path/to/tabix

# Use VCFanno to annotate variants in fetal brain active enhancer regions
$VCFANNO_PATH/vcfanno_linux64.1 $TOML $INPUT_VCF | bgzip > $OUTPUT_VCF
$TABIXPATH/tabix -p vcf $OUTPUT_VCF

