#!/bin/bash

# Annotate gnomad ancestry-specific allele frequencies

# Input and output files
INPUT_VCF=/path/to/input.vcf.gz # Use the VCF generated from script 3_filter_exonic.sh here
OUTPUT_FILE=/path/to/output.vcf.gz

# VCFanno needs a TOML file to describe how annotation should be performed
TOML=/path/to/gnomad_annotation_file.toml

# Paths to tools
VCFANNO_PATH=/path/to/vcfanno
TABIXPATH=/path/to/tabix

# Annotate allele frequencies using VCFanno
$VCFANNO_PATH/vcfanno_linux64.1 $TOML $INPUT_VCF | bgzip > $OUTPUT_FILE
$TABIXPATH/tabix -p vcf $OUTPUT_FILE

