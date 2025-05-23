#!/bin/bash

# Input and output files
INPUT_VCF=/path/to/input.vcf.gz # Use the VCF generated from script 8_annotate_cadd.sh here
OUTPUT_FILE=/path/to/output.vcf.gz

# VCFanno needs a TOML file to describe how annotation should be performed
TOML=/path/to/clinvar_annotation_file.toml

# Paths to tools
VCFANNO_PATH=/path/to/vcfanno
TABIXPATH=/path/to/tabix

# Annotate ClinVar alleles using VCFanno
$VCFANNO_PATH/vcfanno_linux64.1 $TOML $INPUT_VCF | bgzip > $OUTPUT_FILE
$TABIXPATH/tabix -p vcf $OUTPUT_FILE
