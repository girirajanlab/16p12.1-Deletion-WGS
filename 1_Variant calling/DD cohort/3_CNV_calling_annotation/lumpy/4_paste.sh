#!/bin/bash

# Use lumpy to merge VCFs for all samples

# Input and output files
INPUT_PATH=/path/to/input/vcfs # Use the VCFs generated from script 3_regenotype.sh
OUTPUT_PATH=/path/to/output
OUTPUT_NAME=output_prefix

# Reference files
REF=/path/to/hg19/reference.fasta

# Paths to tools
LUMPY=/path/to/lumpy-sv/bin

# Paste calls from all samples into single file
$LUMPY/smoove paste --name $OUTPUT_NAME --outdir $OUTPUT_PATH $INPUT_PATH/*.vcf.gz

