#!/bin/bash

# Merge lumpy SV calls

# Input and output files
LUMPY_VCF_PATH=/path/to/sample/vcfs # Use the VCFs generated from script 1_run_lumpy.sh
OUTPUT_PATH=/path/to/output

# Reference files
REF=/path/to/hg19/reference.fasta

# Paths to tools
LUMPY=/path/to/lumpy-sv/bin

# Merge calls from all samples
$LUMPY/smoove merge --name merged -f $REF --outdir $OUTPUT_PATH $LUMPY_VCF_PATH/*.genotyped.vcf.gz
