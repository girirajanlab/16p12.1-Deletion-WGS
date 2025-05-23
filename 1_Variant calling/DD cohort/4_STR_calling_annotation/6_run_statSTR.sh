#!/bin/bash

# Identify STR expansions across the cohort

# Input and output files
VCFS="/path/to/family/vcfs.vcf,/as/a/comma/separated/list.vcf"
OUTPUT_DIR=/path/to/output/directory

# Tools
MERGESTR=/path/to/mergeSTR
STATSTR=/path/to/statSTR

# Merge all family VCFs into a single file
$MERGESTR/mergeSTR --vcfs $VCFS --out $OUTPUT_DIR/cohort

# Run statSTR to get statistics for the cohort
$STATSTR/statSTR --vcf $OUTPUT_DIR/cohort.vcf.gz \
	--thresh \
	--afreq \
	--acount \
	--mean \
	--mode \
	--var \
	--numcalled \
	--out $OUTPUT_DIR/cohort.statstr


