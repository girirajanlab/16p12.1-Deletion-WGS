#!/bin/bash

# Use lumpy to regenotype samples on the merged call set

# Input and output files
LUMPY_MERGED_VCF=/path/to/merged/calls.vcf # Use the VCFs generated from script 2_merge_calls.sh
INPUT_BAM=/path/to/input.bam # Use the sorted BAM files generated from bwa-mem in script 1_GATK_pipeline.sh
OUTPUT_PATH=/path/to/output
SAMPLE='Sample_ID'

# Reference files
REF=/path/to/hg19/reference.fasta

# Paths to tools
LUMPY=/path/to/lumpy-sv/bin

# Regenotype samples using the merged call set
$LUMPY/smoove genotype -d -x -p 20 --name $SAMPLE --outdir $OUTPUT_PATH --fasta $REF --vcf $LUMPY_MERGED_VCF $INPUT_BAM
