#!/bin/bash

# Use delly to regenotype samples on the merged call set

# Input and output files
DELLY_MERGED_BCF=/path/to/merged/calls.bcf # Use the BCFs generated from script 2_merge_calls.sh
INPUT_BAM=/path/to/input.bam # Use the sorted BAM files generated from bwa-mem in script 1_GATK_pipeline.sh
OUTPUT_FILE=/path/to/output.bcf

# Reference files
REF=/path/to/hg19/reference.fasta

# Paths to tools
DELLY=/path/to/delly_v0.8.2

# Regenotype samples using the merged call set
$DELLY call -g $REF -v $DELLY_MERGED_BCF -o $OUTPUT_FILE $INPUT_BAM
