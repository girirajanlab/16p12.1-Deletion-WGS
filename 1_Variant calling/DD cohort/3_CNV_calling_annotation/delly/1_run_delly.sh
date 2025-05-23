#!/bin/bash

# Call CNVs using delly on a single sample

# Input and output files
INPUT_BAM=/path/to/input.bam # Use the sorted BAM files generated from bwa-mem in script 1_GATK_pipeline.sh
OUTPUT_FILE=/path/to/output.bcf

# Reference files
REF=/path/to/hg19/reference.fasta
EXCLUDE_SITES=/path/to/sites/to/exclude.bed

# Paths to tools
DELLY=/path/to/delly_v0.8.2

# Run delly to call CNVs
$DELLY call \
	-o $OUTPUT_FILE \
	-g $REF \
	-x $EXCLUDE_SITES \
	$INPUT_BAM

