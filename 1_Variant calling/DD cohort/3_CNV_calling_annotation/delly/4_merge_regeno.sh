#!/bin/bash

# Merge BCFs from delly regenotyping

# Input and output files
DELLY_BCF_PATH=/path/to/regenotyped/bcfs/ # Use the BCFs generated from script 3_regenotype.sh
OUTPUT_FILE=/path/to/output.bcf

# Paths to tools
BCFPATH=/path/to/bcftools_v1.12

# Merge calls from all samples
$BCFPATH/bcftools merge -m id -O b -o $OUTPUT_FILE $DELLY_BCF_PATH/*.bcf
$BCFPATH/bcftools index $OUTPUT_FILE
