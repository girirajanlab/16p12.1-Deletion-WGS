#!/bin/bash

# Use bedtools to find 50% reciprocal overlap CNVs from all samples

# Input and output files
INPUT=/path/to/input/file.bed # Use the output of script 12_finalize_calls.py
OUTPUT=/path/to/output/file.bed

# Paths to tools
BEDPATH=/path/to/bedtools_v2.27.1

# Run bedtools
$BEDPATH/bedtools intersect -a $INPUT -b $INPUT -f 0.5 -r -loj > $OUTPUT
