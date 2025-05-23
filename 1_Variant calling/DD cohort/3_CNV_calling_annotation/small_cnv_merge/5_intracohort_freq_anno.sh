#!/bin/bash

# Input and output files
CNV_PATH=/path/to/other/cnvs/ # Use the path to the cnvs from script 4_merge_split.py
OUTPUT_PATH=/path/to/output/files/

# Paths to tools
BEDPATH=/path/to/bedtools_v2.27.1

# Annotate intracohort counts
$BEDPATH/bedtools intersect -a $INPUT_PATH/dels.bed -b $INPUT_PATH/dels.bed -f 0.5 -r -c > $OUTPUT_PATH/intracohort_dels.bed
$BEDPATH/bedtools intersect -a $INPUT_PATH/dups.bed -b $INPUT_PATH/dups.bed -f 0.5 -r -c > $OUTPUT_PATH/intracohort_dups.bed
