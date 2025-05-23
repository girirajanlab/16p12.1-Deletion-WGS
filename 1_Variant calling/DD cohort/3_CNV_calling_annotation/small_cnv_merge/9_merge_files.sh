#!/bin/bash

# Concat files together

# Input and output files
INPUT_PATH=/path/to/cnvs/files/ # Use the path to the files from script 8_gnomadSV_filter.py
OUTPUT_PATH=/path/to/output/files/

cat $INPUT_PATH/gnomadSV_filter_dels.bed > $OUTPUT_PATH/filtered_cnvs.bed
tail -n+2 $INPUT_PATH/gnomadSV_filter_dups.bed >> $OUTPUT_PATH/filtered_cnvs.bed
