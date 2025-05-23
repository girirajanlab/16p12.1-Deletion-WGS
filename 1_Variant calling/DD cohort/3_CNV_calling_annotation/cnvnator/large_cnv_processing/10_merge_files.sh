#!/bin/bash

# Concat files together

# Input and output files
INPUT_PATH=/path/to/cnvs/files/ # Use the path to the files from script 9_gnomadSV_filter.py
OUTPUT_PATH=/path/to/output/files/

cat $INPUT_PATH/cnv_gnomadSV_filter_dels.bed > $OUTPUT_PATH/filtered_cnvs.bed
tail -n+2 $INPUT_PATH/cnv_gnomadSV_filter_dups.bed >> $OUTPUT_PATH/filtered_cnvs.bed

cat $INPUT_PATH/pathogenic_cnv_gnomadSV_filter_dels.bed > $OUTPUT_PATH/filtered_pathogenic_cnvs.bed
tail -n+2 $INPUT_PATH/pathogenic_cnv_gnomadSV_filter_dups.bed >> $OUTPUT_PATH/filtered_pathogenic_cnvs.bed
