#!/bin/bash

# Find overlaps with known STR regions

# Input and output files
INPUT_BED=/path/to/input.bed # Use the file generated from script 9_filter_lowconf_regions.sh
OUTPUT_BED=/path/to/output.bed
TMP_DIR=/path/to/tmp/directory

# Reference files
# This is a reference file of STRs identified by GangSTR from the Gymrek lab
# It can be downloaded here: https://s3.amazonaws.com/gangstr/hg19/genomewide/hg19_ver13_1.bed.gz
# And is described here: https://github.com/gymreklab/GangSTR
STR_REF=/path/to/hg19/strs/hg19_ver13_1.bed

# Paths to tools
BEDPATH=/path/to/bedtools_v2.27.1

# Use bedtools to find either
# (a) CNVs that do not overlap the STR regions
# <-------->                <---------------> CNV
#             <---->                          STR
$BEDPATH/bedtools intersect -wa -a $INPUT_BED -b $STR_REF -v -header > $TMP_DIR/str_no_overlap.bed
# (b) CNVS that overlap the STR regions 100%
# <-----------------------> CNV
#    <--->                  STR
$BEDPATH/bedtools intersect -wa -a $INPUT_BED -b $STR_REF -F 1 -u -header > $TMP_DIR/str_all_overlap.bed

# This second line would also get:
# <---------------------->    CNV
#    <---->            <----> STR
# So run a python script to filter the 100% overlap calls for those with a breakpoint in an STR
python 10.1_str_filter.py $TMP_DIR/str_all_overlap.bed $TMP_DIR

# Merge two files
cat $TMP_DIR/str_no_overlap.bed > $TMP_DIR/str_filter.bed
cat $TMP_DIR/str_all_filter.bed >> $TMP_DIR/str_filter.bed

# Add a header
head -1 $TMP_DIR/str_no_overlap.bed > $OUTPUT_BED
tail -n+2 $TMP_DIR/str_filter.bed | sort | uniq >> $OUTPUT_BED

rm -r $TMP_DIR
