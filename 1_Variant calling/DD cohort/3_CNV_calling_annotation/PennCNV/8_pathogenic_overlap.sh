#!/bin/bash

# Input and output files
INPUT_FILE=/path/to/input.csv # Use the file generated from 7_gnomad_frequency.py
OUTPUTDIR=/path/to/output/directory

# Reference files
PATHOGENIC_CNV=/path/to/known/pathogenic/cnv/regions/pathogenic_cnv_hg19.bed # A bed file of known pathogenic CNV regions

# Paths to tools
BEDPATH=/path/to/bedtools_v2.27.1

# Convert CNV calls into BED format
cut -f 1-3,6,9 -d , $INPUT_FILE | sed 's/,/\t/g' | tail -n+2 > tmp.bed

# Find regions that have a 50% reciprocal overlap with known pathogenic CNVs
# These CNVs will not be subject to the control population frequency filters
# Also separate deletions and duplications
$BEDPATH/bedtools intersect -wa -wb -a tmp.bed -b $PATHOGENIC_CNV -f 0.5 -r | grep "<DEL>" > $OUTPUTDIR/pathogenic_dels.bed

rm tmp.bed
