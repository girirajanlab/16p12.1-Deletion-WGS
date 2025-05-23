#!/bin/bash

# Merge CNV calls from multiple samples together
# And remove calls that have a 50% reciprocal overlap with low-confidence regions

# Input and output files
INPUT_BED=/path/to/list/of/CNV/files.txt # The output from script 3_size_filter.py should be used here
OUTPUT_BED=/path/to/output_bedfile.bed

# Reference files
FILTERED_BED=/path/to/hg19/filtered_regions.bed # File can be downloaded from http://bit.ly/1PDkVPQ and is described in https://www.sciencedirect.com/science/article/pii/S0002929716000690?via%3Dihub#app3
# This file is from Brandler et al. Science 2016 and contains low confidence regions, including:
# 	- Centromeres
#	- Segmental duplications (genomicSuperDups)
#	- Regions with low mappability and 100 bp reads (wgEncodeCrgMapabilityAlign100-mer), and 
#	- Regions subject to somatic V(D)J recombination (parts of antibodies and T cell receptor genes)

# Paths to tools
BEDTOOLSDIR=/path/to/bedtools_v2.27.1

# Use bedtools to remove calls that have a 50% reciprocal overlap with centromeres, SegDups, regions of low mappability, and V(D)J recombination regions from Brandler et al. Science 2016
$BEDTOOLSDIR/bedtools intersect -wa -a $INPUT_BED -b $FILTERED_BED -f 0.5 -r -v -header > $OUTPUT_BED
