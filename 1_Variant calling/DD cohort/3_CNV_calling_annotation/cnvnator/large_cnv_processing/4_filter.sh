#!/bin/bash

# Input and output files
INPUT_FILE=/path/to/input.bed # Use the file generated from 3_split_by_size.py
OUTPUTDIR=/path/to/output/directory

# Reference files
PATHOGENIC_CNV=/path/to/known/pathogenic/cnv/regions/pathogenic_cnv_hg19.bed # A bed file of known pathogenic CNV regions
FILTERED_BED=/path/to/hg19/filtered_regions.bed # File can be downloaded from http://bit.ly/1PDkVPQ and is described in https://www.sciencedirect.com/science/article/pii/S0002929716000690?via%3Dihub#app3
# This file is from Brandler et al. Science 2016 and contains low confidence regions, including:
# 	- Centromeres
#	- Segmental duplications (genomicSuperDups)
#	- Regions with low mappability and 100 bp reads (wgEncodeCrgMapabilityAlign100-mer), and 
#	- Regions subject to somatic V(D)J recombination (parts of antibodies and T cell receptor genes)

# Paths to tools
BEDPATH=/path/to/bedtools_v2.27.1

# Find regions that have a 50% reciprocal overlap with known pathogenic CNVs
# These CNVs will not be subject to the STR filter or the control population frequency filters
# Also separate deletions and duplications
$BEDPATH/bedtools intersect -wa -wb -a $INPUT_FILE -b $PATHOGENIC_CNV -f 0.5 -r | grep "<DEL>" > $OUTPUTDIR/pathogenic_dels.bed
$BEDPATH/bedtools intersect -wa -wb -a $INPUT_FILE -b $PATHOGENIC_CNV -f 0.5 -r | grep "<DUP>" > $OUTPUTDIR/pathogenic_dups.bed

# Find regions that have less than a 50% reciprocal overlap with centromeres, SegDups, regions of low mappability, and V(D)J recombination regions from Brandler et al. Science 2016
$BEDPATH/bedtools intersect -wa -a $INPUT_FILE -b $FILTERED_BED -f 0.5 -r -v -header > $OUTPUTDIR/low_conf_filter.bed

