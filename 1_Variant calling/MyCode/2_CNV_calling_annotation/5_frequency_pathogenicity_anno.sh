#!/bin/bash

# Input and output files
INPUT_DIR="/path/to/input/bed/files" # Use the output directory from script 4_size_filter.py
OUTPUT_DIR="/path/to/output/directory"

# Reference files
CONTROL_CNVS=/path/to/frequency/data/from/control/microarray/cohort/calls.bed
# These data are from Coe et al. Nat Genet 2014 (https://pubmed.ncbi.nlm.nih.gov/25217958/) and contain microarray CNV calls from a control cohort
PATHOGENIC_CNV=/path/to/known/pathogenic/cnv/regions/pathogenic_cnv_hg38.bed # A bed file of known pathogenic CNV regions

# Tools
BEDPATH=/path/to/bedtools_v2.27.1

# Annotate for frequency in control cohort
$BEDPATH/bedtools intersect -a $INPUT_DIR/dels.bed -b $CONTROL_CNVS -f 0.5 -r -c > $OUTPUT_DIR/control_anno_dels.bed
$BEDPATH/bedtools intersect -a $INPUT_DIR/dups.bed -b $CONTROL_CNVS -f 0.5 -r -c > $OUTPUT_DIR/control_anno_dups.bed

# Annotate for intracohort frequency
$BEDPATH/bedtools intersect -a $OUTPUT_DIR/control_anno_dels.bed -b $OUTPUT_DIR/control_anno_dels.bed -f 0.5 -r -c > $OUTPUT_DIR/intracohort_anno_dels.bed
$BEDPATH/bedtools intersect -a $OUTPUT_DIR/control_anno_dups.bed -b $OUTPUT_DIR/control_anno_dups.bed -f 0.5 -r -c > $OUTPUT_DIR/intracohort_anno_dups.bed

# Annotate for overlap with known pathogenic CNV regions
$BEDPATH/bedtools intersect -wa -wb -a $OUTPUT_DIR/intracohort_anno_dels.bed -b $PATHOGENIC_CNV -f 0.5 -r > $OUTPUT_DIR/pathogenic_anno_dels.bed
$BEDPATH/bedtools intersect -wa -wb -a $OUTPUT_DIR/intracohort_anno_dels.bed -b $PATHOGENIC_CNV -f 0.5 -r > $OUTPUT_DIR/pathogenic_anno_dups.bed

