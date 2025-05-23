#!/bin/bash

# Input and output files
PATHOGENIC_PATH=/path/to/pathogenic/cnvs/ # Use the path to the pathogenic cnvs identified in script 4_filter.sh
CNV_PATH=/path/to/other/cnvs/ # Use the path to the cnvs filtered in script 5_str_filter.py
OUTPUT_PATH=/path/to/output/files/

# Reference files
CONTROL_PATH=/path/to/frequency/data/from/control/microarray/cohort/
# These data are from Coe et al. Nat Genet 2014 (https://pubmed.ncbi.nlm.nih.gov/25217958/) and contain microarray CNV calls from a control cohort

# Paths to tools
BEDPATH=/path/to/bedtools_v2.27.1

# Annotate CNVs with counts in a control microarray cohort
$BEDPATH/bedtools intersect -a $CNV_PATH/dels.bed -b $CONTROL_DATA/dels.bed -f 0.5 -r -c > $OUTPUT_PATH/cnv_control_dels.bed
$BEDPATH/bedtools intersect -a $CNV_PATH/dups.bed -b $CONTROL_DATA/dups.bed -f 0.5 -r -c > $OUTPUT_PATH/cnv_control_dups.bed
$BEDPATH/bedtools intersect -a $PATHOGENIC_PATH/pathogenic_dels.bed -b $CONTROL_DATA/dels.bed -f 0.5 -r -c > $OUTPUT_PATH/pathogenic_cnv_control_dels.bed
$BEDPATH/bedtools intersect -a $PATHOGENIC_PATH/pathogenic_dups.bed -b $CONTROL_DATA/dups.bed -f 0.5 -r -c > $OUTPUT_PATH/pathogenic_cnv_control_dups.bed

# Annotate intracohort counts
$BEDPATH/bedtools intersect -a $OUTPUT_PATH/cnv_control_dels.bed -b $OUTPUT_PATH/cnv_control_dels.bed -f 0.5 -r -c > $OUTPUT_PATH/cnv_control_intracohort_dels.bed
$BEDPATH/bedtools intersect -a $OUTPUT_PATH/cnv_control_dups.bed -b $OUTPUT_PATH/cnv_control_dups.bed -f 0.5 -r -c > $OUTPUT_PATH/cnv_control_intracohort_dups.bed
$BEDPATH/bedtools intersect -a $OUTPUT_PATH/pathogenic_cnv_control_dels.bed -b $OUTPUT_PATH/pathogenic_cnv_control_dels.bed -f 0.5 -r -c > $OUTPUT_PATH/pathogenic_cnv_control_intracohort_dels.bed
$BEDPATH/bedtools intersect -a $OUTPUT_PATH/pathogenic_cnv_control_dups.bed -b $OUTPUT_PATH/pathogenic_cnv_control_dups.bed -f 0.5 -r -c > $OUTPUT_PATH/pathogenic_cnv_control_intracohort_dups.bed
