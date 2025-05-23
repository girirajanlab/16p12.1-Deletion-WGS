#!/bin/bash

# Run LDPred2 to calculate PRS

# Input and output files
INPUT_PLINK=/path/to/input/plink/files # Use the output of script 3_vcf2plink.sh
EUR_SAMPS=/path/to/European/samples.fam # Use the output of script 5_combine_peddy_reported.py
OUTPUT_DIR=/path/to/output/directory

# Reference files
ASD_SUMMSTAT=/path/to/Autism/GWAS/summary/statistics/file.txt.gz
INTELL_SUMMSTAT=/path/to/Intelligence/GWAS/summary/statistics/file.txt.gz
EDU_SUMMSTAT=/path/to/Educational_Attainment/GWAS/summary/statistics/file.txt.gz
SCHIZ_SUMMSTAT=/path/to/Schizophrenia/GWAS/summary/statistics/file.txt.gz

# Tools
PLINK=/path/to/PLINK_v1.90

# Subset out European samples from PLINK
$PLINK --bfile $INPUT_PLINK --keep $EUR_SAMPS --make-bed --out $OUTPUT_DIR/EUROPEAN

# Run LDPred2
Rscript helper_scripts/ldpred2.R $OUTPUT_DIR/EUROPEAN $ASD_SUMMSTAT autism $OUTPUT_DIR
Rscript helper_scripts/ldpred2.R $OUTPUT_DIR/EUROPEAN $INTELL_SUMMSTAT intelligence $OUTPUT_DIR
Rscript helper_scripts/ldpred2.R $OUTPUT_DIR/EUROPEAN $EDU_SUMMSTAT educational_attainment $OUTPUT_DIR
Rscript helper_scripts/ldpred2.R $OUTPUT_DIR/EUROPEAN $SCHIZ_SUMMSTAT schizophrenia $OUTPUT_DIR

