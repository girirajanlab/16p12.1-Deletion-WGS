#!/bin/bash

# Input and output files
INPUT_BAM=/path/to/input.bam # Use the sorted BAM files generated from bwa-mem in script 1_GATK_pipeline.sh
OUTPUT_PATH=/path/to/output/dir
OUTPUT_NAME='output_filename'

# Reference files
REF=/path/to/hg19/reference.fasta

# Paths to tools
LUMPY=/path/to/lumpy-sv/bin

# A list of chromosomes to exclude
exclude_chroms=chr6_ssto_hap7,chr6_mcf_hap5,chr6_cox_hap2,chr6_mann_hap4,chr6_apd_hap1,chr6_qbl_hap6,chr6_dbb_hap3,chr17_ctg5_hap1,chr4_ctg9_hap1,chr1_gl000192_random,chrUn_gl000225,chr4_gl000194_random,chr4_gl000193_random,chr9_gl000200_random,chrUn_gl000222,chrUn_gl000212,chr7_gl000195_random,chrUn_gl000223,chrUn_gl000224,chrUn_gl000219,chr17_gl000205_random,chrUn_gl000215,chrUn_gl000216,chrUn_gl000217,chr9_gl000199_random,chrUn_gl000211,chrUn_gl000213,chrUn_gl000220,chrUn_gl000218,chr19_gl000209_random,chrUn_gl000221,chrUn_gl000214,chrUn_gl000228,chrUn_gl000227,chr1_gl000191_random,chr19_gl000208_random,chr9_gl000198_random,chr17_gl000204_random,chrUn_gl000233,chrUn_gl000237,chrUn_gl000230,chrUn_gl000242,chrUn_gl000243,chrUn_gl000241,chrUn_gl000236,chrUn_gl000240,chr17_gl000206_random,chrUn_gl000232,chrUn_gl000234,chr11_gl000202_random,chrUn_gl000238,chrUn_gl000244,chrUn_gl000248,chr8_gl000196_random,chrUn_gl000249,chrUn_gl000246,chr17_gl000203_random,chr8_gl000197_random,chrUn_gl000245,chrUn_gl000247,chr9_gl000201_random,chrUn_gl000235,chrUn_gl000239,chr21_gl000210_random,chrUn_gl000231,chrUn_gl000229,chrUn_gl000226,chr18_gl000207_random

# Run lumpy/smoove to call CNVs
$LUMPY/smoove call --name $OUTPUT_NAME \
 --fasta $REF \
 --outdir $OUTPUT_PATH \
 --excludechroms $exclude_chroms \
 --genotype \
 --duphold \
 $INPUT_BAM
