#!/bin/bash

# Perform QC filters on lumpy calls and split into individual files for each sample

# Input and output files
INPUT_VCF=/path/to/input.vcf.gz # Use the VCFs generated from script 4_paste.sh
FILTER_VCF=/path/to/filtered.vcf
OUTPUT_PATH=/path/to/output/files
SAMPLE_LIST=/path/to/list/of/samples.txt # a file listing all samples processed, one sample per line

# Paths to tools
BCFPATH=/path/to/bcftools_v1.12

# Only keep heterozygous dels with DHFFC < 0.75, heterozygous dups with DHFFC > 1.25, homozygous dels with DHFFC < 0.5, and homozygous dups with DHFFC > 1.5
$BCFPATH/bcftools filter --set-GT . -i '(DHFFC < 0.75 & SVTYPE="DEL" & GT="HET" & GT="ALT") |
	(DHFFC > 1.25 & SVTYPE="DUP" & GT="HET" & GT="ALT") |
	(DHFFC < 0.5 & SVTYPE="DEL" & GT="HOM" & GT="ALT") |
	(DHFFC > 1.5 & SVTYPE="DUP" & GT="HOM" & GT="ALT")' $INPUT_VCF > $FILTER_VCF

# Split VCF into individual VCFs for each sample
while read SAMPLE
do
	$BCFPATH/bcftools view -s $SAMPLE $FILTER_VCF | $BCFPATH/bcftools view -i 'FORMAT/GT="alt"' > $OUTPUT_PATH/$SAMPLE.vcf
done
