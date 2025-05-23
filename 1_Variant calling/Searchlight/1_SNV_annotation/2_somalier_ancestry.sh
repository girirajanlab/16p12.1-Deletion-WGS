#!/bin/bash

# Calculate ancestry using somalier

# Input and output files
INPUT_VCF=/path/to/input/vcf.gz # Use the $OUTPUT_DIR/quality_filter.vcf.gz generated in script 1_annotate_variants.sh
OUTPUT_DIR=/path/to/output/directory
SAMPLE_LIST=/path/to/list/of/samples.txt

# Reference files
REF=/path/to/hg19/reference.fasta
# Reference files can be downloaded from the somalier github page
SITES=/path/to/hg19/sites.vcf.gz # https://github.com/brentp/somalier/releases/tag/v0.2.19
ANC_LABELS=/path/to/ancestry/labels.tsv # https://github.com/brentp/somalier/wiki/ancestry
SOMALIER_GENO=/path/to/somalier/genotype/files/*.somalier # https://github.com/brentp/somalier/wiki/ancestry

# Tools
BCFPATH=/path/to/bcftools_v1.12

# Extract needed sites for somalier
while read SAMPLE
do
	somalier extract -d extracted/ -s $SITES -f $REF --sample-prefix=$SAMPLE $INPUT_VCF
done < $SAMPLE_LIST

# Calculate ancestry
somalier ancestry -o somalier_ancestry --labels $ANC_LABELS $SOMALIER_GENO ++ extracted/*.somalier
