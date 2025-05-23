#!/bin/bash

# Input files
SAMPLE_LIST=/path/to/list/of/samples.csv
VCF_PATH="/path/to/sample/vcfs/" # Use the output directory of script 4_annotate_gnomad.sh here

# Reference files
REF=/path/to/hg38_reference.fasta
# Reference files can be downloaded from the somalier github page
SITES=/path/to/hg38/sites.vcf.gz # https://github.com/brentp/somalier/releases/tag/v0.2.19
ANC_LABELS=/path/to/ancestry/labels.tsv # https://github.com/brentp/somalier/wiki/ancestry
SOMALIER_GENO=/path/to/somalier/genotype/files/*.somalier # https://github.com/brentp/somalier/wiki/ancestry

# Run somalier to determine ancestry on all samples
# Extract needed sites
while read sample
do
    somalier extract -d extracted/ -s $SITES -f $REF --sample-prefix=$sample $VCF_PATH/$sample.vcf.gz
done < $SAMPLE_LIST

# Calculate ancestry
somalier ancestry -o somalier_ancestry --labels $ANC_LABELS $SOMALIER_GENO ++ extracted/*.somalier

