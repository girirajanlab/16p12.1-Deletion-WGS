#!/bin/bash

# Calculate ancestry using somalier

# Input and output files
GVCF_PATH=/path/to/SFARI/SPARK/iWES/input/gvcfs # Path to individual sample gVCFs
# The input gVCFs needed can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
SAMPLE_LIST=/path/to/list/of/samples.txt
OUTPUT_DIR=/path/to/output/directory

# Reference files
REF=/path/to/hg38/reference.fasta
# Reference files can be downloaded from the somalier github page
SITES=/path/to/hg38/sites.vcf.gz # https://github.com/brentp/somalier/releases/tag/v0.2.19
ANC_LABELS=/apth/to/ancestry/labels.tsv # https://github.com/brentp/somalier/wiki/ancestry
SOMALIER_GENO=/path/to/somalier/genotype/files/*.somalier # https://github.com/brentp/somalier/wiki/ancestry

# Tools
BCFPATH=/path/to/bcftools_v1.12
TABIXPATH=/path/to/tabix

# Extract sites needed for somalier
while read SAMPLE
do
	LAST_DIGIT="${SAMPLE: -1}"
	INPUT_VCF=$GVCF_PATH/$LAST_DIGIT/$SAMPLE.gvcf.gz
	$BCFPATH/bcftools view -R $SITES $INPUT_VCF | bgzip > $OUTPUT_DIR/$LAST_DIGIT/$SAMPLE.sites.vcf.gz
	$TABIXPATH/tabix -p vcf $OUTPUT_DIR/$LAST_DIGIT/$SAMPLE.sites.vcf.gz

	somalier extract -d extracted/$LAST_DIGIT/ -s $SITES -f $REF --sample-prefix=$SAMPLE $OUTPUT_DIR/$LAST_DIGIT/$SAMPLE.sites.vcf.gz
done < $SAMPLE_LIST

# Calculate ancestry, splitting samples by last digit
for i in {0..9}
do
	somalier ancestry -o somalier_ancestry_$i --labels $ANC_LABELS $SOMALIER_GENO ++ extracted/$i/*.somalier
done
