#!/bin/bash

# Split pVCF into individual VCF and split and left-normalize variants

# Input files
GVCF_INPUT=/path/to/pvcf.vcf.gz # Use the output of script 2_genotypeGVCF.sh here
SAMPLE_LIST=/path/to/list/of/samples.txt # a file listing all samples processed, one sample per line

# Tools for processing VCFs
VCFSUB=/path/to/vcf-subset
BCFPATH=/path/to/bcftools_v1.3.1

# Reference files
REF=/path/to/hg19/reference.fasta

# Split pVCF into into individual VCFs
while read sample; do $VCFSUB/vcf-subset --exclude-ref -c $sample $GVCF_INPUT > ${sample}.vcf; done < $SAMPLE_LIST

# Split multiple variants and left normalize
# Do this step separately on each individual VCF
while read sample
do
	$BCFPATH/bcftools norm -m-both -o ${sample}_split.vcf ${sample}.vcf
	$BCFPATH/bcftools norm -f $REF -o ${sample}_norm.vcf ${sample}_split.vcf
done < $SAMPLE_LIST
