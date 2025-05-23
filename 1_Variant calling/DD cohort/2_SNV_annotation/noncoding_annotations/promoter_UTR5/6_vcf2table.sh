#!/bin/bash

# Combine all sample VCFs into a single table
# This step will need to be run for all sample VCFs from script 5

# Paths to tools
BCFPATH=/path/to/bcftools_v1.12

# Input and output files
INPUT_VCF=/path/to/input.vcf.gz # Use the VCFs generated from script 5_filter_gnomad.sh here
OUTPUT_FILE=/path/to/output/table.txt

$BCFPATH/bcftools query -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Func.wgEncodeGencodeBasicV19\t%Gene.wgEncodeGencodeBasicV19\t%GeneDetail.wgEncodeGencodeBasicV19\t%gg_AF\t%gg_AF_nfe\t%gg_AF_asj\t%gg_AF_eas\t%gg_AF_amr\t%gg_AF_afr[\t%GT\t%DP\t%AD\t%GQ\t%PL\t%RGQ\t%SB]\n' $INPUT_VCF >> $OUTPUT_FILE
