#!/bin/bash

# Combine all sample VCFs into a single table
# This step will need to be run for all sample VCFs from script 8

# Paths to tools
BCFPATH=/path/to/bcftools_v1.12

# Input and output files
INPUT_VCF=/path/to/input.vcf.gz # Use the VCFs generated from script 8_annotate_cadd.sh here
SAMPLE="sample_id"
OUTPUT_FILE=/path/to/output/table.txt

echo $SAMPLE > samplename.txt
$BCFPATH/bcftools reheader -s samplename.txt $INPUT_VCF | $BCFPATH/bcftools query -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Func.wgEncodeGencodeBasicV38\t%Gene.wgEncodeGencodeBasicV38\t%GeneDetail.wgEncodeGencodeBasicV38\t%ExonicFunc.wgEncodeGencodeBasicV38\t%AAChange.wgEncodeGencodeBasicV38\t%ge_AF\t%ge_AF_nfe\t%ge_AF_asj\t%ge_AF_eas\t%ge_AF_amr\t%ge_AF_afr\t%gg_AF\t%gg_AF_nfe\t%gg_AF_asj\t%gg_AF_eas\t%gg_AF_amr\t%gg_AF_afr\t%CADD_PHRED\t%CADD_RawScore[\t%GT\t%DP\t%AD\t%GQ\t%PL]\n' $INPUT_VCF >> $OUTPUT_FILE
rm samplename.txt
