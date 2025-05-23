#!/bin/bash

# Combine all sample VCFs into a single table
# This step will need to be run for all sample VCFs from script 9

# Paths to tools
BCFPATH=/path/to/bcftools_v1.12

# Input and output files
INPUT_VCF=/path/to/input.vcf.gz # Use the VCFs generated from script 9_annotate_ClinVar.sh here
OUTPUT_FILE=/path/to/output/table.txt

$BCFPATH/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Func.refGene\t%Gene.refGene\t%GeneDetail.refGene\t%ExonicFunc.refGene\t%AAChange.refGene\t%Func.wgEncodeGencodeBasicV19\t%Gene.wgEncodeGencodeBasicV19\t%GeneDetail.wgEncodeGencodeBasicV19\t%ExonicFunc.wgEncodeGencodeBasicV19\t%AAChange.wgEncodeGencodeBasicV19\t%ge_AF\t%ge_AF_nfe\t%ge_AF_asj\t%ge_AF_eas\t%ge_AF_amr\t%ge_AF_afr\t%gg_AF\t%gg_AF_nfe\t%gg_AF_asj\t%gg_AF_eas\t%gg_AF_amr\t%gg_AF_afr\t%CADD_PHRED\t%CADD_RawScore\t%ClinVar_CLNDN\t%ClinVar_CLNDISDB\t%ClinVar_CLNREVSTAT\t%ClinVar_CLNSIG\t%ClinVar_ALLELEID[\t%SAMPLE\t%GT\t%DP\t%AD\t%GQ\t%PL\t%RGQ\t%SB]\n' $INPUT_VCF >> $OUTPUT_FILE

