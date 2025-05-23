#!/bin/bash

# Input and output files
SCRIPT1_OUTPUTDIR=/path/to/output/directory/from/script/1_annotate_variants/ # Use the output directory defined in script 1_annotate_variants.sh
OUTPUT_TABLE=/path/to/output/table.tsv
ANC_FILE=/path/to/ancestry/annotations/file.csv # Use the output from script 3_combine_somalier_reported.py here
SAMPLE=sample_id

# Tools
BCFPATH=/path/to/bcftools_v1.12
TABIXPATH=/path/to/tabix

mkdir -p tmp/$SAMPLE

# Apply population-specific gnomAD filters
EUR=`grep $SAMPLE $ANC_FILE | cut -f 2 -d,`
ASJ=`grep $SAMPLE $ANC_FILE | cut -f 3 -d,`
EAS=`grep $SAMPLE $ANC_FILE | cut -f 4 -d,`
AMR=`grep $SAMPLE $ANC_FILE | cut -f 5 -d,`
AFR=`grep $SAMPLE $ANC_FILE | cut -f 6 -d,`

# For each population, if the sample is annotated for that population,  apply the filter for that population
INVCF=$SCRIPT1_OUTPUTDIR/sample_split/$SAMPLE.vcf.gz
if [ "$EUR" -eq "1" ]
then
	$BCFPATH/bcftools view -i 'ge_AF_nfe<=0.001 | ge_AF_nfe="."' $INVCF | $BCFPATH/bcftools view -i 'gg_AF_nfe<=0.001 | gg_AF_nfe="."' | bgzip > tmp/$SAMPLE/AF_EUR.vcf.gz
	INVCF=tmp/$SAMPLE/AF_EUR.vcf.gz
fi

if [ "$ASJ" -eq "1" ]
then
	$BCFPATH/bcftools view -i 'ge_AF_asj<=0.001 | ge_AF_asj="."' $INVCF | $BCFPATH/bcftools view -i 'gg_AF_asj<=0.001 | gg_AF_asj="."' | bgzip > tmp/$SAMPLE/AF_ASJ.vcf.gz
	INVCF=tmp/$SAMPLE/AF_ASJ.vcf.gz
fi

if [ "$EAS" -eq "1" ]
then
	$BCFPATH/bcftools view -i 'ge_AF_eas<=0.001 | ge_AF_eas="."' $INVCF | $BCFPATH/bcftools view -i 'gg_AF_eas<=0.001 | gg_AF_eas="."' | bgzip > tmp/$SAMPLE/AF_EAS.vcf.gz
	INVCF=tmp/$SAMPLE/AF_EAS.vcf.gz
fi

if [ "$AMR" -eq "1" ]
then
	$BCFPATH/bcftools view -i 'ge_AF_amr<=0.001 | ge_AF_amr="."' $INVCF | $BCFPATH/bcftools view -i 'gg_AF_amr<=0.001 | gg_AF_amr="."' | bgzip > tmp/$SAMPLE/AF_AMR.vcf.gz
	INVCF=tmp/$SAMPLE/AF_AMR.vcf.gz
fi

if [ "$AFR" -eq "1" ]
then
	$BCFPATH/bcftools view -i 'ge_AF_afr<=0.001 | ge_AF_afr="."' $INVCF | $BCFPATH/bcftools view -i 'gg_AF_afr<=0.001 | gg_AF_afr="."' | bgzip > tmp/$SAMPLE/AF_AFR.vcf.gz
	INVCF=tmp/$SAMPLE/AF_AFR.vcf.gz
fi

cp $INVCF $OUTPUT_DIR/filtered/$SAMPLE.vcf.gz
$TABIXPATH/tabix -p vcf $OUTPUT_DIR/filtered/$SAMPLE.vcf.gz

rm -r tmp/$SAMPLE

# Save VCF as table
$BCFPATH/bcftools query -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Func.wgEncodeGencodeBasicV19\t%Gene.wgEncodeGencodeBasicV19\t%GeneDetail.wgEncodeGencodeBasicV19\t%ExonicFunc.wgEncodeGencodeBasicV19\t%AAChange.wgEncodeGencodeBasicV19\t%ge_AF\t%ge_AF_nfe\t%ge_AF_asj\t%ge_AF_eas\t%ge_AF_amr\t%ge_AF_afr\t%gg_AF\t%gg_AF_nfe\t%gg_AF_asj\t%gg_AF_eas\t%gg_AF_amr\t%gg_AF_afr\t%CADD_PHRED\t%CADD_RawScore[\t%GT\t%DP\t%AD\t%GQ\t%PL]\n' $OUTPUT_DIR/filtered/$SAMPLE.vcf.gz >> $OUTPUT_TABLE
