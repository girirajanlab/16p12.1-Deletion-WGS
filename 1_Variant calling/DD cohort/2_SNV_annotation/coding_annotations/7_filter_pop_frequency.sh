#!/bin/bash

# Input and output files
INPUT_VCF=/path/to/input.vcf.gz # Use the VCF generated from script 4_annotate_gnomad.sh here
OUTPUT_VCF=/path/to/output.vcf.gz
ANC_FIle=/path/to/ancestry/annotations/file.csv # Use the output from script 6_combine_somalier_selfreport.py here

# Tools for processing VCFs
BCFPATH=/path/to/bcftools_v1.12
TABIXPATH=/path/to/tabix

# This step is run separately on each sample
SAMPLE=sample_id

mkdir -p tmp/$SAMPLE

# All sample gnomad filter
$BCFPATH/bcftools view -i 'ge_AF<=0.001 | ge_AF="."' $INPUT_VCF | bcftools view -i 'gg_AF<=0.001 | gg_AF="."' | bgzip > tmp/$SAMPLE/AF.vcf.gz

# Population specific gnomad filter
EUR=`grep $SAMPLE $ANC_FILE | cut -f 5 -d,`
ASJ=`grep $SAMPLE $ANC_FILE | cut -f 6 -d,`
EAS=`grep $SAMPLE $ANC_FILE | cut -f 7 -d,`
AMR=`grep $SAMPLE $ANC_FILE | cut -f 8 -d,`
AFR=`grep $SAMPLE $ANC_FILE | cut -f 9 -d,`

# For each population, if the sample is annotated for that population,  apply the filter for that population
INVCF=tmp/$SAMPLE/AF.vcf.gz
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

cp $INVCF $OUTPUT_VCF
$TABIXPATH/tabix -p vcf $OUTPUT_VCF
