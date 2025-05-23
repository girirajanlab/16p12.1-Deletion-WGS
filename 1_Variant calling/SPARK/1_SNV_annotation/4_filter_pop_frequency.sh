#!/bin/bash

# Input and output files
SCRIPT1_OUTPUTDIR=/path/to/output/directory/from/script/1_annotate_variants/ # Use the output directory defined in script 1_annotate_variants.sh
OUTPUT_DIR=/path/to/output/directory
ANC_FILE=/path/to/ancestry/annotations/file.csv # Use the output from script 3_combine_somalier_reported.py here
SAMPLE_LIST=/path/to/list/of/samples.txt
OUTPUT_TABLE=/path/to/output/variant/table.tsv

# Tools
BCFPATH=/path/to/bcftools_v1.12
TABIXPATH=/path/to/tabix

while read SAMPLE
do
	LAST_DIGIT="${SAMPLE: -1}"
	INVCF=$SCRIPT1_OUTPUTDIR/annotated/$LAST_DIGIT/$SAMPLE.vcf.gz

	# Apply population-specific gnomAD filters
	EUR=`grep $SAMPLE $ANC_FILE | cut -f 3 -d,`
	ASJ=`grep $SAMPLE $ANC_FILE | cut -f 4 -d,`
	EAS=`grep $SAMPLE $ANC_FILE | cut -f 5 -d,`
	AMR=`grep $SAMPLE $ANC_FILE | cut -f 6 -d,`
	AFR=`grep $SAMPLE $ANC_FILE | cut -f 7 -d,`

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

	# Save VCF as table
	$BCFPATH/bcftools query -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Func.wgEncodeGencodeBasicV38\t%Gene.wgEncodeGencodeBasicV38\t%GeneDetail.wgEncodeGencodeBasicV38\t%ExonicFunc.wgEncodeGencodeBasicV38\t%AAChange.wgEncodeGencodeBasicV38\t%ge_AF\t%ge_AF_nfe\t%ge_AF_asj\t%ge_AF_eas\t%ge_AF_amr\t%ge_AF_afr\t%gg_AF\t%gg_AF_nfe\t%gg_AF_asj\t%gg_AF_eas\t%gg_AF_amr\t%gg_AF_afr\t%CADD_PHRED\t%CADD_RawScore[\t%GT\t%DP\t%AD\t%GQ\t%PL]\n' $INVCF > $OUTPUT_DIR/tables/$LAST_DIGIT/$SAMPLE.tsv

	# Concat all tables into a single file
	cat $OUTPUT_DIR/tables/$LAST_DIGIT/$SAMPLE.tsv >> $OUTPUT_TABLE

done < $SAMPLE_LIST
