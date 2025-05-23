#!/bin/bash

# Annotate variants with functional consequences and delteriousness metrics

# Input and output files
GVCF_PATH=/path/to/SFARI/SPARK/iWES/input/gvcfs # Path to individual sample gVCFs
# The input gVCFs needed can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
SAMPLE_LIST=/path/to/sample/list.txt
OUTPUT_DIR=/path/to/output/directory

# Reference files
REF=/path/to/hg38/reference.fasta
ANNOVAR_ANNO=/path/to/annovar/annotations
GTOML=/path/to/gnomad_annotation_file.toml
CTOML=/path/to/cadd_annotation_file.toml

# Tools
BCFPATH=/path/to/bcftools_v1.12
ANNOVAR_DIR=/path/to/annovar
VCFANNO_PATH=/path/to/vcfanno
TABIXPATH=/path/to/tabix
JVARKIT=/path/to/jvarkit

# Split multiallelic records and left-normalize variants
# And filter for:
# 1. Non-ref alleles
# 2. QUAL>=50
# 3. Read depth>=8
# 4. Alternative allele depth>0
# 5. Allele balance >=0.25 and <=0.75 or >=0.9
# 6. QUAL/(Alternative allele depth)>=1.5
while read SAMPLE
do
	LAST_DIGIT="${SAMPLE: -1}"
	INPUT_VCF=$GVCF_PATH/$LAST_DIGIT/$SAMPLE.gvcf.gz
	$BCFPATH/bcftools norm -f $REF -m-both $INPUT_VCF | $BCFPATH/bcftools view -i 'GT="alt" & QUAL>=50 & FORMAT/DP>=8 & FORMAT/AD[:1]>0 & (FORMAT/AD[:1])/(FORMAT/DP)>=0.25 & QUAL/(FORMAT/AD[:1])>=1.5' | $BCFPATH/bcftools view -i '(FORMAT/AD[:1])/(FORMAT/DP)<=0.75 | (FORMAT/AD[:1])/(FORMAT/DP)>=0.9' -o $OUTPUT_DIR/quality_filter/$LAST_DIGIT/$SAMPLE.vcf.gz -Oz
	$TABIXPATH/tabix -p vcf $OUTPUT_DIR/quality_filter/$LAST_DIGIT/$SAMPLE.vcf.gz
done < $SAMPLE_LIST

# Merge all files into a single VCF whlie dropping genotype information
# Due to the number of files, first merge by the last two digits, then merge again
for i in `seq 0 99`
do
	printf -v LAST_TWO "%02d" $i
	LAST_DIGIT="${LAST_TWO: -1}"
	$BCFPATH/bcftools merge -m none $OUTPUT_DIR/quality_filter/$LAST_DIGIT/*$LAST_TWO.vcf.gz  | $BCFPATH/bcftools view -G | bgzip > $OUTPUT_DIR/merge1/$LAST_TWO.vcf.gz
	$TABIXPATH/tabix -p vcf $OUTPUT_DIR/merge1/$LAST_TWO.vcf.gz
done

$BCFPATH/bcftools merge -m none $OUTPUT_DIR/merge1/*.vcf.gz | $BCFPATH/bcftools view -G | bgzip > $OUTPUT_DIR/merge.vcf.gz
$TABIXPATH/tabix -p vcf $OUTPUT_DIR/merge.vcf.gz

# Annotate variant effects with ANNOVAR
perl $ANNOVAR_DIR/table_annovar.pl $OUTPUT_DIR/merge.vcf.gz $ANNOVAR_ANNO/humandb \
 -buildver hg38 \
 -out $OUTPUT_DIR/annovar \
 -remove \
 -protocol wgEncodeGencodeBasicV38 \
 -operation g \
 -nastring . \
 -vcfinput \
 -arg '-hgvs'

# Filter for exonic variants
$BCFPATH/bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV38~"exonic" | INFO/Func.wgEncodeGencodeBasicV38~"splicing"' -o $OUTPUT_DIR/filter_exon.vcf.gz -Oz $OUTPUT_DIR/annovar.hg38_multianno.vcf

# Annotate gnomAD frequencies
$VCFANNO_PATH/vcfanno_linux64.1 $GTOML $OUTPUT_DIR/filter_exon.vcf.gz | bgzip > $OUTPUT_DIR/gnomad.vcf.gz

# Remove variants with greater than 0.1% frequency in gnomAD_all
$BCFPATH/bcftools view -i 'ge_AF<=0.001 | ge_AF="."' $OUTPUT_DIR/gnomad.vcf.gz | $BCFPATH/bcftools view -i 'gg_AF<=0.001 | gg_AF="."' | bgzip > $OUTPUT_DIR/population_frequency.vcf.gz

# Annotate CADD
$VCFANNO_PATH/vcfanno_linux64.1 $CTOML $OUTPUT_DIR/population_frequency.vcf.gz | bgzip > $OUTPUT_DIR/cadd.vcf.gz

# Add sample annotations back to variants
# Filters:
# 1. Annotated with ANNOVAR
# 2. Exonic or splicing variant
# 3. Rare or not present in gnomad exome
# 4. Rare or not present in gnomad genome
# The last 3 should be redundant (loci were already filtered for those criteria, so any loci failing those should not have annotation information), but included just in case
while read SAMPLE
do
	LAST_DIGIT="${SAMPLE: -1}"
	INPUT_VCF=$GVCF_PATH/$LAST_DIGIT/$SAMPLE.gvcf.gz
	$BCFPATH/bcftools annotate -a $OUTPUT_DIR/cadd.vcf.gz -c INFO $INPUT_VCF | $BCFPATH/bcftools view -i "ANNOVAR_DATE!='.'" | $BCFPATH/bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV38~"exonic" | INFO/Func.wgEncodeGencodeBasicV38~"splicing"' | $BCFPATH/bcftools view -i 'ge_AF<=0.001 | ge_AF="."' | $BCFPATH/bcftools view -i 'gg_AF<=0.001 | gg_AF="."' | bgzip > $OUTPUT_DIR/annotated/$LAST_DIGIT/$SAMPLE.vcf.gz
	$TABIXPATH/tabix -p vcf $OUTPUT_DIR/annotated/$LAST_DIGIT/$SAMPLE.vcf.gz
done < $SAMPLE_LIST

