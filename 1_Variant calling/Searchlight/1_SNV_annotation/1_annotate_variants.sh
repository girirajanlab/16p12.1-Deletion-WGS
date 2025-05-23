#!/bin/bash

# Annotate variants in a chromsome with functional consequences and delteriousness metrics

# Input and output files
INPUT_VCF=/path/to/SFARI/Searchlight/Daly/Searchlight_exomes.vcf
# This input VCF can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
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
JVARKIT=/path/to/jvarkit

# Add chr annotation to calls and fix AD definition within header
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $INPUT_VCF > $OUTPUT_DIR/ano_fix.vcf
sed -i "s|Number=\.|Number=R|g" $OUTPUT_DIR/ano_fix.vcf
bgzip $OUTPUT_DIR/ano_fix.vcf
tabix $OUTPUT_DIR/ano_fix.vcf.gz

# Split multiallelic records and left-normalize variants
# And filter for:
# 1. Non-ref alleles
# 2. QUAL>=50
# 3. Read depth>=8
# 4. Alternative allele depth>0
# 5. Allele balance >=0.25 and <=0.75 or >=0.9
# 6. QUAL/(Alternative allele depth)>=1.5
$BCFPATH/bcftools norm -f $REF -m-both $OUTPUT_DIR/ano_fix.vcf.gz | bcftools view -i 'GT="alt" & QUAL>=50 & FORMAT/DP>=8 & FORMAT/AD[:1]>0 & (FORMAT/AD[:1])/(FORMAT/DP)>=0.25 & QUAL/(FORMAT/AD[:1])>=1.5' | bcftools view -i '(FORMAT/AD[:1])/(FORMAT/DP)<=0.75 | (FORMAT/AD[:1])/(FORMAT/DP)>=0.9' -o $OUTPUT_DIR/quality_filter.vcf.gz -Oz

# Strip the samples from the VCF
$BCFPATH/bcftools view -G -o $OUTPUT_DIR/samples_stripped.vcf.gz -Oz $OUTPUT_DIR/quality_filter.vcf.gz

# Annotate variant effects with ANNOVAR
perl $ANNOVAR_DIR/table_annovar.pl $OUTPUT_DIR/samples_stripped.vcf.gz $ANNOVAR_ANNO/humandb \
 -buildver hg19 \
 -out $OUTPUT_DIR/annovar \
 -remove \
 -protocol wgEncodeGencodeBasicV19 \
 -operation g \
 -nastring . \
 -vcfinput \
 -arg '-hgvs'

# Filter for exonic variants
$BCFPATH/bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV19~"exonic" | INFO/Func.wgEncodeGencodeBasicV19~"splicing"' -o $OUTPUT_DIR/filter_exon.vcf.gz -Oz $OUTPUT_DIR/annovar.hg19_multianno.vcf

# Annotate gnomAD frequencies
$VCFANNO_PATH/vcfanno_linux64.1 $GTOML $OUTPUT_DIR/filter_exon.vcf.gz | bgzip > $OUTPUT_DIR/gnomad.vcf.gz

# Remove variants with greater than 0.1% frequency in gnomAD_all
$BCFPATH/bcftools view -i 'ge_AF<=0.001 | ge_AF="."' $OUTPUT_DIR/gnomad.vcf.gz | bcftools view -i 'gg_AF<=0.001 | gg_AF="."' | bgzip > $OUTPUT_DIR/population_frequency.vcf.gz

# Annotate CADD
$VCFANNO_PATH/vcfanno_linux64.1 $CTOML $OUTPUT_DIR/population_frequency.vcf.gz | bgzip > $OUTPUT_DIR/cadd.vcf.gz

# Add sample annotations back to variants
$BCFPATH/bcftools index $OUTPUT_DIR/cadd.vcf.gz
$BCFPATH/bcftools isec $OUTPUT_DIR/quality_filter.vcf.gz $OUTPUT_DIR/cadd.vcf.gz -p $OUTPUT_DIR/intersect -n =2 -w 1
bgzip $OUTPUT_DIR/intersect/0000.vcf
$BCFPATH/bcftools index $OUTPUT_DIR/intersect/0000.vcf.gz

$BCFPATH/bcftools merge $OUTPUT_DIR/intersect/0000.vcf.gz $OUTPUT_DIR/cadd.vcf.gz | bgzip > $OUTPUT_DIR/merged.vcf.gz
$TABIXPATH/tabix -p vcf $OUTPUT_DIR/merged.vcf.gz

# Split population VCF into individual VCFs for ancestry filtering
java -Xmx50G -jar $JVARKIT/dist/biostar130456.jar -x -z -p "$OUTPUT_DIR/sample_split/__SAMPLE__.vcf" $OUTPUT_DIR/merged.vcf.gz

