#!/bin/bash

# Annotate variants in a chromsome with functional consequences and delteriousness metrics

# Input and output files
CHROM=chromsome # i.e. 1-22, X, and Y
INPUT_VCF=/path/to/SFARI/SSC/WGS2/Project_SSC_9209samples.JGvariants.2019-06-21/CCDG_9000JG_B01_GRM_WGS_2019-03-21_chr$CHROM.recalibrated_variants.flagged.vcf.gz
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
$BCFPATH/bcftools norm -f $REF -m-both $INPUT_VCF | bcftools view -i 'GT="alt" & QUAL>=50 & FORMAT/DP>=8 & FORMAT/AD[:1]>0 & (FORMAT/AD[:1])/(FORMAT/DP)>=0.25 & QUAL/(FORMAT/AD[:1])>=1.5' | bcftools view -i '(FORMAT/AD[:1])/(FORMAT/DP)<=0.75 | (FORMAT/AD[:1])/(FORMAT/DP)>=0.9' -o $OUTPUT_DIR/$CHROM.quality_filter.vcf.gz -Oz

# Strip the samples from the VCF
$BCFPATH/bcftools view -G -o $OUTPUT_DIR/$CHROM.samples_stripped.vcf.gz -Oz $OUTPUT_DIR/$CHROM.quality_filter.vcf.gz

# Annotate variant effects with ANNOVAR
perl $ANNOVAR_DIR/table_annovar.pl $OUTPUT_DIR/$CHROM.samples_stripped.vcf.gz $ANNOVAR_ANNO/humandb \
 -buildver hg38 \
 -out $OUTPUT_DIR/$CHROM_annovar \
 -remove \
 -protocol wgEncodeGencodeBasicV38 \
 -operation g \
 -nastring . \
 -vcfinput \
 -arg '-hgvs'

# Filter for exonic variants
$BCFPATH/bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV38~"exonic" | INFO/Func.wgEncodeGencodeBasicV38~"splicing"' -o $OUTPUT_DIR/$CHROM.filter_exon.vcf.gz -Oz $OUTPUT_DIR/$CHROM_annovar.hg38_multianno.vcf

# Annotate gnomAD frequencies
$VCFANNO_PATH/vcfanno_linux64.1 $GTOML $OUTPUT_DIR/$CHROM.filter_exon.vcf.gz | bgzip > $OUTPUT_DIR/$CHROM.gnomad.vcf.gz

# Remove variants with greater than 0.1% frequency in gnomAD_all
$BCFPATH/bcftools view -i 'ge_AF<=0.001 | ge_AF="."' $OUTPUT_DIR/$CHROM.gnomad.vcf.gz | bcftools view -i 'gg_AF<=0.001 | gg_AF="."' | bgzip > $OUTPUT_DIR/$CHROM.population_frequency.vcf.gz

# Annotate CADD
$VCFANNO_PATH/vcfanno_linux64.1 $CTOML $OUTPUT_DIR/$CHROM.population_frequency.vcf.gz | bgzip > $OUTPUT_DIR/$CHROM.cadd.vcf.gz

# Add sample annotations back to variants
$BCFPATH/bcftools index $OUTPUT_DIR/$CHROM.cadd.vcf.gz
$BCFPATH/bcftools isec $OUTPUT_DIR/$CHROM.quality_filter.vcf.gz $OUTPUT_DIR/$CHROM.cadd.vcf.gz -p $OUTPUT_DIR/$CHROM -n =2 -w 1
bgzip $OUTPUT_DIR/$CHROM/0000.vcf
$BCFPATH/bcftools index $OUTPUT_DIR/$CHROM/0000.vcf.gz

$BCFPATH/bcftools merge $OUTPUT_DIR/$CHROM/0000.vcf.gz $OUTPUT_DIR/$CHROM.cadd.vcf.gz | bgzip > $OUTPUT_DIR/$CHROM.merged.vcf.gz
$TABIXPATH/tabix -p vcf $OUTPUT_DIR/$CHROM.merged.vcf.gz

# Split population VCF into individual VCFs for ancestry filtering
java -Xmx50G -jar $JVARKIT/dist/biostar130456.jar -x -z -p "$OUTPUT_DIR/$CHROM/sample_split/__SAMPLE__.vcf" $OUTPUT_DIR/$CHROM.merged.vcf.gz

