#!/bin/bash

# Filter out loci that failed population-level filters for a family

# Input and output files
AUTO_VCF=/path/to/input/autosome.vcf.gz # Use the autosomal VCFs generated from script 1_run_GangSTR_dumpSTR.sh here
CHRX_VCF=/path/to/input/chrX.vcf.gz # Use the chrX VCFs generated from script 1_run_GangSTR_dumpSTR.sh here
EXCLUDE_AUTO_LOCI=/path/to/loci/to/exclude/autosome.bed # Use the output of script 3_exclude_loci_bed.sh
EXCLUDE_CHRX_LOCI=/path/to/loci/to/exclude/chrX.bed # Use the output of script 3_exclude_loci_bed.sh
OUTPUT_DIR=/path/to/output/directory

# Tools
VCFTOOLS=/path/to/vcftools
BCFPATH=/path/to/bcftools_v1.3.1
VCF_QUERY=/path/to/vcf_query

# Remove excluded loci from STR calls for family
$VCFTOOLS/vcftools --exclude-bed $EXCLUDE_AUTO_LOCI --gzvcf $AUTO_VCF --recode --recode-INFO-all --out $OUTPUT_DIR/autosome.vcf
cat $OUTPUT_DIR/autosome.vcf | vcf-sort | bgzip -c > $OUTPUT_DIR/autosome.vcf.gz
tabix -p vcf $OUTPUT_DIR/autosome.vcf.gz

$VCFTOOLS/vcftools --exclude-bed $EXCLUDE_CHRX_LOCI --gzvcf $CHRX_VCF --recode --recode-INFO-all --out $OUTPUT_DIR/chrX.vcf
cat $OUTPUT_DIR/chrX.vcf | vcf-sort | bgzip -c > $OUTPUT_DIR/chrX.vcf.gz
tabix -p vcf $OUTPUT_DIR/chrX.vcf.gz

# Get all non-reference calls
$BCFPATH/bcftools view -i 'GT[*]="alt"' $OUTPUT_DIR/autosome.vcf.gz | bgzip > $OUTPUT_DIR/autosome_alt.vcf.gz
tabix -p vcf $OUTPUT_DIR/autosome_alt.vcf.gz

$BCFPATH/bcftools view -i 'GT[*]="alt"' $OUTPUT_DIR/chrX.vcf.gz | bgzip > $OUTPUT_DIR/chrX_alt.vcf.gz
tabix -p vcf $OUTPUT_DIR/chrX_alt.vcf.gz

# Convert to BED format
$VCF_QUERY/vcf-query $OUTPUT_DIR/autosome_alt.vcf.gz -f '%CHROM\t%POS\t%INFO/END\t%INFO/PERIOD\t%INFO/RU\n'  > $OUTPUT_DIR/autosome_alt.bed
$VCF_QUERY/vcf-query $OUTPUT_DIR/chrX_alt.vcf.gz -f '%CHROM\t%POS\t%INFO/END\t%INFO/PERIOD\t%INFO/RU\n'  > $OUTPUT_DIR/chrX_alt.bed

