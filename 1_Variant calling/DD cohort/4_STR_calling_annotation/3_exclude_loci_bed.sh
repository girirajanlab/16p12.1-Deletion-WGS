#!/bin/bash

# Create a bed file of loci that pass population level filters

# Input and output files
AUTO_VCF=/path/to/autosomal.vcf # Use the autosomal VCF generated from script 2_run_mergeSTR_dumpSTR.sh here
CHRX_VCF=/path/to/chrX.vcf # Use the chrX VCF generated from script 2_run_mergeSTR_dumpSTR.sh here
CHRX_HW_VCF=/path/to/chrX_hw.vcf # Use the chrX VCF filtered for Hardy-Weinberg equilibrium generated from script 2_run_mergeSTR_dumpSTR.sh here
OUTPUT_DIR=/path/to/output/directory

# Tools 
VCF_QUERY=/path/to/vcf_query

# Autosomal loci
$VCF_QUERY/vcf-query $AUTO_VCF -f '%CHROM\t%POS\t%INFO/END\t%FILTER\n' |  awk '{if ($4 != "PASS") print $0 ;}' >  $OUTPUT_DIR/auto_loci.bed

# chrX loci
$VCF_QUERY/vcf-query $CHRX_VCF -f '%CHROM\t%POS\t%INFO/END\t%FILTER\n' |  awk '{if ($4 != "PASS") print $0 ;}' >  $OUTPUT_DIR/chrX_filter_loci.bed
$VCF_QUERY/vcf-query $CHRX_HW_VCF -f '%CHROM\t%POS\t%INFO/END\t%FILTER\n' |  awk '{if ($4 != "PASS") print $0 ;}' >  $OUTPUT_DIR/chrX_HW_loci.bed

# Combine the two chrX BED files
cat $OUTPUT_DIR/chrX_filter_loci.bed $OUTPUT_DIR/chrX_HW_loci.bed | sort > $OUTPUT_DIR/chrX_loci.bed
