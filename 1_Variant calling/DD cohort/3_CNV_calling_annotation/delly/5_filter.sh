#!/bin/bash

# Filter delly calls

# Input and output files
INPUT_FILE=/path/to/egenotyped/merged.bcf # Use the BCF generated from script 4_merge_regeno.sh
OUTPUT_PATH=/path/to/output/files
SAMPLE_LIST=/path/to/list/of/samples.txt # a file listing all samples processed, one sample per line

# Paths to tools
BCFPATH=/path/to/bcftools_v1.12
DELLY=/path/to/delly_v0.8.2

# Filter for germline variants with delly
$DELLY filter -f germline -o $OUTPUT_PATH/delly_filter.bcf $INPUT_FILE

# Perform filtering with locus-specific (FILTER==PASS) and sample-specific (FT==PASS) default filters
$BCFPATH/bcftools view -i 'FILTER="PASS" & FORMAT/FT="PASS" & FORMAT/GT="alt"' $OUTPUT_PATH/delly_filter.bcf | bcftools filter --set-GTs . -i 'FORMAT/FT="PASS" & FORMAT/GT="alt"' |  bgzip > $OUTPUT_PATH/filter2.vcf.gz

# Split into individual-sample BCFs
while read SAMPLE
do
	$BCFPATH/bcftools view -s $SAMPLE $OUTPUT_PATH/filter2.vcf.gz | bcftools view -i 'FORMAT/GT="alt"' > $OUTPUT_PATH/$SAMPLE.vcf
done

