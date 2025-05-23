#!/bin/bash

# This script does a few reformatting and filtering steps

# Input and output files
INPUT_TABLE=/path/to/input/table.txt # Use the table generated from script 6_vcf2table.sh here
OUTPUT_TABLE=/path/to/output/table.txt
TMP_FILE=/path/to/temporary/file

# 1. Add a header line
echo -e 'Sample\tChrom\tPos\tRef\tAlt\tQual\tFunc.wgEncodeGencodeBasicV19\tGene.wgEncodeGencodeBasicV19\tGeneDetail.wgEncodeGencodeBasicV19\tgg_AF\tgg_AF_nfe\tgg_AF_asj\tgg_AF_eas\tgg_AF_amr\tgg_AF_afr\tGT\tDP\tAD\tGQ\tPL\tRGQ\tSB' > $TMP_FILE

# 2. Filter out chrM
grep -v 'chrM' $INPUT_TABLE>> $TMP_FILE

# 3. Replace \x3b with ;
# And replace \x3d with :
sed -i 's/\\x3b/;/g' $TMP_FILE
sed -i 's/\\x3d/:/g' $TMP_FILE

# Save
cp $TMP_FILE $OUTPUT_TABLE
rm $TMP_FILE
