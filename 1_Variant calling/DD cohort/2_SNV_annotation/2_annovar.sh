#!/bin/bash

# Use ANNOVAR to annotate variants with affected genes

# Paths to annovar scripts and annotations
ANNOVAR_DIR=/path/to/annovar
ANNOVAR_ANNO=/path/to/annovar/annotations

# Input and output files
INPUT_VCF=/path/to/input.vcf.gz # Use the output of script 1_filter_quality.sh here
OUTPUT_PREFIX=/path/to/output

# Annotate variant effects
perl $ANNOVAR_DIR/table_annovar.pl $INPUT_VCF $ANNOVAR_ANNO/humandb \
 -buildver hg19 \
 -out $OUTPUT_PREFIX \
 -remove \
 -protocol refGene,wgEncodeGencodeBasicV19 \
 -operation g,g \
 -nastring . \
 -vcfinput \
 -arg '-hgvs','-hgvs'

