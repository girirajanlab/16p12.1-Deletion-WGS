#!/bin/bash

# Use ANNOVAR to annotate STR variants with affected genes

# Input and output files
INPUT_TAB=/path/to/expansion/table.tsv # Use the output from script 7_2SD_STR.py
OUTPUT_DIR=/path/to/output/directory

# Paths to annovar scripts and annotations
ANNOVAR_DIR=/path/to/annovar
ANNOVAR_ANNO=/path/to/annovar/annotations

# Format table for ANNOVAR annotation
tail +2 $INPUT_TAB | cut -f1-5 | sort | uniq > $OUTPUT_DIR/annovar_input.bed

# Run ANNOVAR
perl $ANNOVAR_DIR/table_annovar.pl $OUTPUT_DIR/annovar_input.bed $ANNOVAR_ANNO/annovar \
 -buildver hg19 \
 -out $OUTPUT_DIR/STR_annovar \
 -remove \
 -protocol refGene,wgEncodeGencodeBasicV19 \
 -operation g,g \
 -nastring . \
 -arg '-hgvs','-hgvs'
