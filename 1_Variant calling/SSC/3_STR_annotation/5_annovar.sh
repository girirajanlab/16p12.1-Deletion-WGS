#!/bin/bash

# Use ANNOVAR to annotate STR variants with affected genes

# Input and output files
INPUT_TAB=/path/to/expansion/table.tsv # Use the annovar input generated in script 4_anno_locus_info.py
OUTPUT_DIR=/path/to/output/directory

# Paths to annovar scripts and annotations
ANNOVAR_DIR=/path/to/annovar
ANNOVAR_ANNO=/path/to/annovar/annotations

# Run ANNOVAR
perl $ANNOVAR_DIR/table_annovar.pl $OUTPUT_DIR/annovar_input.bed $ANNOVAR_ANNO/annovar \
 -buildver hg38 \
 -out $OUTPUT_DIR/STR_annovar \
 -remove \
 -protocol wgEncodeGencodeBasicV38 \
 -operation g \
 -nastring . \
 -arg '-hgvs'


