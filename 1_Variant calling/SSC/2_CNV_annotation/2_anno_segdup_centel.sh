#!/bin/bash

# Annotate CNV calls with their overlaps with segmental duplications and centromeres and telomeres

# Input and output files
INPUT=/path/to/input/file.txt # Use the output of 1_filter_CNVs.py
OUTPUT_DIR=/path/to/output/files

# Reference files
SEGDUP=/path/to/segmental/duplication/reference.bed
CENTRO_TELO=/path/to/centromere/and/telomere/reference.bed

# Calculate the overlaps of segmental duplications and cetromeres and telomeres
perl scripts/CalculateSDOverlaps.pl $INPUT $SEGDUP $OUTPUT_DIR/segdup.txt
perl scripts/CalculateSDOverlaps.pl $OUTPUT_DIR/segdup.txt $CENTRO_TELO $OUTPUT_DIR/segdup_centro_telo.txt

# Integer outputs were manually converted to percentages
# Calls were manually filtered to remove CNVs with greater than 75% overlap with segmental duplications and centromere/telomere sequences
