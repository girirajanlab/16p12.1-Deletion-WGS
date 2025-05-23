#!/bin/bash

# Annotate CNV calls with their overlaps with segmental duplications and centromeres and telomeres

# Input and output files
INPUT=/path/to/input/file.txt # Use the output of 4_control_frequency.py here
OUTPUT_DIR=/path/to/output/files

# Reference files
SEGDUP=/path/to/segmental/duplication/reference.bed
CENTRO_TELO=/path/to/centromere/and/telomere/reference.bed

# Calculate the overlaps of segmental duplications and cetromeres and telomeres
perl scripts/CalculateSDOverlaps.pl $INPUT $SEGDUP $OUTPUT_DIR/segdup.txt
perl scripts/CalculateSDOverlaps.pl $OUTPUT_DIR/segdup.txt $CENTRO_TELO $OUTPUT_DIR/segdup_centro_telo.txt

# Integer overlaps output by this command were manually converted into percentages

