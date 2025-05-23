#!/bin/bash

# Combine all the tables into a single file

# Input and output files
INPUT_DIR=/path/to/input/directory/tables # Use the directory tables were saved into in script 4_filter_pop_frequency.sh
OUTPUT_TABLE=/path/to/output/table.tsv

cat $INPUT_DIR/*.tsv > $OUTPUT_TABLE
