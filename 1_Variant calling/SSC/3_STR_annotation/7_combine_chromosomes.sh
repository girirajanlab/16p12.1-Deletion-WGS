#!/bin/bash

# Concatenate all chromosome tables into a single file

# Input and output files
INPUT_FILE_LIST=/path/to/list/of/input/filenames.list
OUTPUT_TABLE=/path/to/output/file.tsv

# Concat files
while read FILE
do
	cat $FILE >> $OUTPUT_TABLE
done < $INPUT_FILE_LIST
