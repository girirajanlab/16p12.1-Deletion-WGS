#!/bin/bash

# Separate calls from each caller by sample

# Input and output files
SAMPLE_LIST=/path/to/list/of/samples.csv
OUTPUT_PATH=/path/to/output/files

# CNV caller outputs
CNVNATOR=/path/to/cnvnator/small/cnv/output.bed # Use the output from script cnvnator/small_cnv_processing/5_str_filter.sh
MANTA=/path/to/manta/output.bed # Use the output from script manta/6_str_filter.sh
DELLY=/path/to/delly/output.bed # Use the output from script delly/10_str_filter.sh
LUMPY=/path/to/lumpy/output.bed # Use the output from script lumpy/10_str_filter.sh

while read SAMPLE
do
	grep $SAMPLE $CNVNATOR > $OUTPUT_PATH/$SAMPLE.cnvnator.bed
	grep $SAMPLE $MANTA > $OUTPUT_PATH/$SAMPLE.manta.bed
	grep $SAMPLE $DELLY > $OUTPUT_PATH/$SAMPLE.delly.bed
	grep $SAMPLE $LUMPY > $OUTPUT_PATH/$SAMPLE.lumpy.bed
done < $SAMPLE_LIST
