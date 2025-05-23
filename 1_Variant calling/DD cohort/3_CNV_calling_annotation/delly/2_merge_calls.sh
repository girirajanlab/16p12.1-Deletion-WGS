#!/bin/bash

# Merge delly SV calls

# Input and output files
DELLY_BCF_PATH=/path/to/sample/bcfs/ # Use the BCFs generated from script 1_run_delly.sh
OUTPUT_FILE=/path/to/output.bcf

# Paths to tools
DELLY=/path/to/delly_v0.8.2

# Merge calls from all samples
$DELLY merge -o $OUTPUT_FILE $DELLY_BCF_PATH/*.bcf
