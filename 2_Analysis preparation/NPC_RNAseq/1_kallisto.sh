#!/bin/bash

# Input and Output files
FWD=/path/to/R1.fastq.gz
REV=/path/to/R2.fastq.gz
ADAPTER_FASTA=/path/to/adapter/sequence.fa # File derived from adapter sequences available here: https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
REF=/path/to/GRCh38/reference.fa

OUTPUT_DIR="/path/to/output/files"
NAME="file_prefix"
KIDX_OUTPUT=/path/to/output/kallisto/index.kidx

# Run Trimmomatic
TRIM_R1=${OUTPUT_DIR}/${NAME}_R1_trimmed.fastq.gz
TRIM_R2=${OUTPUT_DIR}/${NAME}_R2_trimmed.fastq.gz
UNPAIR_R1=${OUTPUT_DIR}/${NAME}_R1_unpaired.fastq.gz
UNPAIR_R2=${OUTPUT_DIR}/${NAME}_R2_unpaired.fastq.gz

trimmomatic PE -phred33 -threads 20 $FWD $REV $TRIM_R1 $UNPAIR_R1 $TRIM_R2 $UNPAIR_R2 ILLUMINACLIP:${ADAPTER_FASTA}:2:30:10:7:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Build kallisto index
kallisto index-i $KIDX_OUTPUT $REF

# Run kallisto
kallisto quant -i $KIDX_OUTPUT -o $OUTPUT_DIR -b 100 -t 5 $TRIM_R1 $TRIM_R2

# Before next step, compile all of the abundance values across samples into a single file
