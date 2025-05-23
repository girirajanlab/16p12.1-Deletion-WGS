#!/bin/bash

# Run PennCNV to call CNVs from microarray

# Input and output files
OMNI1_PATH=/path/to/input/OMNI1/final/report/files # Path to Illumina final report files split by sample
OMNI2_PATH=/path/to/input/OMNI2/final/report/files # Path to Illumina final report files split by sample
# Note that there are two input arrays and CNVs were called separately on each array
# These input files can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
OMNI1_TRIO=/path/to/list/of/OMNI1/trio/samples.txt
OMNI2_TRIO=/path/to/list/of/OMNI2/trio/samples.txt
OMNI1_QUAD=/path/to/list/of/OMNI1/quartet/samples.txt
OMNI2_QUAD=/path/to/list/of/OMNI2/quartet/samples.txt
SEXFILE=/path/to/sample/sexes.list # A list of sample sexes for chrX CNV calling
OUTPUT_DIR=/path/to/output/directory

# Tools
PENNCNV_PATH=/path/to/PennCNV-1.0.5

# Compile PFB
$PENNCNV_PATH/compile_pfb.pl `ls $OMNI1_PATH/*.finalReport.txt` --output $OUTPUT_DIR/OMNI1.pfb
$PENNCNV_PATH/compile_pfb.pl `ls $OMNI2_PATH/*.finalReport.txt` --output $OUTPUT_DIR/OMNI2.pfb

# Run PennCNV
for BATCH in OMNI1 OMNI2
do
	if [ $BATCH == 'OMNI1' ]
	then
		FINALREPORT_PATH=$OMNI1_PATH
		TRIO_LIST=$OMNI1_TRIO
		QUAD_LIST=$OMNI1_QUAD
	else
		FINALREPORT_PATH=$OMNI2_PATH
		TRIO_LIST=$OMNI2_TRIO
		QUAD_LIST=$OMNI2_QUAD
	fi

	# Call CNVs per individual
	perl $PENNCNV_PATH/detect_cnv.pl \
		-test \
		-hmm $PENNCNV_PATH/lib/hhall.hmm \
		-pfb $OUTPUT_DIR/$BATCH.pfb \
		$FINALREPORT_PATH/*.finalReport.txt \
		-log $OUTPUT_DIR/logs/${BATCH}_autosome.log \
		-out $OUTPUT_DIR/penncnv_output/${BATCH}_pcnv_autosome.txt

	perl $PENNCNV_PATH/detect_cnv.pl \
		-test --chrX \
		-hmm $PENNCNV_PATH/lib/hhall.hmm \
		-pfb $OUTPUT_DIR/$BATCH.pfb \
		$FINALREPORT_PATH/*.finalReport.txt \
		--sexfile $SEXFILE \
		-log $OUTPUT_DIR/logs/${BATCH}_chrX.log \
		-out $OUTPUT_DIR/penncnv_output/${BATCH}_pcnv_chrX.txt

	# Call CNVs per family
	perl $PENNCNV_PATH/detect_cnv.pl \
		-trio -hmm $PENNCNV_PATH/lib/hhall.hmm \
		-pfb $OUTPUT_DIR/$BATCH.pfb \
		--listfile $TRIO_LIST \
		--out $OUTPUT_DIR/penncnv_output/${BATCH}_trios.txt
	
	perl $PENNCNV_PATH/detect_cnv.pl \
		-quartet -hmm $PENNCNV_PATH/lib/hhall.hmm \
		-pfb $OUTPUT_DIR/$BATCH.pfb \
		--listfile $QUAD_LIST \
		--out $OUTPUT_DIR/penncnv_output/${BATCH}_quads.txt
done

# Any samples that did not pass default quality control criteria in PennCNV were manually removed for downstream analysis using warnings provided in the PennCNV log files

