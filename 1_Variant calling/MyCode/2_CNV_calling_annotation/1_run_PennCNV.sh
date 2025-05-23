#!/bin/bash

# Run PennCNV to call CNVs from microarray

# Note that MyCode samples were run on two arrays and CNV calling is run separately for each array

# Input and output files
OMNI_PATH=/path/to/input/signal/files/for/60K_OMNI/Array # Path to array signal files for MyCode samples genotyped on Omni array
GSA_PATH=/path/to/input/signal/files/for/145K_GSA/Array # Path to array signal files for MyCode samples genotyped on GSA array
# These input files were provided by MyCode
SNP_POS=/path/to/SNP/Position/file.txt
SEXFILE=/path/to/sample/sexes.list # A list of sample sexes for chrX CNV calling
OUTPUT_DIR=/path/to/output/directory

# Tools
PENNCNV_PATH=/path/to/PennCNV-1.0.5

# Create a list of signal files
find $OMNI_PATH/ -name "*.signal.gz" > $OUTPUT_DIR/OMNI_signal.list
find $GSA_PATH/ -name "*.signal.gz" > $OUTPUT_DIR/GSA_signal.list

# Compile PFB
$PENNCNV_PATH/compile_pfb.pl --listfile $OUTPUT_DIR/OMNI_signal.list --snpposfile $SNP_POS --output $OUTPUT_DIR/OMNI.pfb
$PENNCNV_PATH/compile_pfb.pl --listfile $OUTPUT_DIR/GSA_signal.list --snpposfile $SNP_POS --output $OUTPUT_DIR/GSA.pfb

# Make a GC model adjustment file
mkdir $OUTPUT_DIR/GC_files
gunzip -c $PENNCNV_PATH/gc_file/hg38.gc5Base.txt.gz > $OUTPUT_DIR/GC_files/hg38.gc5Base.txt
$PENNCNV_PATH/cal_gc_snp.pl $OUTPUT_DIR/GC_files/hg38.gc5Base.txt $OUTPUT_DIR/OMNI.pfb > $OUTPUT_DIR/GC_files/OMNI.gcmodel
$PENNCNV_PATH/cal_gc_snp.pl $OUTPUT_DIR/GC_files/hg38.gc5Base.txt $OUTPUT_DIR/GSA.pfb > $OUTPUT_DIR/GC_files/GSA.gcmodel

# Run PennCNV
mkdir -p $OUTPUT_DIR/logs
mkdir -p $OUTPUT_DIR/penncnv_output
for BATCH in OMNI GSA
do
	perl $PENNCNV_PATH/detect_cnv.pl \
		-test \
		-hmm $PENNCNV_PATH/lib/hhall.hmm \
		-pfb $OUTPUT_DIR/$BATCH.pfb \
		--list $OUTPUT_DIR/${BATCH}_signal.list \
		-log $OUTPUT_DIR/logs/${BATCH}_autosome.log \
		-out $OUTPUT_DIR/penncnv_output/${BATCH}_pcnv_autosome.txt \
		--gcmodel $OUTPUT_DIR/GC_files/$BATCH.gcmodel
	
	perl $PENNCNV_PATH/detect_cnv.pl \
		-test --chrX \
		-hmm $PENNCNV_PATH/lib/hhall.hmm \
		-pfb $OUTPUT_DIR/$BATCH.pfb \
		--list $OUTPUT_DIR/${BATCH}_signal.list \
		--sexfile $SEXFILE \
		-log $OUTPUT_DIR/logs/${BATCH}_chrX.log \
		-out $OUTPUT_DIR/penncnv_output/${BATCH}_pcnv_chrX.txt
done

# Perform QC of PennCNV calls
mkdir -p $OUTPUT_DIR/QC/logs
mkdir -p $OUTPUT_DIR/QC/qc_pass_samples
mkdir -p $OUTPUT_DIR/QC/qc_sum
mkdir -p $OUTPUT_DIR/QC/goodcnv
for BATCH in OMNI GSA
do
	$PENNCNV_PATH/filter_cnv.pl $OUTPUT_DIR/penncnv_output/${BATCH}_pcnv_autosome.txt \
		--qclogfile $OUTPUT_DIR/QC/logs/autosome_$BATCH.log \
		--qclrrsd 0.35 \
		--qcpassout $OUTPUT_DIR/QC/qc_pass_samples/autosome_$BATCH.qcpass \
		--qcsumout $OUTPUT_DIR/QC/qc_sum/autosome_$BATCH.qcsum \
		--out $OUTPUT_DIR/QC/goodcnv/autosome_$BATCH.goodcnv
	
	$PENNCNV_PATH/filter_cnv.pl $OUTPUT_DIR/penncnv_output/${BATCH}_pcnv_chrX.txt \
		--qclogfile $OUTPUT_DIR/QC/logs/chrX_$BATCH.log \
		--qclrrsd 0.35 \
		--chrX \
		--qcpassout $OUTPUT_DIR/QC/qc_pass_samples/chrX_$BATCH.qcpass \
		--qcsumout $OUTPUT_DIR/QC/qc_sum/chrX_$BATCH.qcsum \
		--out $OUTPUT_DIR/QC/goodcnv/chrX_$BATCH.goodcnv
done
