#!/bin/bash

# Run PennCNV to call CNVs from microarray data

# Input and output files
INPUT_FINALREPORT=/path/to/Illumina/FinalReport.txt
OUTPUT_DIR=/path/to/output/directory
SNP_POS=/path/to/array/SNP_list.csv
SEX_FILE=/path/to/file/with/sample/sexes.tsv
TRIOS=/path/to/file/with/sample/trios.txt # File is in the format Father Mother Proband
QUADS=/path/to/file/with/sample/quads.txt # File is in the format Father Mother Proband Sibling
CHRX_TRIO=/path/to/file/with/chrX/sample/trios.txt # For chrX calling, split all quads into trios

# Paths to tools
$PENNCNV_PATH=/path/to/PennCNV-1.0.5

# Prepare signal intensity files from Illumina Report
$PENNCNV_PATH/split_illumina_report.pl -prefix $OUTPUT_DIR/data -suffix .txt $INPUT_FINALREPORT

# Create a list of SNPs from the array
cut -d ',' -f 2,10,11 $SNP_list > $OUTPUT_DIR/SNP_Positions.csv
sed 's/,/\t/g' $OUTPUT_DIR/SNP_Positions.csv > $OUTPUT_DIR/SNP_Positions.tsv

# Build a PFB
$PENNCNV_PATH/compile_pfb.pl `ls $OUTPUT_DIR/data/*.txt` --snpposfile $OUTPUT_DIR/SNP_Positions.tsv --output Array.pfb

# Call CNVs individually
perl $PENNCNV_PATH/detect_cnv.pl -test -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb Array.pfb $OUTPUT_DIR/data/*.txt -log penncnv_autosome.log -out penncnv_autosome.txt
perl $PENNCNV_PATH/detect_cnv.pl -test --chrx -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb Array.pfb $OUTPUT_DIR/data/*.txt --sexfile $SEX_FILE -log penncnv_chrX.log -out penncnv_chrX.tsv

# Call CNVs by family
perl $PENNCNV_PATH/detect_cnv.pl -trio -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb Array.pfb -cnv penncnv_autosome.txt --listfile $TRIOS -out penncnv_trios.txt
perl $PENNCNV_PATH/detect_cnv.pl -quartet -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb Array.pfb -cnv penncnv_autosome.txt --listfile $QUADS -out penncnv_quads.txt
perl $PENNCNV_PATH/detect_cnv.pl -trio --chrx -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb Array.pfb -cnv penncnv_chrX.tsv --listfile $CHRX_TRIO--sexfile $SEX_FILE -out penncnv_chrX_trios.txt

# Cat autosome and chrX calls together
cat penncnv_autosome.txt > penncnv_all.txt
cat penncnv_chrX.tsv >> penncnv_all.txt

cat penncnv_trios.txt > penncnv_family.txt
cat penncnv_quads.txt >> penncnv_family.txt
cat penncnv_chrX_trios.txt >> penncnv_family.txt

# Any samples that did not pass default quality control criteria in PennCNV were manually removed for downstream analysis using warnings provided in the PennCNV log files

