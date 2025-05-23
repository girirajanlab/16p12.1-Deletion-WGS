#!/bin/bash

# Run PennCNV to call CNVs from microarray

# Input and output files
LRR_DIR=/path/to/LRR/input/txt/files # Path to LRR values for UK Biobank
BAF_DIR=/path/to/LRR/input/txt/files # Path to BAF values for UK Biobank
BIM_DIR=/path/to/input/bim/files # Path to array BIM file
FAM_FILE=/path/to/UKB/FAM/file.fam # Path to array FAM file
SAMPLE_LIST=/path/to/list/of/sample_ids.txt # A list of sample IDs
SEXFILE=/path/to/sample/sexes.list # A list of sample sexes for chrX CNV calling
OUTPUT_DIR=/path/to/output/directory

# Tools
PENNCNV_PATH=/path/to/PennCNV-1.0.5

# Split files
# Also create a file of SNP names and positions
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	lines=$(wc -l $INPUT_DIR/ukb_snp_chr${chrom}_v2.bim | cut -f 1 -d ' ')
	num_lines=$((($lines+2)/3))

	# Split LRR
	split -l $num_lines -d $INPUT_DIR/ukb_l2r_chr${chrom}_v2.txt $OUTPUT_DIR/split_files/l2r_chr${chrom}_
	# Split BAF
	split -l $num_lines -d $INPUT_DIR/ukb_baf_chr${chrom}_v2.txt $OUTPUT_DIR/split_files/baf_chr${chrom}_

	# Save SNP names
	cut -f 2 $INPUT_DIR/ukb_snp_chr${chrom}_v2.bim > tmp_name
	cut -f 1 $INPUT_DIR/ukb_snp_chr${chrom}_v2.bim > tmp_chr
	cut -f 4 $INPUT_DIR/ukb_snp_chr${chrom}_v2.bim > tmp_pos
	paste tmp_name tmp_chr tmp_pos >> $OUTPUT_DIR/UKB_SNP_Positions.list
	sed -i 's/\t23\t/\tX\t/' $OUTPUT_DIR/UKB_SNP_Positions.list
	rm tmp_*
done

# Split columns in files to create separate signal files for each sample
find $OUTPUT_DIR/split_files/ -name "*_chr*" > $OUTPUT_DIR/split_files.list
while read SPLIT
do
	Rscript helper_scripts/subset_cols.R $SPLIT $FAM_FILE $OUTPUT_DIR
done < $OUTPUT_DIR/split_files.list

# Concat split file by sample
cut -f 1 $OUTPUT_DIR/UKB_SNP_Positions.list > tmp_names.list
while read SAMPLE
do
	BAF_OUTPUT=$OUTPUT_DIR/BAF_LRR_files/$SAMPLE/${SAMPLE}_baf.txt
	LRR_OUTPUT=$OUTPUT_DIR/BAF_LRR_files/$SAMPLE/${SAMPLE}_lrr.txt
	for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
	do
		cat $OUTPUT_DIR/sample_split/$SAMPLE/${SAMPLE}_baf_chr${chrom}_00.txt > $BAF_OUTPUT
		cat $OUTPUT_DIR/sample_split/$SAMPLE/${SAMPLE}_baf_chr${chrom}_01.txt >> $BAF_OUTPUT
		cat $OUTPUT_DIR/sample_split/$SAMPLE/${SAMPLE}_baf_chr${chrom}_02.txt >> $BAF_OUTPUT

		cat $OUTPUT_DIR/sample_split/$SAMPLE/${SAMPLE}_lrr_chr${chrom}_00.txt > $LRR_OUTPUT
		cat $OUTPUT_DIR/sample_split/$SAMPLE/${SAMPLE}_lrr_chr${chrom}_01.txt >> $LRR_OUTPUT
		cat $OUTPUT_DIR/sample_split/$SAMPLE/${SAMPLE}_lrr_chr${chrom}_02.txt >> $LRR_OUTPUT
	done

	paste tmp_names.list $BAF_OUTPUT $LRR_OUTPUT | gzip -c > $OUTPUT_DIR/signal_files/${SAMPLE}_input_signal.txt.gz

done < $SAMPLE_LIST
rm tmp_names.list

# Save a list of input files
find $OUTPUT_DIR/signal_files/ -name "*_input_signal.txt.gz" > $OUTPUT_DIR/input_signal.list

# PennCNV can use use system commands in addition to strict filenames
# So add a command to unzip the files before piping them into PennCNV
sed -i 's/^/`gunzip -c /' $OUTPUT_DIR/input_signal.list
sed -i 's/$/`/' $OUTPUT_DIR/input_signal.list

# Make a PFB file
# PennCNV cannot use gzipped files to make a PFB file (although it can use them for CNV calling)
# Pick 1000 files, unzip them, and use them as the input to make a PFB
mkdir -p $OUTPUT_DIR/pfb_sample_files
shuf -f 1000 $OUTPUT_DIR/input_signal.list > $OUTPUT_DIR/tmp_pfb_input_files.list
for i in {1..1000}
do
	FILE=`head -n $i $OUTPUT_DIR/tmp_pfb_input_files.list | tail -1`
	SAMPLE=`head -n $i head -n $i $OUTPUT_DIR/tmp_pfb_input_files.list | tail -1 | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '_'`

	echo -e "Name\t"$SAMPLE".B Allele Freq" > $OUTPUT_DIR/pfb_sample_files/$SAMPLE'_pfb_input.txt'
	gunzip -c $FILE | tail -n +2 | awk '{ print $1 "\t" $4 }' >> $OUTPUT_DIR/pfb_sample_files/$SAMPLE'_pfb_input.txt'
done
rm $OUTPUT_DIR/tmp_pfb_input_files.list
find $OUTPUT_DIR/pfb_sample_files/ -name -name "*_pfb_input.txt" > $OUTPUT_DIR/pfb_input_files.list

# Compile PFB
$PENNCNV_PATH/compile_pfb.pl --listfile $OUTPUT_DIR/pfb_input_files.list --snpposfile $OUTPUT_DIR/UKB_SNP_Positions.list --output $OUTPUT_DIR/UKB.pfb

# Create lists of 1000 files to run concurrently in PennCNV
mkdir -p $OUTPUT_DIR/pcnv_input_lists
total_lines=`wc -l $OUTPUT_DIR/input_signal.list | cut -f 1 -d ' '`
num_batches=$(( $(($total_lines/1000))+1 ))
for (( i = 0; i <= $num_batches; i++ ))
do
	line=$(( $i * 1000 ))
	if [ $line -gt $total_lines ]
	then
		last_line=$(( $(($i -1)) * 1000 ))
		nlines=$(( $total_lines - $last_line ))

		tail -n $nlines $OUTPUT_DIR/input_signal.list > $OUTPUT_DIR/pcnv_input_lists/$i'_pcnv_input.list'
	else
		head -n $line $OUTPUT_DIR/input_signal.list | tail -1000 > $OUTPUT_DIR/pcnv_input_lists/$i'_pcnv_input.list'
	fi
done

# Run PennCNV
mkdir -p $OUTPUT_DIR/logs
mkdir -p $OUTPUT_DIR/penncnv_output
for (( i = 0; i <= $num_batches; i++ ))
do
	INPUT=$OUTPUT_DIR/pcnv_input_lists/$i'_pcnv_input.list'
	perl $PENNCNV_PATH/detect_cnv.pl \
		-test \
		-hmm $PENNCNV_PATH/lib/hhall.hmm \
		-pfb $OUTPUT_DIR/UKB.pfb \
		--list $INPUT \
		-log $OUTPUT_DIR/logs/$i'_autosome.log' \
		-out $OUTPUT_DIR/penncnv_output/$i'_pcnv_autosome.txt'
	
	perl $PENNCNV_PATH/detect_cnv.pl \
		-test --chrX \
		-hmm $PENNCNV_PATH/lib/hhall.hmm \
		-pfb $OUTPUT_DIR/UKB.pfb \
		--list $INPUT \
		--sexfile $SEXFILE \
		-log $OUTPUT_DIR/logs/$i'_chrX.log' \
		-out $OUTPUT_DIR/penncnv_output/$i'_pcnv_chrX.txt'
done

# Perform QC of PennCNV calls
mkdir -p $OUTPUT_DIR/QC/logs
mkdir -p $OUTPUT_DIR/QC/qc_pass_samples
mkdir -p $OUTPUT_DIR/QC/qc_sum
mkdir -p $OUTPUT_DIR/QC/goodcnv
for (( i = 0; i <= $num_batches; i++ ))
do
	$PENNCNV_PATH/filter_cnv.pl $OUTPUT_DIR/penncnv_output/$i'_pcnv_autosome.txt' \
		--qclogfile $OUTPUT_DIR/QC/logs/autosome_$i.log \
		--qclrrsd 0.35 \
		--qcnumcnv 30 \
		--qcwf 0.03 \
		--qcpassout $OUTPUT_DIR/QC/qc_pass_samples/autosome_$i.qcpass \
		--qcsumout $OUTPUT_DIR/QC/qc_sum/autosome_$i.qcsum \
		--out $OUTPUT_DIR/QC/goodcnv/autosome_$i.goodcnv
	
	$PENNCNV_PATH/filter_cnv.pl $OUTPUT_DIR/penncnv_output/$i'_pcnv_chrX.txt' \
		--qclogfile $OUTPUT_DIR/QC/logs/chrX_$i.log \
		--qclrrsd 0.35 \
		--qcnumcnv 30 \
		--qcwf 0.03 \
		--chrX \
		--qcpassout $OUTPUT_DIR/QC/qc_pass_samples/chrX_$i.qcpass \
		--qcsumout $OUTPUT_DIR/QC/qc_sum/chrX_$i.qcsum \
		--out $OUTPUT_DIR/QC/goodcnv/chrX_$i.goodcnv
done

