#!/bin/bash

# Run PennCNV to call CNVs from microarray

# Input and output files
SIGNAL_PATH=/path/to/SPARK/iWGS_v1.1/genotyping/signal # Path to array signal files for SPARK samples
ARRAY_RESOURCE=/path/to/SPARK/iWGS_v1.1/resources/InfiniumCoreExome-24v1-1_A2.GRCh38.csv # Array information table
# The input files needed can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
SEXFILE=/path/to/sample/sexes.list # A list of sample sexes for chrX CNV calling
OUTPUT_DIR=/path/to/output/directory

# Tools
PENNCNV_PATH=/path/to/PennCNV-1.0.5

# Create a list of signal files
find $SIGNAL_PATH/ -name "*.signal.gz" > $OUTPUT_DIR/SPARK_signal.list

# Create a PFB file
# PennCNV cannot use gzipped signal files for this, so unzip the signal files for 1000 random samlpes to create a PFB
mkdir -p $OUTPUT_DIR/pfb_sample_files
shuf -f 1000 $OUTPUT_DIR/SPARK_signal.list > $OUTPUT_DIR/tmp_pfb_input_files.list
for i in {1..1000}
do
	FILE=`head -n $i $OUTPUT_DIR/tmp_pfb_input_files.list | tail -1`
	SAMPLE=`head -n $i head -n $i $OUTPUT_DIR/tmp_pfb_input_files.list | tail -1 | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.'`

	echo -e "Name\t"$SAMPLE".B Allele Freq" > $OUTPUT_DIR/pfb_sample_files/$SAMPLE'_pfb_input.txt'
	gunzip -c $FILE | tail -n +2 | awk '{ print $1 "\t" $4 }' >> $OUTPUT_DIR/pfb_sample_files/$SAMPLE'_pfb_input.txt'
done
rm $OUTPUT_DIR/tmp_pfb_input_files.list
find $OUTPUT_DIR/pfb_sample_files/ -name -name "*_pfb_input.txt" > $OUTPUT_DIR/pfb_input_files.list

# Create a SNP position file
echo -e 'Name\tChr\tPos' > $OUTPUT_DIR/SPARK_SNP_Positions.list
cut -f 2,10,11 -d , $ARRAY_RESOURCE | tail -n+9 | sed 's/,/\t/g' >> $OUTPUT_DIR/SPARK_SNP_Positions.list

# Compile PFB
$PENNCNV_PATH/compile_pfb.pl --listfile $OUTPUT_DIR/pfb_input_files.list --snpposfile $OUTPUT_DIR/SPARK_SNP_Positions.list --output $OUTPUT_DIR/SPARK.pfb

# Make a GC model adjustment file
mkdir $OUTPUT_DIR/GC_files
gunzip -c $PENNCNV_PATH/gc_file/hg38.gc5Base.txt.gz > $OUTPUT_DIR/GC_files/hg38.gc5Base.txt
$PENNCNV_PATH/cal_gc_snp.pl $OUTPUT_DIR/GC_files/hg38.gc5Base.txt $OUTPUT_DIR/SPARK.pfb > $OUTPUT_DIR/GC_files/SPARK.hg38.gcmodel

# Create lists of 100 files to run concurrently in PennCNV
mkdir -p $OUTPUT_DIR/pcnv_input_lists
total_lines=`wc -l $OUTPUT_DIR/SPARK_signal.list | cut -f 1 -d ' '`
num_batches=$(( $(($total_lines/100))+1 ))
for (( i = 0; i <= $num_batches; i++ ))
do
	line=$(( $i * 100 ))
	if [ $line -gt $total_lines ]
	then
		last_line=$(( $(($i -1)) * 100 ))
		nlines=$(( $total_lines - $last_line ))

		tail -n $nlines $OUTPUT_DIR/SPARK_signal.list > tmp1
		yes '`echo "Name\tChr\tPosition\t.B Allele Freq\t.Log R Ratio" ; (gunzip -c' | head -n $nlines > tmp1
		yes '| tail -n +2)`' | head -n $nlines > tmp3
		paste -d " " tmp1 tmp2 tmp3 > $OUTPUT_DIR/pcnv_input_lists/$i'_pcnv_input.list'
	else
		head -n $line $OUTPUT_DIR/SPARK_signal.list | tail -100 > tmp2
		yes '`echo "Name\tChr\tPosition\t.B Allele Freq\t.Log R Ratio" ; (gunzip -c' | head -n 100 > tmp1
		yes '| tail -n +2)`' | head -n 100 > tmp3
		paste -d " " tmp1 tmp2 tmp3 > $OUTPUT_DIR/pcnv_input_lists/$i'_pcnv_input.list'
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
		-pfb $OUTPUT_DIR/SPARK.pfb \
		--list $INPUT \
		-log $OUTPUT_DIR/logs/$i'_autosome.log' \
		-out $OUTPUT_DIR/penncnv_output/$i'_pcnv_autosome.txt' \
		--gcmodel $OUTPUT_DIR/GC_files/SPARK.hg38.gcmodel
	
	perl $PENNCNV_PATH/detect_cnv.pl \
		-test --chrX \
		-hmm $PENNCNV_PATH/lib/hhall.hmm \
		-pfb $OUTPUT_DIR/SPARK.pfb \
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
		--qcpassout $OUTPUT_DIR/QC/qc_pass_samples/autosome_$i.qcpass \
		--qcsumout $OUTPUT_DIR/QC/qc_sum/autosome_$i.qcsum \
		--out $OUTPUT_DIR/QC/goodcnv/autosome_$i.goodcnv
	
	$PENNCNV_PATH/filter_cnv.pl $OUTPUT_DIR/penncnv_output/$i'_pcnv_chrX.txt' \
		--qclogfile $OUTPUT_DIR/QC/logs/chrX_$i.log \
		--qclrrsd 0.35 \
		--chrX \
		--qcpassout $OUTPUT_DIR/QC/qc_pass_samples/chrX_$i.qcpass \
		--qcsumout $OUTPUT_DIR/QC/qc_sum/chrX_$i.qcsum \
		--out $OUTPUT_DIR/QC/goodcnv/chrX_$i.goodcnv
done

