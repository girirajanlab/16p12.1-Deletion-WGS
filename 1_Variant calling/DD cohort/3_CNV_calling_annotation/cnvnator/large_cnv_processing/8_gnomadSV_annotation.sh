#!/bin/bash

# Annotate the gnomAD SV frequency

# Input and output files
INPUT_PATH=/path/to/cnvs/files/ # Use the path to the files from script 7_frequency_filter.py
OUTPUT_PATH=/path/to/output/files/

# Reference files
GNOMADSV_PATH=/path/to/gnomadsv/sites/gnomad_v2.1_sv.sites/files
# Folder contains gnomAD SV sites split by deletions and duplications
# GnomADSV sites can be downloaded from gnomAD

# Paths to tools
BEDPATH=/path/to/bedtools_v2.27.1

# Create a header line
gnomad_head=`head -1 $GNOMADSV_PATH/gnomad_v2.1_sv.sites.bed`
echo  -e 'Chr\tStart\tEnd\tType\tName\tLength\tSample\tMicroarray_count\tIntracohort_count\tmicroarray_freq\t'$gnomad_head > header.txt
sed -i 's/ \+/\t/g' header.txt

cat header.txt > $OUTPUT_PATH/cnv_gnomadSV_anno_dels.bed
cat header.txt > $OUTPUT_PATH/pathogenic_cnv_gnomadSV_anno_dels.bed
cat header.txt > $OUTPUT_PATH/cnv_gnomadSV_anno_dups.bed
cat header.txt > $OUTPUT_PATH/pathogenic_cnv_gnomadSV_anno_dups.bed

rm header.txt

# Because of the way the gnomAD SVs are labelled, this will be a more complicated annotation than the one used for microarray controls and intracohort
# Build a lookup table of CNVs and all gnomADSV CNVs it is associated with
$BEDPATH/bedtools intersect -a $INPUT_PATH/cnv_frequency_filter_dels.bed -b $GNOMADSV_PATH/gnomad_v2.1_sv.dels.bed -loj -f 0.5 -r >> $OUTPUT_PATH/cnv_gnomadSV_anno_dels.bed
$BEDPATH/bedtools intersect -a $INPUT_PATH/pathogenic_cnv_frequency_filter_dels.bed -b $GNOMADSV_PATH/gnomad_v2.1_sv.dels.bed -loj -f 0.5 -r >> $OUTPUT_PATH/pathogenic_cnv_gnomadSV_anno_dels.bed

$BEDPATH/bedtools intersect -a $INPUT_PATH/cnv_frequency_filter_dups.bed -b $GNOMADSV_PATH/gnomad_v2.1_sv.dups.bed -loj -f 0.5 -r >> $OUTPUT_PATH/cnv_gnomadSV_anno_dups.bed
$BEDPATH/bedtools intersect -a $INPUT_PATH/pathogenic_cnv_frequency_filter_dups.bed -b $GNOMADSV_PATH/gnomad_v2.1_sv.dups.bed -loj -f 0.5 -r >> $OUTPUT_PATH/pathogenic_cnv_gnomadSV_anno_dups.bed

