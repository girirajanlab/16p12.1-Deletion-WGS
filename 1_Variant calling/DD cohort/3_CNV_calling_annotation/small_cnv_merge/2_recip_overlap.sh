#!/bin/bash

# Check calls from each caller to identify calls with 50% reciprocal overlap between callers

# Input and output files
SAMPLE=sample_id
INPUT_PATH=/path/to/input/files # Use the files generated from script 1_separate_cnvs.sh
OUTPUT_PATH=/path/to/output/files

CNVNATOR=$INPUT_PATH/$SAMPLE.cnvnator.bed
MANTA=$INPUT_PATH/$SAMPLE.manta.bed
DELLY=$INPUT_PATH/$SAMPLE.delly.bed
LUMPY=$INPUT_PATH/$SAMPLE.lumpy.bed

# Paths to tools
BEDPATH=/path/to/bedtools_v2.27.1

# Use bedtools to annotate the number of calls with 50% reciprocal overlaps from every caller
# CNVnator-Manta
$BEDPATH/bedtools intersect -a $CNVNATOR -b $MANTA -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.cnvnator.manta.bed
# Add Delly
$BEDPATH/bedtools intersect -a $OUTDIR/$SAMPLE.cnvnator.manta.bed -b $DELLY -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.cnvnator.manta.delly.bed
rm $OUTDIR/$SAMPLE.cnvnator.manta.bed
# Add Lumpy
$BEDPATH/bedtools intersect -a $OUTDIR/$SAMPLE.cnvnator.manta.delly.bed -b $LUMPY -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.cnvnator.manta.delly.lumpy.bed
rm $OUTDIR/$SAMPLE.cnvnator.manta.delly.bed

# Manta-CNVnator
$BEDPATH/bedtools intersect -a $MANTA -b $CNVNATOR -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.manta.cnvnator.bed
# Add Delly
$BEDPATH/bedtools intersect -a $OUTDIR/$SAMPLE.manta.cnvnator.bed -b $DELLY -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.manta.cnvnator.delly.bed
rm $OUTDIR/$SAMPLE.manta.cnvnator.bed
# Add Lumpy
$BEDPATH/bedtools intersect -a $OUTDIR/$SAMPLE.manta.cnvnator.delly.bed -b $LUMPY -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.manta.cnvnator.delly.lumpy.bed
rm $OUTDIR/$SAMPLE.manta.cnvnator.delly.bed

# Delly-CNVnator
$BEDPATH/bedtools intersect -a $DELLY -b $CNVNATOR -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.delly.cnvnator.bed
# Add Manta
$BEDPATH/bedtools intersect -a $OUTDIR/$SAMPLE.delly.cnvnator.bed -b $MANTA -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.delly.cnvnator.manta.bed
rm $OUTDIR/$SAMPLE.delly.cnvnator.bed
# Add Lumpy
$BEDPATH/bedtools intersect -a $OUTDIR/$SAMPLE.delly.cnvnator.manta.bed -b $LUMPY -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.delly.cnvnator.manta.lumpy.bed
rm $OUTDIR/$SAMPLE.delly.cnvnator.manta.bed

# Lumpy-CNVnator
$BEDPATH/bedtools intersect -a $LUMPY -b $CNVNATOR -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.lumpy.cnvnator.bed
# Add Manta
$BEDPATH/bedtools intersect -a $OUTDIR/$SAMPLE.lumpy.cnvnator.bed -b $MANTA -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.lumpy.cnvnator.manta.bed
rm $OUTDIR/$SAMPLE.lumpy.cnvnator.bed
# Add Delly
$BEDPATH/bedtools intersect -a $OUTDIR/$SAMPLE.lumpy.cnvnator.manta.bed -b $DELLY -wa -f 0.5 -r -c > $OUTDIR/$SAMPLE.delly.cnvnator.manta.lumpy.bed
rm $OUTDIR/$SAMPLE.lumpy.cnvnator.manta.bed

# Remove cany calls that are only present in one caller (counts for other 3 callers are 0)
grep -vP '\t0\t0\t0' $OUTDIR/$SAMPLE.cnvnator.manta.delly.lumpy.bed > $OUTDIR/$SAMPLE.combined.bed
grep -vP '\t0\t0\t0' $OUTDIR/$SAMPLE.manta.cnvnator.delly.lumpy.bed >> $OUTDIR/$SAMPLE.combined.bed
grep -vP '\t0\t0\t0' $OUTDIR/$SAMPLE.delly.cnvnator.manta.lumpy.bed >> $OUTDIR/$SAMPLE.combined.bed
grep -vP '\t0\t0\t0' $OUTDIR/$SAMPLE.delly.cnvnator.manta.lumpy.bed >> $OUTDIR/$SAMPLE.combined.bed
