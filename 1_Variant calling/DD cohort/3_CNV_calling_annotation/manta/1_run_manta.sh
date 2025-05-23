#!/bin/bash

# Input and output files
INPUT_BAM=/path/to/input.bam # Use the sorted BAM files generated from bwa-mem in script 1_GATK_pipeline.sh
OUTPUT_FILE=/path/to/output.vcf

# Reference files
REF=/path/to/hg19/reference.fasta
CHROMS=/path/to/chromosome_bedfile.bed.gz

# Paths to tools
MANTADIR=/path/to/manta_v1.6.0

# Temporary files
TMPDIR=/path/to/tmp/directory

# Run Manta
mkdir -p $TMPDIR

$MANTADIR/bin/configManta.py \
	--bam=$INPUT_BAM \
	--referenceFasta=$REF \
	--runDir=$TMPDIR \
	--callRegions=$CHROMS

# There is a bug with the line isEmail = isLocalSmtp()
# Replace with isEmail = False
cat $TMPDIR/runWorkflow.py | sed 's/isEmail = isLocalSmtp()/isEmail = False/g' > $TMPDIR/runWorkflow.noemail.py
chmod +x $TMPDIR/runWorkflow.noemail.py
$TMPDIR/runWorkflow.noemail.py -j 40 -g 40

# Copy the diplodSV VCF as the final output
gunzip --keep $TMPDIR/results/variants/diploidSV.vcf.gz
mv $TMPDIR/results/variants/diploidSV.vcf $OUTPUT_FILE

