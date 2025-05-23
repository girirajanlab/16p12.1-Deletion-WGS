#!/bin/bash

# Tools for FASTQ processing
TRIMMOMATIC=/path/to/trimmomatic_v.036
BWADIR=/path/to/bwa_v0.7.13
SAMDIR=/path/to/samtools_v1.8
GATKDIR=/path/to/GATK_v3.7
PICARD=/path/to/picard_tools_v2.9.0

# Reference files
REF=/path/to/hg19_reference
# GATK resource data for hg19 reference can be accessed from a Google Cloud bucket: gs://gatk-legacy-bundles
DBSNP=/path/to/dbsnp_vcf
OMNI=/path/to/1000G_omni_vcf
HAPMAP=/path/to/hapmap_vcf
MILLS=/path/to/mills_and_1000G_gold_standard_indels_vcf
G1K=/path/to/1000G_phase1_indesl_vcf

# Read locations
FWD=/path/to/R1.fastq.gz
REV=/path/to/R2.fastq.gz

java -jar $TRIMMOMATIC/trimmomatic-0.36.jar PE -threads 10 -phred33 $FWD $REV R1_trim.fq.gz R1_unpaired.fq.gz R2_trim.fq.gz R2_unpaired.fq.gz LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20

$BWADIR/bwa mem -t 10 -R "@RG\tID:\tSM:\tPL:illumina\tLB:lib1" $REF R1_trim.fq.gz R2_trim.fq.gz > out.sam
$SAMDIR/samtools view -@ 10 -h -m 15G -b -T $REF -o out.bam out.sam
$SAMDIR/samtools sort -@ 10 -m 14G -o sorted.bam out.bam
$SAMDIR/samtools index sorted.bam
rm out.sam

java -d64 -Xmx150g -Djava.io.tmpdir=tmp -jar $PICARD/picard.jar MarkDuplicates \
 I=sorted.bam \
 O=rmdup.bam \
 REMOVE_DUPLICATES=true \
 M=metrics.txt

$SAMDIR/samtools index rmdup.bam

java -d64 -Xmx150g -jar $GATKDIR/GenomeAnalysisTK.jar \
 -T BaseRecalibrator \
 -R $REF \
 -nct 10 \
 -I rmdup.bam \
 -knownSites $DBSNP \
 -knownSites $MILLS \
 -knownSites $G1K \
 -o BQSR_report.table

java -d64 -Xmx150g -jar $GATKDIR/GenomeAnalysisTK.jar \
 -T PrintReads \
 -R $REF \
 -nct 10 \
 -I rmdup.bam \
 -BQSR BQSR_report.table \
 -o rmdup_BQSR.bam

java -d64 -Xmx150g -jar $GATKDIR/GenomeAnalysisTK.jar \
 -T HaplotypeCaller \
 -R $REF \
 --dbsnp $DBSNP \
 -I rmdup_BQSR.bam \
 -o gvcf.g.vcf \
 --genotyping_mode DISCOVERY \
 -mbq 20 \
 --emitRefConfidence GVCF

