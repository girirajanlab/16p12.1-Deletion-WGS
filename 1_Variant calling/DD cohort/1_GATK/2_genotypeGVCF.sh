#!/bin/bash

# Apply variant recalibration to GATK processed variants

# Tools for merging GVCFs
TRIMMOMATIC=/path/to/trimmomatic_v.036
BWADIR=/path/to/bwa_v0.7.13
SAMDIR=/path/to/samtools_v1.8
GATKDIR=/path/to/GATK_v3.7
GATK4DIR=/path/to/GATK_v4.0.11.0 # GATK v4 is only used for IndexFeatureFile
PICARD=/path/to/picard_tools_v2.9.0
VCFTOOLS=/path/to/vcftools

# Reference file locations
REF=/path/to/hg19_reference
# GATK resource data for hg19 reference can be accessed from a Google Cloud bucket: gs://gatk-legacy-bundles
DBSNP=/path/to/dbsnp_vcf
OMNI=/path/to/1000G_omni_vcf
HAPMAP=/path/to/hapmap_vcf
MILLS=/path/to/mills_and_1000G_gold_standard_indels_vcf
G1K=/path/to/1000G_phase1_indesl_vcf

# Input files
INPUT_VCFS=/path/to/gatk_output # path to the output files from script 1_GATK_pipeline.sh

# Combine GVCF files from all samples into a single file
perl $VCFTOOLS/src/perl/vcf-concat $INPUT_VCFS/*.vcf | perl $VCFTOOLS/src/perl/vcf-sort -t ./tmp | bgzip -c > genotype_gvcf.vcf.gz
$GATK4DIR/gatk --java-options "-Xmx12g -Xms12g" IndexFeatureFile -F genotype_gvcf.vcf.gz

# Perform and apply variant recalibration in multiple stages
java -d64 -Xmx1200g -jar $GATKDIR/GenomeAnalysisTK.jar \
 -T VariantRecalibrator \
 -R $REF \
 -input genotype_gvcf.vcf.gz \
 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
 -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI \
 -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1K \
 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
 -an DP \
 -an QD \
 -an FS \
 -an SOR \
 -an MQ \
 -an MQRankSum \
 -an ReadPosRankSum \
 -an InbreedingCoeff \
 --mode SNP \
 --recal_file recalibrate_SNP.recal \
 --tranches_file recalibrate_SNP.tranches \
 --rscript_file recalibrate_SNP_plots.R \
 -nt 100

java -d64 -Xmx1200g -Djava.io.tmpdir=./tmp -jar $GATKDIR/GenomeAnalysisTK.jar \
 -T ApplyRecalibration \
 -R $REF \
 -input genotype_gvcf.vcf.gz \
 -mode SNP \
 --ts_filter_level 99.5 \
 -recalFile recalibrate_SNP.recal \
 -tranchesFile recalibrate_SNP.tranches \
 -o genotype_raw_recal.vcf \
 -nt 100

java -d64 -Xmx1200g -jar $GATKDIR/GenomeAnalysisTK.jar \
 -T VariantRecalibrator \
 -R $REF \
 -input genotype_raw_recal.vcf \
 --resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS \
 --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
 -an QD \
 -an DP \
 -an FS \
 -an SOR \
 -an MQRankSum \
 -an ReadPosRankSum \
 -an InbreedingCoeff \
 -mode INDEL \
 --maxGaussians 4 \
 -recalFile recalibrate_INDEL.recal \
 -tranchesFile recalibrate_INDEL.tranches \
 -rscriptFile recalibrate_INDEL_plots.R \
 -nt 100

java -d64 -Xmx1200g -Djava.io.tmpdir=./tmp -jar $GATKDIR/GenomeAnalysisTK.jar \
 -T ApplyRecalibration \
 -R $REF \
 -input genotype_raw_recal.vcf \
 -mode INDEL \
 --ts_filter_level 99.0 \
 -recalFile recalibrate_INDEL.recal \
 -tranchesFile recalibrate_INDEL.tranches \
 -o genotype_recal.vcf \
 -nt 100

java -d64 -Xmx1200g -jar $GATKDIR/GenomeAnalysisTK.jar \
 -T VariantFiltration \
 -R $REF \
 -o genotype_final.vcf \
 --variant genotype_recal.vcf \
 --clusterWindowSize 10 \
 --clusterSize 3 \
 --filterExpression "DP < 5 " --filterName "LowCoverage" \
 --filterExpression "QUAL < 30.0 " --filterName "VeryLowQual" \
 --filterExpression "QUAL > 30.0 && QUAL < 50.0 " --filterName "LowQual" \
 --filterExpression "QD < 1.5 " --filterName "LowQD"
