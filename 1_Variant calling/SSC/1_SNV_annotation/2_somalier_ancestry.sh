#!/bin/bash

# Calculate ancestry using somalier

# Input and output files
OUTPUT_DIR=/path/to/output/directory
SAMPLE_LIST=/path/to/list/of/samples.txt

# Reference files
REF=/path/to/hg38/reference.fasta
# Reference files can be downloaded from the somalier github page
SITES=/path/to/hg38/sites.vcf.gz # https://github.com/brentp/somalier/releases/tag/v0.2.19
ANC_LABELS=/path/to/ancestry/labels.tsv # https://github.com/brentp/somalier/wiki/ancestry
SOMALIER_GENO=/path/to/somalier/genotype/files/*.somalier # https://github.com/brentp/somalier/wiki/ancestry

# Tools
BCFPATH=/path/to/bcftools_v1.12
TABIXPATH=/path/to/tabix
JVARKIT=/path/to/jvarkit

# Subset sites and split samples on chromsome VCFs
for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y;
do
	INPUT_VCF=/path/to/SFARI/SSC/WGS2/Project_SSC_9209samples.JGvariants.2019-06-21/CCDG_9000JG_B01_GRM_WGS_2019-03-21_chr${CHROM}.recalibrated_variants.flagged.vcf.gz
	# This input VCF can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/

	# Subset out relevant sites from population VCF
	$BCFPATH/bcftools view -R $SITES $INPUT_VCF | bgzip > $OUTPUT_DIR/$CHROM.subset.vcf.gz
	$TABIXPATH/tabix -p vcf $OUTPUT_DIR/$CHROM.subset.vcf.gz

	# Split population VCF into individual VCFs for ancestry filtering
	java -Xmx50G -jar $JVARKIT/dist/biostar130456.jar -x -z -p "$OUTPUT_DIR/$CHROM/sample_split/__SAMPLE__.vcf" $OUTPUT_DIR/$CHROM.subset.vcf.gz
done

# Merge single-chromsome VCFs by sample and extract needed sites for somalier
while read $SAMPLE
do
	touch tmp/$SAMPLE/file_list.txt
	for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y;
	do
		echo $OUTPUT_DIR/$CHROM/sample_split/$SAMPLE.vcf >> tmp/$SAMPLE/file_list.txt
	done

	# Merge individual chromosome VCFs
	$BCFPATH/bcftools concat -f tmp/$SAMPLE/file_list.txt | bgzip > $OUTPUT_DIR/$SAMPLE.vcf.gz
	$TABIXPATH/tabix -p vcf $OUTPUT_DIR/$SAMPLE.vcf.gz

	# Extract sites for somalier
    somalier extract -d extracted/ -s $SITES -f $REF --sample-prefix=$SAMPLE $OUTPUT_DIR/$SAMPLE.vcf.gz

	# Clean up
	rm tmp/$SAMPLE/file_list.txt

done < $SAMPLE_LIST

# Calculate ancestry
somalier ancestry -o somalier_ancestry --labels $ANC_LABELS $SOMALIER_GENO ++ extracted/*.somalier

