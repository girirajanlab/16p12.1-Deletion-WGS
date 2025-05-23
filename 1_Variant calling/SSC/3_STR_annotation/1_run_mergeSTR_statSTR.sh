#!/bin/bash

# Merge and calculate statistics on STR calls from Mitra et al. Nature 2021 (PubMed: https://pmc.ncbi.nlm.nih.gov/articles/PMC7810352/)

# Input and output files
STR_PATH=/path/to/str/vcf/directory # STR VCFs are from Mitra et al. Nature 2021 and can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
OUTPUT_DIR=/path/to/output/directory

# Tools
MERGESTR=/path/to/mergeSTR
STATSTR=/path/to/statSTR

# Merge VCFs per chromosome
# VCFs are separated by chromosome with 5 phases per chromosome
for i in `seq 1 23`
do
	chrom=$i
	if [[ $i -eq 23 ]]
	then
		chrom=X
	fi
	vcfs=vcfs="${STR_PATH}/phase1_${chrom}.sorted.vcf.gz,${STR_PATH}/phase2_${chrom}.sorted.vcf.gz,${STR_PATH}/phase3_1_${chrom}.sorted.vcf.gz,${STR_PATH}/phase3_2_${chrom}.sorted.vcf.gz,${STR_PATH}/phase4_${chrom}.sorted.vcf.gz"
	$MERGESTR/mergeSTR --vcfs $vcfs --out $OUTPUT_DIR/merge/chr${chrom}
	bgzip $OUTPUT_DIR/merge/chr${chrom}.vcf
	tabix -p vcf $OUTPUT_DIR/merge/chr${chrom}.vcf.gz
done

# Get statistics by chromosome
for i in `seq 1 23`
do
	chrom=$i
	if [[ $i -eq 23 ]]
	then
		chrom=X
	fi
	$STATSTR/statSTR --vcf $OUTPUT_DIR/merge/chr${chrom}.vcf \
		--thresh \
		--afreq \
		--acount \
		--mean \
		--mode \
		--var \
		--numcalled \
		--out $OUTPUT_DIR/statstr/chr{$chrom}
done

