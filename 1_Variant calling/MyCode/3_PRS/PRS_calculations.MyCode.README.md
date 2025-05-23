This README describes the scripts used for PRS calculation in the MyCode cohort

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- _1_signal2lgen.py_: Convert signal files to LGEN format
- _2_prepare_for_imputation.sh_: Prepare array data for imputation
- _3_vcf2plink.sh_: Convert imputed VCFs into PLINK format
- _4_peddy_ancestry.sh_: Calculate ancestry using peddy
- _5_combine_peddy_reported.py_: Identify European samples from self reported and peddy-calculated ancestry
- _6_run_ldpred2.sh_: Calculate PRS using LDpred2
- _7_merge_scores.py_: Merge all PRS into a single file
- helper_scripts/
	- _heterozygosity_filter.R_: Filter samples with high heterozygosity
	- _ldpred2.R_: Run LDpred2 to calculate PRS
		