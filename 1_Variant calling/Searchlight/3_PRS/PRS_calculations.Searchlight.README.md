This README describes the scripts used for PRS calculations in the Simons Searchlight cohort

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- _1_prepare_for_imputation.sh_: Prepare array data for imputation
- _2_vcf2plink.sh_: Convert imputed VCFs into PLINK format
- _3_peddy_ancestry.sh_: Calculate ancestry using peddy
- _4_combine_peddy_reported.py_: Identify European samples from reported and peddy-calculated ancestry
- _5_run_ldpred2.sh_: Calculate PRS using LDpred2
- _6_merge_scores.py_: Merge all PRS into a single file
- helper_scripts/
	- _heterozygosity_filter.R_: Filter samples with high heterozygosity
	- _ldpred2.R_: Run LDpred2 to calculate PRS
