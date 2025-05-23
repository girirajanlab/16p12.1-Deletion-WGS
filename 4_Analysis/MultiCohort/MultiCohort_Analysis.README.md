This README describes the scripts used for analysis of across cohorts

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- MetaAnalysis/
	- _1_compile_data.py_: Compile t-test data from DD cohort, UKB, AoU, and SPARK
	- _2_meta_analysis.R_: Perform meta-analysis of the t-test results
	- _3_plot_meta.py_: Plot results of meta-analysis
- Phenotype Prevalence/
	- _1_phenotype_frequencies.py_: Calculate differences in phenotype frequencies across cohorts of 16p12.1 deletion carriers
- Searchlight SSC/
	- _1_regression_heatmap.py_: Plot regression results from SSC and Searchlight probands together