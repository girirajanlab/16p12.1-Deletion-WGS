This README describes scripts for compiling and integrating genetic and phenotypic information for analyses

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- DD cohort/
	- _1_make_table_s1a.py_: Compile and organize genetic and phenotypic data for individuals in the DD cohort
	- _2_make_genelists.py_: Create lists of genes with secondary variants in all probands and probands with specific phenotypic domains
	- _3_genelist_by_proband.py_: Create lists of genes with secondary variants for each proband
- MyCode/
	- _1_parse_ICD10.py_: Parse ICD10 medical record data
	- _2_condense_data.py_: Compile and organize genetic and phenotypic data for individuals in MyCode
- Searchlight/
	- _1_compile_data.py_: Compile and organize genetic and phenotypic data for individuals in Searchlight
	- _2_make_genelists.py_: Create lists of secondary variant genes in Searchlight probands
	- _3_gather_SRS.py_: Compile SRS data on all Searchlight 16p11.2 CNV probands
- SPARK/
	- _1_condense_data.py_: Compile and organize genetic and phenotypic data for individuals in SPARK
	- _2_get_controls.py_: Identify age and sex matched control samples in SPARK cohort
- SSC/
	- _1_identify_first_hits.py_: Identify primary CNVs and SNVs for probands in SSC
	- _2_compile_data.py_: Compile and organize genetic and phenotypic data by primary variant for probands in SSC
	- _3_merge_data.py_: Compile data from all primary variants in SSC into a single file
	- _4_make_genelists.py_: Create lists of genes with secondary variants for probands with each primary variant
- UKB/
	- _1_parse_burden.ipynb_: Parse SNV and CNV burden data from UK Biobank samples
	- _2_gather_sample_data.ipynb_: Gather phenotypic information from UK Biobank samples
	- _3_identify_samples.ipynb_: Identify age and sex matched control samples in UK Biobank
	- _4_parse_phenotypes.ipynb_: Parse questionnaire phenotype data for 16p12.1 deletion samples
	- _5_parse_ICD10.ipynb_: Parse ICD10 data for 16p12.1 deletion samples and controls
	- _6_compile_data.ipynb_: Compile and organize genetic and phenotypic data for 16p12.1 deletion carriers in the UK Biobank
	- _7_gather_phewas.ipynb_: Create input files for PheWAS
