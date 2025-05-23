import pandas as pd
import numpy as np

# Create lists of the genes affected by secondary variants in the SSC sub-cohorts

# Input and output files
COHORTS="/path/to/cohort/file.csv" # Use the cohort output from script 1_identify_first_hits.py
SNVS="/path/to/SNV/variants.csv" # Use the SNV output of script 1_identify_first_hits.py
CNVS="/path/to/CNV/variants.csv" # Use the CNV output of script 1_identify_first_hits.py
STRS="/path/to/STR/variants.csv" # Use the output of script 1_Variant calling/SSC/3_STR_annotation/10_annotate_loeuf.py

OUTPUT_DIR="/path/to/output/directory" # Output will be one file per cohort containing a list of genes affected by secondary variants

# Find secondary variant genes
cohort=pd.read_csv(COHORTS)
cohorts=['DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'No primary variant']
for co in cohorts:
	df=cohort[cohort[co]][['Sample', co]].copy()
	# SNV
	snvs=pd.read_csv(SNVS)
	if co=='DBD Tier 1 SNVs':
		snvs=snvs[~snvs.primary]
	snvs=snvs[snvs.Sample.isin(df.Sample.to_list())]
	snvs=snvs[['Gene_id_']]
	
	# CNV
	cnvs=pd.read_csv(CNVS)
	if co=='Large rare deletions':
		cnvs=cnvs[(~cnvs.primary) | (cnvs['Del/Dup']!='Del')]
	if co=='Large rare duplications':
		cnvs=cnvs[(~cnvs.primary) | (cnvs['Del/Dup']!='Dup')]
	cnvs=cnvs[cnvs.patientID.isin(df.Sample.to_list())]
	cnvs=cnvs[cnvs['Del/Dup'].isin(['Del', 'Dup'])]
	cnvs=cnvs[['Gene_id_']]
	
	# STRs
	str=pd.read_csv(STRS)
	str=str[str.Sample.isin(df.Sample.to_list())]
	str=str[['Gene_id_']]
	
	all_genes=sorted(list(set(snvs.Gene_id_.to_list()+cnvs.Gene_id_.to_list()+str.Gene_id_.to_list())))
	all_genes=[i for i in all_genes if "ENSGR" not in i]

	print(all_genes[0:10])

	with open(f'{OUTPUT_DIR}/{co}_genes.csv', 'w') as outfile:
		outfile.write('\n'.join(all_genes))
		outfile.write('\n')
