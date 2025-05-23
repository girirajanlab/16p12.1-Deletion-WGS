import pandas as pd
import numpy as np

# Input and output files
INPUT_PATH="/path/to/cohort/data/tables/" # Use the OUTPUT_PATH from script 2_compile_data.py
OUTPUT_FILE="/path/to/output/file.csv" # This file will contain data for the whole cohort, merged into a single file

# Merge the data from all SSC cohorts together
# If a sample has multiple first hits, remove all the first hits from the burden table
cohorts=['DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'No primary variant']
df=pd.DataFrame()
for co in cohorts:
	codf=pd.read_csv(f'{INPUT_PATH}/{co}.csv')
	codf[co]=True
	cols=[i for i in codf.columns.to_list() if i!='Cohort']
	df=pd.concat([df, codf[cols]])

for co in cohorts:
	df.loc[df[co].isnull(), co]=False

pheno_cols=['Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Social responsiveness (SRS)',
			'Repetitive behavior (RBS-R)', 'Coordination disorder (DCDQ)', 'BMI z-score']
prs_cols=['SCZ PRS', 'Intelligence PRS', 'Education PRS', 'Autism PRS']

var_cols=['All coding SNVs', 'Missense', 'LOF', 'Splice', 'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)',
			'Genes del.', 'Genes dup.', 'Genes del. (LF)', 'Genes dup. (LF)',
			'STRs', 'STRs (LF)']
			
agg_dict={}
for vc in var_cols:
	agg_dict[vc]='min'
for co in cohorts:
	agg_dict[co]='any'

df.fillna(1000, inplace=True)

df.sort_values(by='Sample', inplace=True)
df=df.groupby(['Sample', 'Sex', 'Age']+pheno_cols+prs_cols).agg(agg_dict)
df.reset_index(inplace=True)
df.sort_values(by='Sample', inplace=True)

df=df.replace(1000, np.nan)

df['Any variant']=df[cohorts[0:3]].any(axis=1)

df=df[['Sample', 'DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'Any variant', 'Sex', 'Age',
			'All coding SNVs', 'Missense', 'LOF', 'Splice', 'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)',
			'Genes del.', 'Genes dup.', 'Genes del. (LF)', 'Genes dup. (LF)',
			'STRs', 'STRs (LF)',
			'SCZ PRS', 'Intelligence PRS', 'Education PRS', 'Autism PRS',
			'Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Social responsiveness (SRS)',
			'Repetitive behavior (RBS-R)', 'Coordination disorder (DCDQ)', 'BMI z-score']]

# Save
df.to_csv(OUTPUT_FILE, index=False)
