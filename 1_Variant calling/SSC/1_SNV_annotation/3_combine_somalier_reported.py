import pandas as pd
import numpy as np

# Overlap predicted and self-reported ancestry data

# Input and output files
rep_path="/path/to/SSC/core/descriptive/files" # SSC phenotypic data can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
som_output="/path/to/somalier/ancestry/output.tsv" # Use the output of script 2_somalier_ancestry.sh here
output_file="/path/to/output/ancestry/table.csv"

# Parse reported data
groups=['proband', 'mother', 'father', 'sibling', 'sibling2', 'twin']
rep=pd.DataFrame()
for g in groups:
	df=pd.read_csv(f'{rep_path}/{g}_ssc_core_descriptive.csv')
	df=df[['individual', 'race', 'ethnicity']]
	rep=pd.concat([rep, df])

race=['white', 'asian', 'african-amer', 'native-american', 'native-hawaiian']
anc=['EUR', 'EAS', 'AFR', 'AMR', 'OCN']
for i in range(5):
	rep['rep_'+anc[i]]=0
	rep.loc[rep.race==race[i], 'rep_'+anc[i]]=1
rep.loc[rep.ethnicity=='hispanic', 'rep_AMR']=1

# Load somalier-predicted ancestry
pred=pd.read_csv(som_output, sep='\t')
pred['Sample']=pred['#sample_id']

# Merge
df=pd.merge(rep, pred, on='Sample', how='right')
df=df[['Sample', 'predicted_ancestry', 'race', 'ethnicity', 'rep_EAS', 'rep_AFR', 'rep_AMR', 'rep_OCN', 'rep_EUR']]
df.fillna('.', inplace=True)

groups=['EUR', 'ASJ', 'EAS', 'AMR', 'AFR']
for g in groups:
        df[g]=0
        df.loc[df.predicted_ancestry.str.contains(g), g]=1
        if 'rep_'+g in df.columns.to_list():
                df.loc[df['rep_'+g]==1, g]=1

# Reformat
df=df[['Sample']+groups+['predicted_ancestry', 'race', 'ethnicity']]

# Save
df.to_csv(output_file, index=False)

