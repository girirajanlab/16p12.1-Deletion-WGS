import pandas as pd

# Overlap predicted and self-reported ancestry data

# Input and output files
rep_file="/path/to/SPARK/phenotypic/data/individuals_registration.csv" # SPARK phenotypic data can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
som_path="/path/to/somalier/ancestry/directory" # Use the directory script 2_somalier_ancestry.sh output files to
output_file="/path/to/output/ancestry/table.csv"

# Parse reported data
rep=pd.read_csv(rep_file)
rep['Sample']=rep.subject_sp_id

race_cols=[i for i in rep.columns.to_list() if 'race' in i]
rep=rep[['Sample']+race_cols+['hispanic']]
race_cols=['race_asian', 'race_african_amer', 'race_native_amer', 'race_native_hawaiian', 'race_white', 'hispanic']
race_short=['EAS', 'AFR', 'AMR', 'OCN', 'EUR', 'AMR']
for i in range(6):
        if i!=5:
                rep['rep_'+race_short[i]]=0
        rep.loc[rep[race_cols[i]]==1, 'rep_'+race_short[i]]=1

# Load somalier-predicted ancestry
somalier=pd.DataFrame()
for i in range(10):
        subsom=pd.read_csv(f'{som_path}/somalier_ancestry_{i}.somalier-ancestry.tsv', sep='\t')
        subsom=subsom[subsom['#sample_id'].str.contains('SP')]
        subsom['Sample']='SP'+subsom['#sample_id'].str.split('SP', expand=True)[1]
        somalier=pd.concat([somalier, subsom])

# Merge
df=pd.merge(rep, somalier, on='Sample', how='right')
df=df[['Sample', 'predicted_ancestry', 'rep_EAS', 'rep_AFR', 'rep_AMR', 'rep_OCN', 'rep_EUR']+race_cols]
df.fillna('.', inplace=True)

groups=['EUR', 'ASJ', 'EAS', 'AMR', 'AFR']
for g in groups:
        df[g]=0
        df.loc[df.predicted_ancestry.str.contains(g), g]=1
        if 'rep_'+g in df.columns.to_list():
                df.loc[df['rep_'+g]==1, g]=1

# Reformat
df['last_digit']=df.Sample.str[-1]
df=df[['Sample', 'last_digit']+groups+['predicted_ancestry']+race_cols]

# Save
df.to_csv(output_file, index=False)
