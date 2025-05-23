import pandas as pd

# Input and output files
reported_ancestry="/path/to/reported/ancestry.csv" # Reported ancestry data was provided by MyCode
som_output="/path/to/somalier/ancestry/output.tsv" # Use the output of script 5_somalier_ancestry.sh here
output_file="/path/to/output/ancestry/table.csv"

# Parse reported data
rep=pd.read_csv(reported_ancestry)
rep=rep[['SEQN_ID', 'PT_RACE_1', 'PT_ETHNICITY']]

rep['anc']='unknown'
rep.loc[rep.PT_RACE_1=='White', 'anc']='EUR'
rep.loc[rep.PT_RACE_1=='Black Or African American', 'anc']='AFR'
rep.loc[rep.PT_ETHNICITY=='Hispanic or Latino', 'anc']='AMR'

rep['Sample']=rep.SEQN_ID

# Load somalier-predicted ancestry
somalier=pd.read_csv(som_output, sep='\t')
somalier=somalier[somalier['#sample_id'].str.contains('GHS')]
somalier['Sample']='GHS'+somalier['#sample_id'].str.split('GHS', expand=True)[1]

# Merge
df=pd.merge(rep, somalier, on='Sample', how='right')
df=df[['Sample', 'anc', 'predicted_ancestry']]
df.fillna('.', inplace=True)

groups=['EUR', 'ASJ', 'EAS', 'AMR', 'AFR']
for g in groups:
        df[g]=0
        df.loc[df.anc.str.contains(g), g]=1
        df.loc[df.predicted_ancestry.str.contains(g), g]=1

# Save
df.to_csv(output_file, index=False)

