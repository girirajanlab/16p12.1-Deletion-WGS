import pandas as pd
import numpy as np

# Overlap predicted and self-reported ancestry data

# Input and output files
rep_path="/path/to/Searchlight/svip_phase_2_subjects.csv" # Searchlight phenotypic data can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
som_output="/path/to/somalier/ancestry/output.tsv" # Use the output of script 2_somalier_ancestry.sh here
sample_list="/path/to/list/of/Searchlight/samples.txt"
output_file="/path/to/output/ancestry/table.csv"

# Parse reported data
rep=pd.read_csv(rep_path)
rep=rep[['sfari_id', 'race']]

rep['Sample']=rep.sfari_id

rep['rep_EUR']=0
rep.loc[rep.race=='white', 'rep_EUR']=1
rep['rep_OCN']=0
rep.loc[rep.race=='pacific-islander', 'rep_OCN']=1

rep['som_form']=rep.Sample+rep.Sample

# Load somalier-predicted ancestry
somalier=pd.read_csv(som_output, sep='\t')
ssamples=pd.read_csv(sample_list, sep='\t', header=None, names=['Sample'])
ssamples['som_form']=ssamples.Sample+ssamples.Sample
somalier=somalier[somalier['#sample_id'].isin(ssamples.som_form.to_list())]
somalier['som_form']=somalier['#sample_id'].str.replace('.', '-')

# Merge
df=pd.merge(rep, somalier, on='som_form', how='right')
df['Sample']=df.som_form.apply(lambda x: x[:int((len(x)/2))])
df=df[['Sample', 'predicted_ancestry', 'race', 'rep_OCN', 'rep_EUR']]
df.fillna('.', inplace=True)

groups=['EUR', 'ASJ', 'EAS', 'AMR', 'AFR']
for g in groups:
        df[g]=0
        df.loc[df.predicted_ancestry.str.contains(g), g]=1
        if 'rep_'+g in df.columns.to_list():
                df.loc[df['rep_'+g]==1, g]=1

# Reformat
df['Sample']=df.Sample.str.replace('-', '.')
df=df[['Sample']+groups+['predicted_ancestry', 'race']]

# Save
df.to_csv(output_file, index=False)

