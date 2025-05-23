import pandas as pd

# Self reported ancestry
selfreport='/path/to/self_reported_ancestry.csv'
output_file='/path/to/output/ancestry.csv'

# Match the somalier output ancestry annotations with self reported ancestry annotations
somalier=pd.read_csv('somalier_output/somalier_ancestry.somalier-ancestry.tsv', sep='\t')
somalier=somalier[somalier['#sample_id'].str.contains('SG')]
somalier['IID']='SG'+somalier['#sample_id'].str.split('SG', expand=True)[1]

known=pd.read_csv(self_report)

# Convert the known ancestry information into broader classes
anc_map={'Caucasian':'EUR', 'British':'EUR', 'White':'EUR', 'European':'EUR', 'English':'EUR', 'Australian/Italian':'EUR',
			'Aboriginal Australian':'OCE',
			'Indian':'SAS',
			'White/Jewish':'EUR/ASJ',
			'African':'AFR',
			'Chinese':'EAS', 'Asian':'EAS',
			'Chinese/European':'EUR/EAS',
			'Caucasian/American Indian':'EUR/AMR',
			'Caucasian/Lebanse':'EUR/NEA',
			'Lebanese':'NEA'}
known['self_short']=known.self_reported_ancestry.map(anc_map)

# Use union of self report and somalier to determine which ancestry filters need to be applied
df=pd.merge(known, somalier, on='IID', how='inner')
df=df[['IID', 'self_reported_ancestry', 'self_short', 'predicted_ancestry']]
df.fillna('.', inplace=True)

groups=['EUR', 'ASJ', 'EAS', 'AMR', 'AFR']
for g in groups:
	df[g]=0
	df.loc[df.self_short.str.contains(g), g]=1
	df.loc[df.predicted_ancestry.str.contains(g), g]=1

# Save
df.to_csv(output_file, index=False)