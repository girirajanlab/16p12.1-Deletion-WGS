import pandas as pd

# Calculate the frequencies of each phenotype in the DD cohort adults and children

# Input and output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
OUTPUT_TAB="/path/to/output/prevalence/data.csv"

# Load data
df=pd.read_csv(TABS1A)
df=df[(df['16p12.1 deletion']=='Carrier')]

df['Cohort']='DD (adults)'
df.loc[df.Relationship.isin(['Proband', 'Sibling', 'Cousin']), 'Cohort']='DD (children)'

df=df[['Cohort',
		'Depression (Questionnaire)', 'Anxiety (Questionnaire)', 'Sleep trouble (Questionnaire)', 'Psychosis (Questionnaire)', 'Addiction (Questionnaire)', 'Mood lability (Questionnaire)',
		'Birth/pregnancy complications', 'Preterm birth', 'Microcephaly', 'Macrocephaly', 'Strabismus',
		'Seizures', 'Heart defects', 'Hearing loss', 'Vision problems', 'Feeding problems',
		'Obesity', 'ID/DD', 'Motor delay', 'Speech delay', 'Language disorder', 'Aide in school',
		'Learning disability', 'ASD', 'ADHD', 'OCD', 'Schizophrenia', 'BPD', 'Depression', 'Anxiety', 'Sleep trouble', 'PDD']]

# Convert to counts
posdf=df.groupby('Cohort').agg('sum').transpose()
totdf=df.groupby('Cohort').agg(lambda x: len([i for i in x if i==i])).transpose()

countdf=pd.merge(posdf, totdf, right_index=True, left_index=True, suffixes=['_Count', '_Total'])

countdf['Phenotype']=countdf.index.to_list()

countdf=countdf.melt(id_vars='Phenotype')

countdf['column']=countdf.Cohort.str.split('_', expand=True)[1]
countdf['Cohort']=countdf.Cohort.str.split('_', expand=True)[0]

countdf=pd.pivot(countdf, index=['Phenotype', 'Cohort'], columns='column', values='value')
countdf.reset_index(inplace=True)

countdf=countdf[((countdf.Cohort=='DD (adults)') & (countdf.Phenotype.str.contains('Questionnaire'))) | ((countdf.Cohort=='DD (children)') & (~countdf.Phenotype.str.contains('Questionnaire')))]

countdf['Without']=countdf.Total-countdf.Count
countdf[['Count', 'Without', 'Total']]=countdf[['Count', 'Without', 'Total']].astype(int)

countdf=countdf[['Phenotype', 'Cohort', 'Count', 'Without', 'Total']]
countdf.Phenotype=countdf.Phenotype.str.split(' \\(', expand=True)[0]

# Clean up and save
countdf.sort_values(by=['Cohort', 'Phenotype'], inplace=True)

countdf.to_csv(OUTPUT_TAB, index=False)
