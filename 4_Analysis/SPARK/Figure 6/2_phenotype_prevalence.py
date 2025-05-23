import pandas as pd

# Calculate the frequencies of each phenotype in SPARK

# Input and output files
SPARK_DATA="/path/to/SPARK/data.csv" # Use the output of script For_GitHub\3_Data preparation\SPARK\1_condense_data.py
OUTPUT_FREQ"/path/to/output/frequencies.csv"

# Calculate frequencies
df=pd.read_csv(SPARK_DATA)
df=df[['Birth/pregnancy complications', 'Preterm birth', 'Microcephaly', 'Macrocephaly', 'Strabismus',
		'Seizures', 'Heart defects', 'Hearing loss', 'Vision problems', 'Feeding problems',
		'Obesity', 'ID/DD', 'Motor delay', 'Speech delay', 'Language disorder', 'Aide in school',
		'Learning disability', 'ASD', 'ADHD', 'OCD', 'Schizophrenia', 'BPD', 'Depression', 'Anxiety', 'Sleep trouble', 'PDD']]

# Convert to counts
posdf=pd.DataFrame(df.sum(axis=0))
posdf.columns=['Count']
totdf=pd.DataFrame((~df.isnull()).sum())
totdf.columns=['Total']

countdf=pd.merge(posdf, totdf, right_index=True, left_index=True)
countdf['Without']=countdf.Total-countdf.Count
countdf=countdf.astype(int)

# Clean up
countdf['Cohort']='SPARK'
countdf['Phenotype']=countdf.index.to_list()

countdf=countdf[['Phenotype', 'Cohort', 'Count', 'Without', 'Total']]

countdf.to_csv(OUTPUT_FREQ, index=False)
