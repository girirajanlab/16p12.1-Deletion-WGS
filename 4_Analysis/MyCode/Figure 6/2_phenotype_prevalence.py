import pandas as pd

# Calculate the frequencies of each phenotype in MyCode

# Input and output files
MYCODE_DATA="/path/to/MyCode/data/files.csv" # Use the output of script 3_Data preparation\MyCode\2_condense_data.py
OUTPUT_FREQ"/path/to/output/frequencies.csv"

# Load data
df=pd.read_csv(MYCODE_DATA)
df=df[['Sleep trouble (ICD10)', 'Addiction (ICD10)', 'Depression (ICD10)', 'Anxiety (ICD10)', 'Psychosis (ICD10)']]
df.columns=['Sleep trouble', 'Addiction', 'Depression', 'Anxiety', 'Psychosis']

# Convert to counts
posdf=pd.DataFrame(df.sum(axis=0))
posdf.columns=['Count']
totdf=pd.DataFrame((~df.isnull()).sum())
totdf.columns=['Total']

countdf=pd.merge(posdf, totdf, right_index=True, left_index=True)
countdf['Without']=countdf.Total-countdf.Count
countdf=countdf.astype(int)

countdf['Cohort']='MyCode'
countdf['Phenotype']=countdf.index.to_list()

countdf=countdf[['Phenotype', 'Cohort', 'Count', 'Without', 'Total']]

countdf.to_csv(OUTPUT_FREQ, index=False)
