import pandas as pd
# Note that you will also need openpyxl available in your environment to read Excel files

# Filter CNV calls from Sanders et al. Neuron 2015 (PubMed: https://pubmed.ncbi.nlm.nih.gov/26402605/)

# Input and output files
TABS2='path/to/Sanders_Neuron_2015/Table_S2.xlsx'
TABS3='path/to/Sanders_Neuron_2015/Table_S3.xlsx'
# These files are Supplemental Tables S2 and S3 from Sanders et al. Neuron 2015
# They represent de novo CNVs and inherited CNVs, respectively, identified in SSC families from microarray data
OUTPUT_BED="/path/to/output/table.bed"

# Load in data
tab2=pd.read_excel(TABS2)
tab3=pd.read_excel(TABS3)

# Adjust formatting for consistency across tables
tab3['PatientID']=tab3.patientID
tab2['Inheritance']='de novo'

# Combine files
df=pd.concat([tab2, tab3])

print(df)
print(df.columns.to_list())

# Remove any calls <50kb
df=df[df.Size>=50000]

# Remove any inherited calls that encompass less than 5 SNPs (unless de novo)
df=df[(df.SNPs>=5) | (df.Inheritance=='de novo')]

# Remove any calls without hg19 coordinates
df=df[df['Start(hg19)']!='.']

# Save only needed columns
df=df[['Chr', 'Start(hg19)', 'End(hg19)', 'PatientID', 'Del/Dup', 'Inheritance']]
df.to_csv(OUTPUT_BED, sep='\t', index=False)
