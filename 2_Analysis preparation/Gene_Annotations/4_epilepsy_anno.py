#!bin/python
import pandas as pd

# Annotate with epilepsy genes derived from Wang et al. Seizure 2017
df = pd.read_csv('intermediate_tables/3_szdb_annotations.csv')

epi =pd.read_excel('Data_Files/Wang_Seizure_2017_epilepsy_gene_list.xlsx', header = None, names = ['Gene', 'index'])
    
df['Wang_Epilepsy']=0
df.loc[df.gene_symbol.isin(epi.Gene.to_list()), 'Wang_Epilepsy']=1
print(df.Wang_Epilepsy.value_counts())

# Save to final file
df.to_csv('intermediate_tables/4_wang_epilepsy_annotations.csv', index = False)