import pandas as pd

# Add annotations to STRING database download

# Input and Output files
STRING_DB="/path/to/STRING/DB/9606.protein.links.full.v12.0.onlyAB.txt.gz" # STRING interaction data was downloaded from https://stringdb-downloads.org/download/stream/protein.links.v12.0/9606.protein.links.v12.0.onlyAB.txt.gz
BIOMART_MAP="/path/to/gene/protein/map.txt" # Gene and protein stable IDs were downloaded from BioMart: https://grch37.ensembl.org/biomart/martview/e5b591d25b647e380ce50cd8155adbdb

OUTPUT="/path/to/output/STRING/network.csv.gz"

# Load STRING data
df=pd.read_csv(STRING_DB, sep=' ', compression='gzip')
df['protein1']=df.protein1.str.split('9606.', expand=True)[1]
df['protein2']=df.protein2.str.split('9606.', expand=True)[1]

# Map proteins to genes
biomart_map=pd.read_csv(BIOMART_MAP, sep='\t')
biomart_map=biomart_map[~biomart_map['Protein stable ID'].isnull()]
biomart_map.index=biomart_map['Protein stable ID']

protein_map=biomart_map['Gene stable ID'].to_dict()

df['Gene1']=df.protein1.map(protein_map)
df['Gene2']=df.protein2.map(protein_map)

df=df[~df[['Gene1', 'Gene2']].isnull().any(axis=1)]

name_map=biomart_map['Gene name'].to_dict()
df['Gene1_Symbol']=df.protein1.map(name_map)
df['Gene2_Symbol']=df.protein2.map(name_map)

# Re-calculate the interaction score to keep only certain associations
print(df)
cols=['fusion', 'cooccurence', 'homology', 'coexpression', 'coexpression_transferred', 'experiments', 'experiments_transferred', 'database', 'database_transferred']
df=df[['Gene1', 'Gene1_Symbol', 'protein1', 'Gene2', 'Gene2_Symbol', 'protein2']+cols+['combined_score']]
prior=0.041
for c in cols:
	df[c]=df[c]/1000
	df[f'{c}_calc']=df[c]
	df.loc[df[f'{c}_calc']<prior, f'{c}_calc']=prior
	df[f'{c}_noprior']=(df[f'{c}_calc']-prior)/(1-prior)
	df[f'{c}_subtot']=(1-df[f'{c}_noprior'])
df['s_tot_noprior']=1-(df[[f'{c}_subtot' for c in cols]].prod(axis=1))
df['new_combined_score']=df.s_tot_noprior+(prior*(1-df.s_tot_noprior))
df['combined_score']=df.combined_score/1000
df=df[['Gene1', 'Gene1_Symbol', 'protein1', 'Gene2', 'Gene2_Symbol', 'protein2']+cols+['combined_score', 'new_combined_score']]
df[cols+['combined_score', 'new_combined_score']]=(df[cols+['combined_score', 'new_combined_score']]*1000).astype(int)

# Save
df.to_csv(OUTPUT, index=False, compression='gzip')
