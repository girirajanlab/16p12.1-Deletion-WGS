import pandas as pd
from collections import Counter

# Get connectovity of all genes in STRING DB

# Input and output files
GENE_ANNO="/path/to/gene/annotations.csv" # Use the output of script 2_Analysis preparation\Gene_Annotations\5_add_loeuf.py
STRING="/path/to/output/STRING/network.csv.gz" # Use the output of script 1_parse_string.py

CONNECTIVITY_OUT="/path/to/output/network/connectivity/file.csv"

# Load all genes
background=pd.read_csv(GENE_ANNO)

# Filter the background gene list to include only those present in STRING
string=pd.read_csv(STRING, compression='gzip')
string_genes=string.Gene1.to_list()+string.Gene2.to_list()
background=background[background.GeneID.isin(string_genes)]

# Filter STRING to highest confidence interactions
string.new_combined_score=string.new_combined_score/1000
string=string[string.new_combined_score>=0.9]

# Get a counter of the total number of connections a gene is involved in
connections=dict(Counter(string.Gene1.to_list()+string.Gene2.to_list()))
background['Connections']=background.GeneID.map(connections)
background.fillna(0, inplace=True)
background['Connections']=background.Connections.astype(int)

# Assign connectivity quartiles
background['Quartile']=0
quants=[0.25, 0.5, 0.75, 1]
for q in reversed(quants):
	quantile=background.Connections.quantile(q=q)
	background.loc[background.Connections<=quantile, 'Quartile']=q

# Save
background.to_csv(CONNECTIVITY_OUT, index=False)