import pandas as pd
import numpy as np

import networkx as nx

from collections import Counter
from gprofiler import GProfiler

from upsetplot import plot
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype']=42

import os

# Identify all connector genes across all probands

# Input and output files
STRING="/path/to/annotated/string/data.csv.gz" # Use the output of script 2_Analysis preparation\STRING_Network\1_parse_string.py
BIOMART_MAP="/path/to/gene/protein/map.txt" # Gene and protein stable IDs were downloaded from BioMart: https://grch37.ensembl.org/biomart/martview/e5b591d25b647e380ce50cd8155adbdb
GENE_LIST_DIR="/path/to/proband/gene/list/directory" # Use the output directory of script 3_Data preparation\Searchlight\2_make_genelists.py for 16p11.2 deletion individual proband gene lists

OUTPUT_PATH_DIR="/path/to/output/directory/for/paths/"
OUTPUT_GML="/path/to/output/file.graphml" # This output file was inmported into Cytoscape to create Fig 4B
TOP_CONNECTORS="/path/to/top/two/connector/genes.csv" # First and second degree connector genes for 16p11.2 deletion genes

# Load STRING relationships
string=pd.read_csv(STRING, compression='gzip')
string.new_combined_score=string.new_combined_score/1000
string=string[string.new_combined_score>=0.9]
string['weight']=1/string.new_combined_score

# Load BioMART map
biomart_map=pd.read_csv(BIOMART_MAP, sep='\t')
biomart_map.index=biomart_map['Gene stable ID']
name_map=biomart_map['Gene name'].to_dict()

# Load 16p12 genes
genes_16p11={'BOLA2':'ENSG00000183336',
			'SLX1B':'ENSG00000181625',
			'SULT1A4':'ENSG00000213648',
			'SPN':'ENSG00000197471',
			'QPRT':'ENSG00000103485',
			'C16orf54':'ENSG00000185905',
			'ZG16':'ENSG00000174992',
			'MAZ':'ENSG00000103495',
			'PRRT2':'ENSG00000167371',
			'PAGR1':'ENSG00000185928',
			'MVP':'ENSG00000013364',
			'CDIPT':'ENSG00000103502',
			'SEZ6L2':'ENSG00000174938',
			'ASPHD1':'ENSG00000174939',
			'KCTD13':'ENSG00000174943',
			'TMEM219':'ENSG00000149932',
			'TAOK2':'ENSG00000149930',
			'HIRIP3':'ENSG00000149929',
			'INO80E':'ENSG00000169592',
			'DOC2A':'ENSG00000149927',
			'C16orf92':'ENSG00000167194',
			'FAM57B':'ENSG00000149926',
			'ALDOA':'ENSG00000149925',
			'SLC7A5P1':'ENSG00000260727',
			'PPP4C':'ENSG00000149923',
			'TBX6':'ENSG00000149922',
			'YPEL3':'ENSG00000090238',
			'GDPD3':'ENSG00000102886',
			'MAPK3':'ENSG00000102882',
			'CORO1A':'ENSG00000102879',
			'SLX1A':'ENSG00000132207',
			'SULT1A3':'ENSG00000261052'}

# Load probands
probands=[i[0].split('/')[8] for i in os.walk(GENE_LIST_DIR) if i[0].split('/')[8]!='']

# Helper functions
# Save paths between 16p12.1 genes and secondary variant genes
def save_paths(G, sh_genes, genes_16p11=genes_16p11):
	fh_genes=list(genes_16p11.values())
	outdf=pd.DataFrame()
	for fg in fh_genes:
		fg_nodes=[]
		if not G.has_node(fg):
			continue
		for sh in sh_genes:
			if not G.has_node(sh):
				continue
			if not nx.has_path(G, fg, sh):
				continue
			sp=nx.shortest_path(G, source=fg, target=sh, weight='weight')
			add_space=[np.nan]*(17-len(sp))
			fg_nodes.append(sp[0:(len(sp)-1)]+add_space+[sp[-1], len(sp)])
		fgdf=pd.DataFrame(fg_nodes, columns=['Primary']+[f'C{i}' for i in range(1, 16)]+['Secondary', 'Length'])
		outdf=pd.concat([outdf, fgdf])
	return outdf
# Make network
def make_network(df):
	G = nx.Graph()
	gt_palette={'16p11.2':'#4A9AE0', 'Connector':'#6E6E6E'}
	for gt in ['16p11.2', 'Connector']:
		subdf1=df[df.Gene1_genetype==gt][['Gene1', 'Gene1_Symbol']].copy()
		subdf2=df[df.Gene2_genetype==gt][['Gene2', 'Gene2_Symbol']].copy()
		subdf1.columns=['Gene', 'Gene_Symbol']
		subdf2.columns=['Gene', 'Gene_Symbol']
		subdf=pd.concat([subdf1, subdf2])
		subdf.drop_duplicates(inplace=True)
		G.add_nodes_from(subdf[['Gene', 'Gene_Symbol']].apply(lambda row: (row.Gene, {'Symbol':row.Gene_Symbol, 'color':gt_palette[gt]}), axis=1).to_list())
	G.add_edges_from(df[['Gene1', 'Gene2', 'weight']].apply(lambda row: (row.Gene1, row.Gene2, {'weight':row.weight}), axis=1).to_list())
	return G

G = nx.Graph()
G.add_nodes_from(string.Gene1.to_list()+string.Gene2.to_list())
G.add_edges_from(string[['Gene1', 'Gene2', 'weight']].apply(lambda row: (row.Gene1, row.Gene2, {'weight':row.weight}), axis=1).to_list())

g16p11=list(genes_16p11.values())

for pro in probands:
	sh_genes=pd.read_csv(f'{GENE_LIST_DIR}/{pro}/All_variants.csv', header=None, names=['Gene_id']).Gene_id.to_list()
	
	# Get full pathways between 16p12.1 genes and all second hits in probands
	countdf=save_paths(G, sh_genes)
	countdf['Sample']=pro
	outdf=pd.concat([outdf, countdf])
	
# Drop any all-NA columns
outdf=outdf.dropna(axis=1, how='all')

# Save
outdf.to_csv(f'{OUTPUT_PATH_DIR}/16p11_proband_paths.csv', index=False)

# Save a version with mapped names
cols=[i for i in outdf.columns.to_list() if i not in ['Sample', 'Length']]
symbol_out=outdf.copy()
for c in cols:
	symbol_out[c]=symbol_out[c].map(name_map)
symbol_out.to_csv(f'{OUTPUT_PATH_DIR}/16p11_proband_paths_symbol.csv', index=False)

# Make networks showing the first connector leading from each 16p11.2 del gene
gene_list=[i for i in list(outdf.Primary.unique())+list(outdf.C1.unique()) if i==i]
substring=string[(string.Gene1.isin(gene_list)) & (string.Gene2.isin(gene_list))].copy()
for g in ['Gene1', 'Gene2']:
		substring[g+'_genetype']='Connector'
		substring.loc[substring[g].isin(list(outdf.Primary.unique())), g+'_genetype']='16p12.1'
G=make_network(substring)

# Annotate genes with size based on the number of pathways they are present in
size_dict=dict(Counter(outdf.Primary.to_list()+outdf.C1.to_list()))
nx.set_node_attributes(G, size_dict, 'Size')

# Save graphml for plotting in Cytoscape
nx.write_graphml(G, OUTPUT_GML)

# Reformat connectors as presented in Table S3A
df=outdf[['Primary', 'C1', 'C2']].copy()
df.drop_duplicates(inplace=True)
df=df.melt(id_vars='Primary')
df.drop_duplicates(inplace=True)
df['Degree']=df.variable.map({'C1':'1st degree', 'C2':'2nd degree'})

# Drop missing
df=df[(~df.Primary.isnull()) & (~df.value.isnull())]

# Map names
df['Primary_symbol']=df.Primary.map(name_map)

df['Primary_degree']=df.Primary_symbol+' '+df.Degree
df['placeholder']=True

df=df[['value', 'Primary_degree', 'placeholder']]
df.drop_duplicates(inplace=True)
df.reset_index(inplace=True, drop=True)

df=df.pivot(index='value', columns='Primary_degree', values='placeholder')
df.fillna(False, inplace=True)

# Clean up
df['Ensembl ID']=df.index.to_list()
df['Gene Symbol']=df['Ensembl ID'].map(name_map)

# Save
df.to_csv(TOP_CONNECTORS, index=False)
