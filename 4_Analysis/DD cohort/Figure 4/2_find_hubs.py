import pandas as pd
import numpy as np
import networkx as nx
from collections import Counter
import os

# Create individualized networks for probands and pathways
# Each network will contain only the 16p12.1 deletion genes, second hit genes, and any pathway genes

# Input and output files
STRING="/path/to/annotated/string/data.csv.gz" # Use the output of script 2_Analysis preparation\STRING_Network\1_parse_string.py
BIOMART_MAP="/path/to/gene/protein/map.txt" # Gene and protein stable IDs were downloaded from BioMart: https://grch37.ensembl.org/biomart/martview/e5b591d25b647e380ce50cd8155adbdb
GENE_LIST_DIR="/path/to/proband/gene/list/directory" # Use the output of script 3_Data preparation\DD cohort\3_genelist_by_proband.py

OUTPUT_GML_DIR="/path/to/output/GML/directory/" # A directory for output network graphml files per proband. GML files for three probands were plotted in Cytoscape to produce Fig 4D
OUTPUT_CONNECTORS="/path/to/output/connector/gene/file.csv" # A dataframe containing all connector genes and whether they act as hub genes

# Load STRING relationships
string=pd.read_csv(STRING, compression='gzip')
string.new_combined_score=string.new_combined_score/1000
string=string[string.new_combined_score>=0.9]
string['weight']=1/string.new_combined_score

# Load BioMART map
biomart_map=pd.read_csv(BIOMART_MAP, sep='\t')
biomart_map.index=biomart_map['Gene stable ID']
name_map=biomart_map['Gene name'].to_dict()
# Update name for MOSMO
name_map['ENSG00000185716']='MOSMO'

# Load 16p12 genes
genes_16p12={'UQCRC2':'ENSG00000140740',
				'PDZD9':'ENSG00000155714',
				'MOSMO':'ENSG00000185716',
				'VWA3A':'ENSG00000175267',
				'EEF2K':'ENSG00000103319',
				'POLR3E':'ENSG00000058600',
				'CDR2':'ENSG00000140743'}

# Load probands
probands=[i[0].split('/')[8] for i in os.walk(GENE_LIST_DIR) if i[0].split('/')[8]!='']

# Create networks for each proband

# Helper functions
# Make network
def make_network(df):
	G = nx.Graph()
	gt_palette={'16p12.1':'#B4654A', 'Second hit':'#3A6F71', 'Connector':'#6E6E6E'}
	for gt in ['16p12.1', 'Second hit', 'Connector']:
		subdf1=df[df.Gene1_genetype==gt][['Gene1', 'Gene1_Symbol']].copy()
		subdf2=df[df.Gene2_genetype==gt][['Gene2', 'Gene2_Symbol']].copy()
		subdf1.columns=['Gene', 'Gene_Symbol']
		subdf2.columns=['Gene', 'Gene_Symbol']
		subdf=pd.concat([subdf1, subdf2])
		subdf.drop_duplicates(inplace=True)
		G.add_nodes_from(subdf[['Gene', 'Gene_Symbol']].apply(lambda row: (row.Gene, {'Symbol':row.Gene_Symbol, 'color':gt_palette[gt]}), axis=1).to_list())
	G.add_edges_from(df[['Gene1', 'Gene2', 'weight']].apply(lambda row: (row.Gene1, row.Gene2, {'weight':row.weight}), axis=1).to_list())
	return G
# Restrict network to the shortest paths with 16p12.1 genes and second hits
def restrict_network(G, sh_genes, genes_16p12=genes_16p12):
	fh_genes=list(genes_16p12.values())
	keep_nodes=[]
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
			fg_nodes+=sp
		fgdf=pd.DataFrame.from_dict(dict(Counter(fg_nodes)), orient='index', columns=['Shortest_pathways'])
		fgdf['Primary_gene']=fg
		fgdf['Connector_gene']=fgdf.index.to_list()
		outdf=pd.concat([outdf, fgdf])
		keep_nodes+=fg_nodes
	adf=pd.DataFrame.from_dict(dict(Counter(keep_nodes)), orient='index', columns=['Shortest_pathways'])
	adf['Primary_gene']='All'
	adf['Connector_gene']=adf.index.to_list()
	outdf=pd.concat([outdf, adf])
	keep_nodes=list(set(keep_nodes))
	subG=nx.Graph(G.subgraph(keep_nodes))
	return subG, outdf

outdf=pd.DataFrame()
sh_genes_all=[]
for pro in probands:
	g16p12=list(genes_16p12.values())
	sh_genes=pd.read_csv(f'{GENE_LIST_DIR}/{pro}/All_variants.csv', header=None, names=['Gene_id']).Gene_id.to_list()
	sh_genes_all+=sh_genes
	
	# Make a network of just the 16p12 genes and second hit genes
	substring=string.copy()
	for g in ['Gene1', 'Gene2']:
		substring[g+'_genetype']='Connector'
		substring.loc[substring[g].isin(sh_genes), g+'_genetype']='Second hit'
		substring.loc[substring[g].isin(g16p12), g+'_genetype']='16p12.1'
	
	# Make network
	G=make_network(substring)
	# Restrict network to only pathways that contain 16p12.1 del genes and a second hit
	subG, countdf=restrict_network(G, sh_genes)
	countdf['Sample']=pro
	
	nodes=list(subG.nodes)
	size_dict=dict(zip(countdf[countdf.Primary_gene=='All'].Connector_gene.to_list(), countdf[countdf.Primary_gene=='All'].Shortest_pathways.to_list()))
	sizes=list(size_dict.values())
	top_95=10
	if len(sizes)>0:
		top_95=np.percentile(sizes, 95)
	
	countdf['Hub']=(countdf.Shortest_pathways>=top_95)
	# Define separate thresholds for hubs in each primary variant
	pvs=list(countdf.Primary_gene.unique())
	for pv in pvs:
		if pv=='All':
			continue
		sizes=countdf[countdf.Primary_gene==pv].Shortest_pathways.to_list()
		top_95=10
		if len(sizes)>0:
			top_95=np.percentile(sizes, 95)
		if top_95<5: # If top 5% is 4 or less, will include every gene on pathways to a single second hit from all four 16p12.1 genes
			top_95=5
		countdf.loc[(countdf.Primary_gene==pv) & (countdf.Shortest_pathways>=top_95), 'Hub']=True
	outdf=pd.concat([outdf, countdf])
	
	name_dict={}
	top_95=10
	sizes=list(size_dict.values())
	if len(sizes)>0:
		top_95=np.percentile(sizes, 95)
	for n in nodes:
		if n not in size_dict.keys():
			size_dict[n]=0
		# Give genes in the top 95% of connectivity names for graphing
		if size_dict[n]<top_95:
			name_dict[n]=''
		else:
			name_dict[n]=name_map[n]
	nx.set_node_attributes(subG, name_dict, 'Node_name')
	nx.set_node_attributes(subG, size_dict, 'Size')
	
	# Save a GraphML file of the network for Cytoscape
	gt_palette={'16p12.1':'#B4654A', 'Second hit':'#3A6F71', 'Connector':'#6E6E6E'}
	for i in range(2):
		geneset=[g16p12, sh_genes][i]
		genetype=['16p12.1', 'Second hit'][i]
		color=gt_palette[genetype]
		missing_genes=[i for i in geneset if i not in list(subG.nodes)]
		if len(missing_genes)>0:
			subG.add_nodes_from([(i, {'Symbol':name_map[i], 'color':color, 'Size':0}) for i in missing_genes])
	nx.write_graphml(subG, f'{OUTPUT_GML_DIR}{pro}_connector_network.graphml')

# Map gene names for convenience
outdf['Primary_gene_symbol']=outdf.Primary_gene.map(name_map)
outdf['Connector_gene_symbol']=outdf.Connector_gene.map(name_map)
outdf.fillna('', inplace=True)

# Save
outdf[['Primary_gene', 'Primary_gene_symbol', 'Connector_gene', 'Connector_gene_symbol', 'Sample', 'Shortest_pathways', 'Hub']].to_csv(OUTPUT_CONNECTORS, index=False)

# GMT files from this script were imported into Cytoscape for plotting
