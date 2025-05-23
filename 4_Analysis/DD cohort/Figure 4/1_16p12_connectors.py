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
GENE_LIST_DIR="/path/to/proband/gene/list/directory" # Use the output of script 3_Data preparation\DD cohort\3_genelist_by_proband.py
BACKGROUND="/path/to/string/background/genes.csv" # Use the output of script 2_Analysis preparation\STRING_Network\2_string_connectivity.py
GMT="/path/to/Human_GOBP_AllPathways_noPFOCR_no_GO_iea_March_01_2025_ensembl.gmt" # Downloaded from http://download.baderlab.org/EM_Genesets/current_release/Human/ensembl/Human_GOBP_AllPathways_noPFOCR_no_GO_iea_March_01_2025_ensembl.gmt

OUTPUT_PATH_DIR="/path/to/output/directory/for/paths/"
OUTPUT_GML="/path/to/output/file.graphml" # This output file was inmported into Cytoscape to create Fig 4A
UPSET_FIG="/path/to/output/upset/plot.pdf" # Upset plot presented in Fig 4B
TOP_CONNECTORS="/path/to/top/two/connector/genes.csv" # First and second degree connector genes for 16p12.1 deletion genes as presented in Table S3A
OUTPUT_GO="/path/to/GO/enrichment/output.csv" # GO enrichment data presented in Table S3C
OUTPUT_GMT="/path/to/output/GMT/files/" # A directory for GMT files for use in EnrichmentMap in Cytoscape. These files are used to create Fig 4D

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

# Helper functions
# Save paths between 16p12.1 genes and secondary variant genes
def save_paths(G, sh_genes, genes_16p12=genes_16p12):
	fh_genes=list(genes_16p12.values())
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
	gt_palette={'16p12.1':'#B4654A', 'Connector':'#6E6E6E'}
	for gt in ['16p12.1', 'Connector']:
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

for pro in probands:
	g16p12=list(genes_16p12.values())
	sh_genes=pd.read_csv(f'{GENE_LIST_DIR}/{pro}/All_variants.csv', header=None, names=['Gene_id']).Gene_id.to_list()
	
	# Get full pathways between 16p12.1 genes and all second hits in probands
	countdf=save_paths(G, sh_genes)
	countdf['Sample']=pro
	outdf=pd.concat([outdf, countdf])
	
# Drop any all-NA columns
outdf=outdf.dropna(axis=1, how='all')

# Save
outdf.to_csv(f'{OUTPUT_PATH_DIR}/16p12_proband_paths.csv', index=False)

# Save a version with mapped names
cols=[i for i in outdf.columns.to_list() if i not in ['Sample', 'Length']]
symbol_out=outdf.copy()
for c in cols:
	symbol_out[c]=symbol_out[c].map(name_map)
symbol_out.to_csv(f'{OUTPUT_PATH_DIR}/16p12_proband_paths_symbol.csv', index=False)

# Make networks showing the first connector leading from each 16p12.1 del gene
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

# Make an Upset plot of the data
conndf=outdf[[i for i in outdf.columns.to_list() if i not in ['Secondary', 'Length', 'Sample']]].melt(id_vars='Primary', var_name='Connection_degree', value_name='Gene')
conndf.Connection_degree=conndf.Connection_degree.str.split('C', expand=True)[1].astype(int)
conndf=conndf[~conndf.Gene.isnull()]

# Restrict to top 3 genes
conndf=conndf[conndf.Connection_degree<=2]

conndf['Primary Gene']=conndf.Primary.map({v: k for k, v in genes_16p12.items()})
conndf['Gene_deg']=conndf['Primary Gene']+' '+conndf.Connection_degree.astype(str)

conndf=conndf[['Gene_deg', 'Gene']]
conndf['placeholder']=1
conndf.drop_duplicates(inplace=True)

conndf=conndf.pivot(index='Gene', columns='Gene_deg', values='placeholder')
conndf.fillna(0, inplace=True)
conndf=conndf.astype(bool)

cols=conndf.columns.to_list()
conndf.reset_index(inplace=True)
conndf=conndf.groupby(cols).size()

plt.savefig(UPSET_FIG)

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

# Run gProfiler on unique connector genes for each 16p12.1 gene
gp=GProfiler(return_dataframe=True)

# Load background genes
backdf=pd.read_csv(BACKGROUND)
# Remove any genes from the background that have less than 2 connections (as they could not be connector genes, by definition)
backdf=backdf[backdf.Connections>=2]
background=backdf.GeneID.to_list()

go_res=pd.DataFrame()
for pg in primary_genes:
	# Unique connectors
	genes=list(df[df.Primary_gene_symbol==pg]['Connector_gene'].unique())
	other_genes=list(df[df.Primary_gene_symbol!=pg]['Connector_gene'].unique())
	uniq_genes=list(set(genes)-set(other_genes))
	
	go_result=gp.profile(organism='hsapiens', query=uniq_genes, background=background, sources=['GO:BP'])
	go_result['Primary_gene_symbol']=pg
	go_res=pd.concat([go_res, go_result])

# Save
go_res.to_csv(OUTPUT_GO, index=False)

# Create GMT files for Cytoscape
go_res=go_res[go_res.source=='GO:BP']

for gene in ['EEF2K', 'MOSMO', 'POLR3E', 'UQCRC2']:
	
	# Save GO GMT file for cytoscape
	outfile=open(f'{OUTPUT_GMT}/{gene}.gmt', 'w')
	gmt=open(GMT, 'r')
	for line in gmt.readlines():
		anno=line.split('\t')[0]
		if "%GOBP%" not in anno:
			continue
		term='GO:'+anno.split('GO:')[1]
		if term in go_res[go_res.Primary_gene_symbol==gene].native.to_list():
			outfile.write(line)
	outfile.close()
	gmt.close()

# GMT files were imported into the EnrichmentMap plugin within Cytoscape to create Fig 4D