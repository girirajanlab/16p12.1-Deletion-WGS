import pandas as pd
import numpy as np
from gprofiler import GProfiler
import re

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

plt.rcParams['pdf.fonttype']=42

# Perform gene ontology analysis using gProfiler

# Input and output files
GENE_LIST_DIR="/path/to/proband/gene/list/directory" # Use the output of script 3_Data preparation\DD cohort\2_make_genelists.py
GENE_ANNO="/path/to/gene/annotations.csv" # Use the output of script 2_Analysis preparation\Gene_Annotations\5_add_loeuf.py
GMT="/path/to/Human_GOBP_AllPathways_noPFOCR_no_GO_iea_March_01_2025_ensembl.gmt" # Downloaded from http://download.baderlab.org/EM_Genesets/current_release/Human/ensembl/Human_GOBP_AllPathways_noPFOCR_no_GO_iea_March_01_2025_ensembl.gmt
BRITE_08907="/path/to/KEGG/BRITE/HIERACHY/br08901.keg" # BRITE hierarchy file downloaded from https://www.kegg.jp/kegg-bin/download_htext?htext=br08901&format=htext&filedir=

GO_OUTPUT="/path/to/output/GO/KEGG/enrichments.csv" # Data presented in Table S4B
GMT_OUTPUT_DIR="/path/to/GMT/output/directory/" # A directory for GMT files for use in EnrichmentMap in Cytoscape. These files are used to create Fig 5C
KEGG_FIG="/path/to/output/KEGG/waterfall/plot.pdf" # Figure shown in Fig S4E

# Run gProfiler
child_domains=['Behavioral features', 'Psychiatric features', 'Nervous System Abnormalities', 'Congenital Anomalies', 'GrowthSkeletal Defects']
gp=GProfiler(return_dataframe=True)

# Load background genes
background=sorted(list(set(pd.read_csv(GENE_ANNO).GeneID.to_list())))

outdf=pd.DataFrame(columns=['native', 'name'])
for pheno in child_domains:
	genes=sorted(list(set(pd.read_csv(f'{GENE_LIST_DIR}/proband_phenotypes/{pheno}/All_variants.csv', header=None, names=['Gene']).Gene.to_list())))
	
	go_result=gp.profile(organism='hsapiens', query=genes, background=background, sources=['GO:BP', 'KEGG'])
	go_result['Phenotype']=pheno
	
	# Save GO GMT file for cytoscape
	outfile=open(f'{GMT_OUTPUT_DIR}/{pheno}.gmt', 'w')
	gmt=open(GMT, 'r')
	for line in gmt.readlines():
		anno=line.split('\t')[0]
		if "%GOBP%" not in anno:
			continue
		term='GO:'+anno.split('GO:')[1]
		if term in go_result[go_result.source=='GO:BP'].native.to_list():
			outfile.write(line)
	outfile.close()
	gmt.close()	
	
	outdf=pd.concat([outdf, go_result])
# Save
outdf.to_csv(GO_OUTPUT, index=False)

# Make a plot of GO terms in Cytoscape using the GMT files using the EnrichmentMap plugin

# Make a waterfall plot of KEGG terms
# Annotate KEGG terms with parents from KEGG hierarchy 
kegg_terms=outdf[outdf.source=='KEGG'].copy()
kegg_terms['ID']=kegg_terms.native.str.split('KEGG:', expand=True)[1]

hier_out=[]
with open(BRITE_08907, 'r') as hierarchy:
	a=''
	b=''
	for line in hierarchy.readlines():
		if line[0]=='A':
			a=line.split('>')[1].split('<')[0]
		if line[0]=='B':
			b=line.split('B ')[1].strip()
		if line[0]=='C':
			kegg_id=re.search(r"\d{5}", line).group()
			kegg_term=line.split(kegg_id+' ')[1].strip()
			hier_out.append([a, b, kegg_id, kegg_term])

brite_hierarchy=pd.DataFrame(hier_out, columns=['A', 'B', 'ID', 'Term'])
brite_hierarchy.index=brite_hierarchy.ID.to_list()

kegg_terms=pd.merge(kegg_terms, brite_hierarchy, on='ID', how='left')

# PLot terms as waterfall with colors
palette={"Metabolism":"#D88C8C",
		"Genetic Information Processing":"#7A9EAB",
		"Environmental Information Processing":"#D8B08C",
		"Cellular Processes":"#88B2A9",
		"Organismal Systems":"#B7A6D7",
		"Human Diseases":"#5A6B70",
		"Drug Development":"#D68BA3"}
kegg_terms['neglog10p']=(-1)*np.log10(kegg_terms.p_value)
max_p=kegg_terms.neglog10p.max()
max_x=kegg_terms.Phenotype.value_counts().max()

fig, axs=plt.subplots(ncols=2, nrows=len(child_domains), width_ratios=[5, 1],
						sharex=True, sharey=True,
						figsize=(6, 6))
for idx, cd in enumerate(child_domains):
	ax=axs[idx, 0]
	# Order terms
	subdf=kegg_terms[kegg_terms.Phenotype==cd].copy()
	subdf.sort_values(by='p_value', ascending=True, inplace=True)
	subdf.reset_index(inplace=True, drop=True)
	subdf['x']=subdf.index.to_list()
	
	print(subdf)
	
	# Plot
	sns.scatterplot(data=subdf, x='x', y='neglog10p', hue='A', palette=palette, ax=ax, legend=False)
	
	# Clean up figure
	ax.set_ylim(0, max_p+1)
	ax.set_xlim(-0.5, max_x+0.5)
	
	ax.set_ylabel('')
	ax.set_xlabel('')
	ax.set_title(cd)
	
	if idx==2:
		ax.set_ylabel('-log10(p value)')

# Add custom legend
legend_elem=[]
for key in palette:
	legend_elem.append(Line2D([0], [0], marker='o', color='white', markerfacecolor=palette[key], label=key, markersize=10))
axs[0, 1].legend(handles=legend_elem, loc='center left', fontsize=8)

for i in range(5):
	axs[i, 1].axis('off')

plt.tight_layout()
	
plt.savefig(KEGG_FIG)
plt.close()