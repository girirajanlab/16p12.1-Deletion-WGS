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
GENE_LIST_DIR="/path/to/proband/gene/list/directory" # Use the output directory of script 3_Data preparation\Searchlight\2_make_genelists.py
GENE_ANNO="/path/to/gene/annotations.csv" # Use the output of script 2_Analysis preparation\Gene_Annotations\5_add_loeuf.py
GMT="/path/to/Human_GOBP_AllPathways_noPFOCR_no_GO_iea_March_01_2025_ensembl.gmt" # Downloaded from http://download.baderlab.org/EM_Genesets/current_release/Human/ensembl/Human_GOBP_AllPathways_noPFOCR_no_GO_iea_March_01_2025_ensembl.gmt

GO_OUTPUT="/path/to/output/GO/KEGG/enrichments.csv" # Data presented in Table S6D
GMT_OUTPUT_DIR="/path/to/GMT/output/directory/" # A directory for GMT files for use in EnrichmentMap in Cytoscape. These files are used to create Fig 7D

# Run gProfiler
cohorts=['16p11.2 deletion', '16p11.2 duplication']
background=sorted(list(set(pd.read_csv(GENE_ANNO).GeneID.to_list())))

gp=GProfiler(return_dataframe=True)
outdf=pd.DataFrame(columns=['native', 'name'])
for co in cohorts:
	genes=sorted(list(set(pd.read_csv(f'{GENE_LIST_DIR}/{co}_genes.csv', header=None, names=['Gene']).Gene.to_list())))
	
	go_result=gp.profile(organism='hsapiens', query=genes, background=background, sources=['GO:BP', 'KEGG'])
	go_result['Cohort']=co
	
	# Save GMT file for cytoscape
	outfile=open(f'{GMT_OUTPUT_DIR}/{co}.gmt', 'w')
	gmt=open(GMT, 'r')
	for line in gmt.readlines():
		anno=line.split('\t')[0]
		if "%GOBP%" not in anno:
			continue
		term='GO:'+anno.split('GO:')[1]
		if term in go_result.native.to_list():
			outfile.write(line)
	outfile.close()
	gmt.close()	
	
	outdf=pd.concat([outdf, go_result])
# Save
outdf.to_csv(GO_OUTPUT, index=False)
