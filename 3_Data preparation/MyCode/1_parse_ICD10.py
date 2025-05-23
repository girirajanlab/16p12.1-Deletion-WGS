import pandas as pd

# Parse ICD10 data

# Input and output files
ICD10_DATA="/paht/to/input/ICD10/data.csv" # This file contains medical record billing codes for the MyCode cohort

CODING_19="/path/to/UK_Biobank/data_coding/19.tsv" # UKB ICD10 phenotypes are encoded using Data-coding 19 - we will use this file to determine the hierarchy for ICD10 codes in MyCode

CHAPTER_OUTPUT="/path/to/ICD10/Chapter/output/table.csv"
PHENOTYPE_OUTPUT="/path/to/ICD10/specific/phenotype/output/table.csv"

# Load data
icd=pd.read_csv(ICD10_DATA)

# Restrict to INPATIENT codes
icd=icd[icd['type']!='OUTPATIENT']

# Restrict to ICD10 codes (remove ICD9 codes)
icd=icd[~icd.ICD10.isnull()]

# Annotate using tiers from the hierarchy from the UK Biobank coding 19
# Map coding to ICD10 codes
coding=pd.read_csv(CODING_19_PATH, sep='\t')
coding['ICD10']=coding.meaning.str.split(' ', expand=True)[0]
# Annotate terms hierarchical categories
# Chapter (i.e. "Chapter I Certain infectious and parasitic diseases")
# Block (i.e. "A00-A09 Intestinal infectious diseases")
# Sub-block (i.e. "A00 Cholera")
coding['label']='.'
coding.loc[coding.coding.str.contains('Chapter'), 'label']='Chapter'
coding.loc[coding.coding.str.contains('Block'), 'label']='Block'
coding.loc[~(coding.meaning.str.contains('\\.')) & (~coding.coding.str.contains('Block')) & (~coding.coding.str.contains('Chapter')), 'label']='Sub-block'

# Reorganize coding as hierarchy
hierarchy=nx.DiGraph()
hierarchy.add_nodes_from(coding.node_id.to_list())
# Add parent-child edges
parent=[tuple(r) for r in coding[coding.parent_id!=0][['parent_id', 'node_id']].to_numpy()]
hierarchy.add_edges_from(parent)

# Assign each node a chapter, block, and sub-block
chapters=[]
blocks=[]
subblocks=[]
for node in coding.node_id.to_list():
    preds=[node]+list(hierarchy.predecessors(node))
    while True:
        old_preds=preds.copy()
        for p in old_preds:
            preds+=list(hierarchy.predecessors(p))
        preds=sorted(list(set(preds)))
        if len(preds)==len(old_preds):
            break
    
    for i in range(3):
        lst=[chapters, blocks, subblocks][i]
        lab=['Chapter', 'Block', 'Sub-block'][i]
        subdf=coding[(coding.node_id.isin(preds)) & (coding.label==lab)]
        if subdf.shape[0]>0:
            lst.append(subdf.node_id.to_list()[0])
        else:
            lst.append('.')
    
coding['Chapter']=chapters
coding['Block']=blocks
coding['Sub-block']=subblocks

# Assign chapters to MyCode sample ICD10 codes
icd['node']=icd.ICD10.map(dict(zip(coding.ICD10.to_list(), coding.node_id.to_list())))
icd['Sub-block']=icd.ICD10.map(dict(zip(coding.ICD10.to_list(), coding['Sub-block'].to_list())))
icd['Block']=icd.ICD10.map(dict(zip(coding.ICD10.to_list(), coding.Block.to_list())))
icd['Chapter']=icd.ICD10.map(dict(zip(coding.ICD10.to_list(), coding.Chapter.to_list())))

# Save a sample x chapter table
chaps=list(icd.Chapter.unique())
new_df=pd.DataFrame(0, index=list(icd.Sample.unique()), columns=[g+'_'+str(int(float(c))) for c in chaps])
for c in chaps:
	new_df.loc[new_df.index.isin(icd[icd.Chapter==c]['Sample'].to_list()), g+'_'+str(int(float(c)))]=1

new_df['Sample']=new_df.index
new_df=new_df[['Sample']+[g+'_'+str(int(float(c))) for c in chaps]]

# Save
new_df.to_csv(CHAPTER_OUTPUT, index=False)

# Select data to be used for comparing the prevalence of disorders between cohorts
# ICD10 codes we need:
# Node ID - meaning [type]
# 29110 - F32 Depressive episode [Sub-block]
# 29940 - F51 Nonorganic sleep disorders [Sub-block]
# R45.86 [node - not available]
# 29420 - F41 Other anxiety disorders [Sub-block]
# 27560 - F10.2 Dependence syndrome (Alcohol) [node]
# [27670, 27780, 27890, 28000, 28110, 28220, 28440, 28550] -  F[11-16, 18-19].2 Dependence syndrome (Drugs) [node]
# 780 - F20-29 Schizophrenia, schizotypal and delusional disorders [Block]
nodes={27560:'addiction', 27670:'addiction', 27780:'addiction', 27890:'addiction',
	28000:'addiction', 28110:'addiction', 28220:'addiction', 28440:'addiction', 28550:'addiction'}
subblock={29940:'sleep', 29110:'depression', 29420:'anxiety'}
block={780:'psychosis'}

icd_interp=pd.DataFrame(0, index=list(icd.Sample.unique()), columns=['sleep', 'addiction', 'depression', 'anxiety', 'psychosis'])

for n in nodes.keys():
	samps=icd[icd.node==str(n)]['Sample'].to_list()
	icd_interp.loc[samps, nodes[n]]=1
for sb in subblock.keys():
	samps=icd[icd['Sub-block']==str(sb)]['Sample'].to_list()
	icd_interp.loc[samps, subblock[sb]]=1
samps=icd[icd.Block=='780']['Sample'].to_list()
icd_interp.loc[samps, 'psychosis']=1

# Save
icd_interp.to_csv(PHENOTYPE_OUTPUT, index=False)

