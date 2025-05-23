import pandas as pd
import numpy as np

# Make lists of secondary variant genes for probands

# Input and output files
SNVS="/path/to/SNV/variants.csv" # Use the output of script 1_Variant calling/DD cohort/2_SNV_annotation/coding_annotations/14_loeuf_scores.py
CNVS="/path/to/CNV/calls.txt" # Use the output from script 1_Variant calling/DD cohort/3_CNV_calling_annotation/merge_all_cnvs/2_annotate_loeuf.py
STRS="/path/to/STR/calls.csv" # Use the output from script 1_Variant calling/DD cohort/4_STR_calling_annotation/11_annotate_loeuf.py
TABS1A="/path/to/TableS1A.csv" # Use the output of script 1_make_table_s1a.py

pheno_names={'Behavioral features (Child domain)':'Behavioral features', 'Psychiatric features (Child domain)':'Psychiatric features', 'Nervous System Abnormalities (Child domain)':'Nervous System Abnormalities',
			'Congenital Anomalies (Child domain)':'Congenital Anomalies', 'Growth/Skeletal Defects (Child domain)':'GrowthSkeletal Defects'}

PROBAND_LIST="/path/to/output/list/of/probands.csv"
OUTPUT_DIR="/path/to/output/directory/" # Gene lists will be saved into two directories within this folder:
# proband_genes/ - will contain genes across all probands broken up by variant type
# proband_phenotypes/ - will contain genes in probands with each phenotypic domain, with a sub-directory for each phenotypic domain and separate files for each variant type
# Load data
df=pd.read_csv(TABS1A)

# Restrict to probands with WGS data
df=df[(df.Relationship=='Proband') & (df.WGS=='X')]

# Save list of probands
df[['Sample']].to_csv(PROBAND_LIST, index=False)

# SNVs
snvs=pd.read_csv(SNVS)
snvs=snvs[snvs.Sample.isin(df.Sample.to_list())]
snvs['Gene_id_']=snvs['Gene_id_'].str.split(';')
snvs=snvs.explode('Gene_id_')

for mt in ['missense', 'lof', 'splice']:
	sub=snvs[snvs.Mut_type==mt]
	genes=sorted(list(sub['Gene_id_'].unique()))
	
	# Save genes to list
	name=mt.title()
	if mt=='lof':
		name='LOF'
	with open(f'{OUTPUT_DIR}/proband_genes/'+name+'.csv', 'w') as outfile:
		outfile.write('\n'.join(genes))
		outfile.write('\n')

	for key in pheno_names.keys():
		subk=sub[sub.Sample.isin(df[df[key]>0].Sample.to_list())]
		genes=sorted(list(subk['Gene_id_'].unique()))
		with open(f'{OUTPUT_DIR}/proband_phenotypes/'+pheno_names[key]+'/'+name+'.csv', 'w') as outfile:
			outfile.write('\n'.join(genes))
			outfile.write('\n')

	# LF
	sub=sub[sub.LOEUF<=0.35]
	genes=sorted(list(sub['Gene_id_'].unique()))
	with open(f'{OUTPUT_DIR}/proband_genes/'+name+' (LF).csv', 'w') as outfile:
		outfile.write('\n'.join(genes))
		outfile.write('\n')

	for key in pheno_names.keys():
		subk=sub[sub.Sample.isin(df[df[key]>0].Sample.to_list())]
		genes=sorted(list(subk['Gene_id_'].unique()))
		with open(f'{OUTPUT_DIR}/proband_phenotypes/'+pheno_names[key]+'/'+name+' (LF).csv', 'w') as outfile:
			outfile.write('\n'.join(genes))
			outfile.write('\n')
	
genes=sorted(list(snvs['Gene_id_'].unique()))
with open(f'{OUTPUT_DIR}/proband_genes/All coding SNVs.csv', 'w') as outfile:
	outfile.write('\n'.join(genes))
	outfile.write('\n')
for key in pheno_names.keys():
	sub=snvs[snvs.Sample.isin(df[df[key]>0].Sample.to_list())]
	genes=sorted(list(sub['Gene_id_'].unique()))
	with open(f'{OUTPUT_DIR}/proband_phenotypes/'+pheno_names[key]+'/All coding SNVs.csv', 'w') as outfile:
		outfile.write('\n'.join(genes))
		outfile.write('\n')

genes=sorted(list(snvs[snvs.LOEUF<=0.35]['Gene_id_'].unique()))
with open(f'{OUTPUT_DIR}/proband_genes/All coding SNVs (LF).csv', 'w') as outfile:
	outfile.write('\n'.join(genes))
	outfile.write('\n')
for key in pheno_names.keys():
	sub=snvs[(snvs.Sample.isin(df[df[key]>0].Sample.to_list())) & (snvs.LOEUF<=0.35)]
	genes=sorted(list(sub['Gene_id_'].unique()))
	with open(f'{OUTPUT_DIR}/proband_phenotypes/'+pheno_names[key]+'/All coding SNVs (LF).csv', 'w') as outfile:
		outfile.write('\n'.join(genes))
		outfile.write('\n')

# CNVs
cnvs=pd.read_csv(CNVS, sep='\t')
cnvs=cnvs[cnvs.Sample.isin(df.Sample.to_list())]
cnvs['Gene_ID']=cnvs['Gene_ID'].str.split(' ')
cnvs=cnvs.explode('Gene_ID')
cnvs.replace('.', np.nan, inplace=True)
cnvs.LOEUF=cnvs.LOEUF.astype(float)

cnv_names={'DEL':'Genes del', 'DUP':'Genes dup'}
for cn in cnv_names.keys():
	sub=cnvs[cnvs['Type']==cn]
	genes=sorted(list(sub['Gene_ID'].unique()))
	
	with open(f'{OUTPUT_DIR}/proband_genes/'+cnv_names[cn]+'.csv', 'w') as outfile:
		outfile.write('\n'.join(genes))
		outfile.write('\n')
	
	for key in pheno_names.keys():
		subk=sub[sub.Sample.isin(df[df[key]>0].Sample.to_list())]
		genes=sorted(list(subk['Gene_ID'].unique()))
		with open(f'{OUTPUT_DIR}/proband_phenotypes/'+pheno_names[key]+'/'+cnv_names[cn]+'.csv', 'w') as outfile:
			outfile.write('\n'.join(genes))
			outfile.write('\n')

	# LF
	sub=sub[sub.LOEUF<=0.35]
	genes=sorted(list(sub['Gene_ID'].unique()))
	with open(f'{OUTPUT_DIR}/proband_genes/'+cnv_names[cn]+' (LF).csv', 'w') as outfile:
		outfile.write('\n'.join(genes))
		outfile.write('\n')

	for key in pheno_names.keys():
		subk=sub[sub.Sample.isin(df[df[key]>0].Sample.to_list())]
		genes=sorted(list(subk['Gene_ID'].unique()))
		with open(f'{OUTPUT_DIR}/proband_phenotypes/'+pheno_names[key]+'/'+cnv_names[cn]+' (LF).csv', 'w') as outfile:
			outfile.write('\n'.join(genes))
			outfile.write('\n')

# STRs
strs=pd.read_csv(STRS, index_col=0)
strs=strs[strs.Sample.isin(df.Sample.to_list())]
strs['Gene_id_']=strs['Gene_id_'].str.split(';')
strs=strs.explode('Gene_id_')

genes=sorted(list(strs['Gene_id_'].unique()))
with open(f'{OUTPUT_DIR}/proband_genes/STRs.csv', 'w') as outfile:
	outfile.write('\n'.join(genes))
	outfile.write('\n')
for key in pheno_names.keys():
	sub=strs[strs.Sample.isin(df[df[key]>0].Sample.to_list())]
	genes=sorted(list(sub['Gene_id_'].unique()))
	with open(f'{OUTPUT_DIR}/proband_phenotypes/'+pheno_names[key]+'/STRs.csv', 'w') as outfile:
		outfile.write('\n'.join(genes))
		outfile.write('\n')

genes=sorted(list(strs[strs.LOEUF<=0.35]['Gene_id_'].unique()))
with open(f'{OUTPUT_DIR}/proband_genes/STRs (LF).csv', 'w') as outfile:
	outfile.write('\n'.join(genes))
	outfile.write('\n')
for key in pheno_names.keys():
	sub=strs[(strs.Sample.isin(df[df[key]>0].Sample.to_list())) & (strs.LOEUF<=0.35)]
	genes=sorted(list(sub['Gene_id_'].unique()))
	with open(f'{OUTPUT_DIR}/proband_phenotypes/'+pheno_names[key]+'/STRs (LF).csv', 'w') as outfile:
		outfile.write('\n'.join(genes))
		outfile.write('\n')

# Concat all files into single one
df=pd.DataFrame()
vars=['All coding SNVs', 'Missense', 'LOF', 'Splice', 'Genes del', 'Genes dup', 'STRs', 'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del (LF)', 'Genes dup (LF)', 'STRs (LF)']
for var in vars:
	vdf=pd.read_csv(f'{OUTPUT_DIR}/proband_genes/'+var+'.csv', header=None, names=['Gene'])
	df=pd.concat([df, vdf])
df.drop_duplicates(inplace=True)

# Save
df.to_csv(f'{OUTPUT_DIR}/proband_genes/All_variants.csv', index=None, header=None)

for key in pheno_names.keys():
	df=pd.DataFrame()
	for var in vars:
		vdf=pd.read_csv(f'{OUTPUT_DIR}/proband_phenotypes/'+pheno_names[key]+'/'+var+'.csv', header=None, names=['Gene'])
		df=pd.concat([df, vdf])
	df.drop_duplicates(inplace=True)
	df.to_csv(f'{OUTPUT_DIR}/proband_phenotypes/'+pheno_names[key]+'/All_variants.csv', index=None, header=None)
