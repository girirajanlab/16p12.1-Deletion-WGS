import pandas as pd
import numpy as np

import subprocess

# Make lists of secondary variant genes for probands

# Input and output files
SNVS="/path/to/SNV/variants.csv" # Use the output of script 1_Variant calling/DD cohort/2_SNV_annotation/coding_annotations/14_loeuf_scores.py
CNVS="/path/to/CNV/calls.txt" # Use the output from script 1_Variant calling/DD cohort/3_CNV_calling_annotation/merge_all_cnvs/2_annotate_loeuf.py
STRS="/path/to/STR/calls.csv" # Use the output from script 1_Variant calling/DD cohort/4_STR_calling_annotation/11_annotate_loeuf.py
TABS1A="/path/to/TableS1A.csv" # Use the output of script 1_make_table_s1a.py

pheno_names={'Behavioral features (Child domain)':'Behavioral features', 'Psychiatric features (Child domain)':'Psychiatric features', 'Nervous System Abnormalities (Child domain)':'Nervous System Abnormalities',
			'Congenital Anomalies (Child domain)':'Congenital Anomalies', 'Growth/Skeletal Defects (Child domain)':'GrowthSkeletal Defects'}

PROBAND_LIST="/path/to/output/list/of/probands.csv"
OUTPUT_DIR="/path/to/output/directory/" # Directories for each proband will be saved here with separate files for genes affected by each variant type

# Load data
df=pd.read_csv(TABS1A)
df=df[(df.Relationship=='Proband') & (df['Estonian Biobank Sample'].isnull())]
child_domains=['Behavioral features (Child domain)', 'Psychiatric features (Child domain)', 'Nervous System Abnormalities (Child domain)', 'Congenital Anomalies (Child domain)', 'Growth/Skeletal Defects (Child domain)']
df=df[['Sample', 'All coding SNVs', 'Genes del.']+child_domains]

probands=df.Sample.to_list()

for pro in probands:
	cmd=['mkdir', f'{OUTPUT_DIR}\\{pro}']
	subprocess.run(cmd, shell=True)

# SNVs
snvs=pd.read_csv(SNVS)
snvs=snvs[snvs.Sample.isin(df.Sample.to_list())]
snvs['Gene_id_']=snvs['Gene_id_'].str.split(';')
snvs=snvs.explode('Gene_id_')

for pro in probands:
	for mt in ['missense', 'lof', 'splice']:
		sub=snvs[(snvs.Mut_type==mt) & (snvs.Sample==pro)]
		genes=sorted(list(sub['Gene_id_'].unique()))
		
		# Save genes to list
		name=mt.title()
		if mt=='lof':
			name='LOF'
		with open(f'{OUTPUT_DIR}/{pro}/{name}.csv', 'w') as outfile:
			outfile.write('\n'.join(genes))
			outfile.write('\n')

		# LF
		sub=sub[sub.LOEUF<=0.35]
		genes=sorted(list(sub['Gene_id_'].unique()))
		with open(f'{OUTPUT_DIR}/{pro}/{name} (LF).csv', 'w') as outfile:
			outfile.write('\n'.join(genes))
			outfile.write('\n')

for pro in probands:
	genes=sorted(list(snvs[snvs.Sample==pro]['Gene_id_'].unique()))
	with open(f'{OUTPUT_DIR}/{pro}/All coding SNVs.csv', 'w') as outfile:
		outfile.write('\n'.join(genes))
		outfile.write('\n')

	genes=sorted(list(snvs[(snvs.Sample==pro) & (snvs.LOEUF<=0.35)]['Gene_id_'].unique()))
	with open(f'{OUTPUT_DIR}/{pro}/All coding SNVs (LF).csv', 'w') as outfile:
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
for pro in probands:
	for cn in cnv_names.keys():
		sub=cnvs[(cnvs['Type']==cn) & (cnvs.Sample==pro)]
		genes=sorted(list(sub['Gene_ID'].unique()))
		
		with open(f'{OUTPUT_DIR}/{pro}/{cnv_names[cn]}.csv', 'w') as outfile:
			outfile.write('\n'.join(genes))
			outfile.write('\n')

		# LF
		sub=sub[sub.LOEUF<=0.35]
		genes=sorted(list(sub['Gene_ID'].unique()))
		with open(f'{OUTPUT_DIR}/{pro}/{cnv_names[cn]} (LF).csv', 'w') as outfile:
			outfile.write('\n'.join(genes))
			outfile.write('\n')

# STRs
strs=pd.read_csv(STRS, index_col=0)
strs=strs[strs.Sample.isin(df.Sample.to_list())]
strs['Gene_id_']=strs['Gene_id_'].str.split(';')
strs=strs.explode('Gene_id_')

for pro in probands:
	genes=sorted(list(strs[strs.Sample==pro]['Gene_id_'].unique()))
	write_name=smap[pro]
	with open(f'{OUTPUT_DIR}/{pro}/STRs.csv', 'w') as outfile:
		outfile.write('\n'.join(genes))
		outfile.write('\n')

	genes=sorted(list(strs[(strs.Sample==pro) & (strs.LOEUF<=0.35)]['Gene_id_'].unique()))
	with open(f'{OUTPUT_DIR}/{pro}/STRs (LF).csv', 'w') as outfile:
		outfile.write('\n'.join(genes))
		outfile.write('\n')

# Concat all files into single one
vars=['All coding SNVs', 'Missense', 'LOF', 'Splice', 'Genes del', 'Genes dup', 'STRs', 'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del (LF)', 'Genes dup (LF)', 'STRs (LF)']
for pro in probands:
	df=pd.DataFrame()
	for var in vars:
		vdf=pd.read_csv(f'{OUTPUT_DIR}/{pro}/'+var+'.csv', header=None, names=['Gene'])
		df=pd.concat([df, vdf])
	df.drop_duplicates(inplace=True)

	# Save
	df.to_csv(f'{OUTPUT_DIR}/{pro}/All_variants.csv', index=None, header=None)