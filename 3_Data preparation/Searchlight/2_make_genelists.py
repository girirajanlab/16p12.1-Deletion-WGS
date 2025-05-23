import pandas as pd
import numpy as np

import subprocess

# Make lists of secondary variant genes for 16p11.2 CNV probands

# Input and output files
SNVS="/path/to/SNV/variants.csv" # Use the output of script 1_Variant calling/Searchlight/1_SNV_annotation/7_loeuf_scores.py
CNVS="/path/to/CNV/calls.txt" # Use the output from script 1_Variant calling/Searchlight/2_CNV_annotation/8_annotate_loeuf.py
COHORT_DATA_DIR="/path/to/cohort/tables/" # Use the output directory from script 1_compile_data.py

OUTPUT_DIR="/path/to/output/directory/" # Output files will be saved here
# . - will contain genes across all 16p11.2 CNV probands broken up by primary variant
# by_proband/ - contains two sub-directories (one for each CNV), then a directory for each proband with separate files of secondary variant genes broken up by variant type

cohorts = ['16p11.2 deletion', '16p11.2 duplication']

# Variants by cohort
for cohort in cohorts:
	df=pd.read_csv(f'{COHORT_DATA_DIR}/{cohort}.csv')[['Sample', 'Cohort']]
	df.drop_duplicates(inplace=True)
	
	probands=df.Sample.to_list()
	
	for pro in probands:
		cmd=['mkdir', OUTPUT_DIR+'\\by_proband\\'+cohort+'\\'+pro]
		subprocess.run(cmd, shell=True)

	# SNVs
	snv=pd.read_csv(SNVS)
	snv=snv[snv.Sample.isin(df.Sample.to_list())]
	snv['Gene_id']=snv['Gene_id'].str.split(';')
	snv=snv.explode('Gene_id')
	snv.Gene_id=snv.Gene_id.str.split('.', expand=True)[0]

	snv_genes=list(snv.Gene_id.unique())

	for pro in probands:
		for mt in ['missense', 'lof', 'splice']:
			sub=snv[(snv.Mut_type==mt) & (snv.Sample==pro)]
			genes=sorted(list(sub['Gene_id'].unique()))
			
			# Save genes to list
			name=mt.title()
			if mt=='lof':
				name='LOF'
			with open(f'{OUTPUT_DIR}/by_proband/{cohort}/{pro}/{name}.csv', 'w') as outfile:
				outfile.write('\n'.join(genes))
				outfile.write('\n')

			# LF
			sub=sub[sub.LOEUF<=0.35]
			genes=sorted(list(sub['Gene_id'].unique()))
			with open(f'{OUTPUT_DIR}/by_proband/{cohort}/{pro}/{name} (LF).csv', 'w') as outfile:
				outfile.write('\n'.join(genes))
				outfile.write('\n')

	for pro in probands:
		genes=sorted(list(snv[snv.Sample==pro]['Gene_id'].unique()))
		with open(f'{OUTPUT_DIR}/by_proband/{cohort}/{pro}/All coding SNVs.csv', 'w') as outfile:
			outfile.write('\n'.join(genes))
			outfile.write('\n')

		genes=sorted(list(snv[(snv.Sample==pro) & (snv.LOEUF<=0.35)]['Gene_id'].unique()))
		with open(f'{OUTPUT_DIR}by_proband/{cohort}/{pro}/All coding SNVs (LF).csv', 'w') as outfile:
			outfile.write('\n'.join(genes))
			outfile.write('\n')

	# CNVs
	cnv=pd.read_csv(CNVS)
	cnv=cnv[cnv.Sample.isin(df.Sample.to_list())]
	cnv=cnv[cnv.CNV_Type.isin(['Del', 'Dup'])]
	cnv['Gene_id_']=cnv['Gene_id_'].str.split(';')
	cnv=cnv.explode('Gene_id_')
	cnv_genes=list(cnv.Gene_id_.unique())

	cnv_names={'Del':'Genes del', 'Dup':'Genes dup'}
	for pro in probands:
		for cn in cnv_names.keys():
			sub=cnv[(cnv['CNV_Type']==cn) & (cnv.Sample==pro)]
			genes=sorted(list(sub['Gene_id_'].unique()))
			
			with open(f'{OUTPUT_DIR}/by_proband/{cohort}/{pro}/{cnv_names[cn]}.csv', 'w') as outfile:
				outfile.write('\n'.join(genes))
				outfile.write('\n')

			# LF
			sub=sub[sub.LOEUF<=0.35]
			genes=sorted(list(sub['Gene_id_'].unique()))
			with open(f'{OUTPUT_DIR}/by_proband/{cohort}/{pro}/{cnv_names[cn]} (LF).csv', 'w') as outfile:
				outfile.write('\n'.join(genes))
				outfile.write('\n')

	all_genes=sorted(list(set(snv_genes+cnv_genes)))
	all_genes=[i for i in all_genes if "ENSGR" not in i]

	print(all_genes[0:10])

	with open(f'{OUTPUT_DIR}/{cohort}_genes.csv', 'w') as outfile:
		outfile.write('\n'.join(all_genes))
		outfile.write('\n')
	
	# Concat all variant files into single file per proband
	vars=['All coding SNVs', 'Missense', 'LOF', 'Splice', 'Genes del', 'Genes dup', 'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del (LF)', 'Genes dup (LF)']
	for pro in probands:
		df=pd.DataFrame()
		for var in vars:
			vdf=pd.read_csv(f'{OUTPUT_DIR}/by_proband/{cohort}/{pro}/{var}.csv', header=None, names=['Gene'])
			df=pd.concat([df, vdf])
		df.drop_duplicates(inplace=True)

		# Save
		df.to_csv(f'{OUTPUT_DIR}/by_proband/{cohort}/{pro}/All_variants.csv', index=None, header=None)