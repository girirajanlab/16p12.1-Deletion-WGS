import pandas as pd
import numpy as np

# Condense genetic and phenotypic data for SPARK samples

# Input and output files
SAMPLE_TABLE="/path/to/SPARK/sample/information.csv" # This file contains sample information, including sample ID, family ID, sex, and role in the family

SNVS="/path/to/SNV/burden/table.csv" # Use the output of script 1_Variant calling/SPARK/1_SNV_annotation/8_burden_table.py
CNVS="/path/to/CNV/variants.csv" # Use the output of script 1_Variant calling/SPARK/2_CNV_calling_annotation/9_annotate_loeuf.py
PRS="/path/to/PRS/scores.csv" # Use the output of script 1_Variant calling/SPARK/3_PRS/6_merge_scores.py

PHENOTYPE_PATH="/path/to/SPARK Collection version5/phenotype/files/" # This should be the path to three files: basic_medical_screening.csv, background_history_child.csv, and individuals_registration.csv
# These phenotype files can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/

OUTPUT_TABLE="/path/to/output/SPARK/data/table.csv"

# Load samples
df=pd.read_csv(SAMPLE_TABLE)[['Sample', 'Family', 'Sex', 'role']]

# Restrict to kids
df=df[df.role.isin(['Proband', 'Unaffected Sibling', 'Affected Sibling'])]

# Convert sex annotations
df['Sex']=df.Sex.map({'male':'M', 'female':'F'})

# SNVs
snvs=pd.read_csv(SNVS)
df=pd.merge(df, snvs, on='Sample', how='left')

# CNVs
cnv=pd.read_csv(CNVS)
# Remove 16p12.1 deletion calls
cnv=cnv[cnv.Sample.isin(df.Sample.to_list())]
cnv=cnv[~((cnv.Type=='DEL') & (cnv.Pathogenic=='16p12.1'))]

cnv.replace('.', np.nan, inplace=True)
cnv.LOEUF=cnv.LOEUF.astype(float)

df['Genes_del']=df.Sample.map(cnv[cnv.Type=='DEL'].Sample.value_counts().to_dict())
df['Genes_dup']=df.Sample.map(cnv[cnv.Type=='DUP'].Sample.value_counts().to_dict())

df['Genes_del_LF']=df.Sample.map(cnv[(cnv.Type=='DEL') & (cnv.LOEUF<=0.35)].Sample.value_counts().to_dict())
df['Genes_dup_LF']=df.Sample.map(cnv[(cnv.Type=='DUP') & (cnv.LOEUF<=0.35)].Sample.value_counts().to_dict())

# Fill any NA CNV data with 0
df=df.fillna({'Genes_del':0, 'Genes_dup':0, 'Genes_del_LF':0, 'Genes_dup_LF':0})

# PRS
prs=pd.read_csv(PRS)
prs.columns=['FID', 'Sample', 'autism_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'SCZ_PRS']
df=pd.merge(df, prs[['Sample', 'autism_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'SCZ_PRS']], on='Sample', how='left')

# Clean up columns
df=df[['Sample', 'Family', 'Sex',
                'All_coding_SNVs', 'LOF', 'Missense', 'Splice', 'All_coding_SNVs_LF', 'LOF_LF', 'Missense_LF', 'Splice_LF',
                'Genes_del', 'Genes_dup', 'Genes_del_LF', 'Genes_dup_LF',
                'autism_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'SCZ_PRS']]
df.columns=['Sample', 'Family', 'Sex',
                        'All coding SNVs', 'LOF', 'Missense', 'Splice', 'All coding SNVs (LF)', 'LOF (LF)', 'Missense (LF)', 'Splice (LF)',
                        'Genes del.', 'Genes dup.', 'Genes del. (LF)', 'Genes dup. (LF)',
                        'Autism PRS', 'Intelligence PRS', 'Education PRS', 'SCZ PRS']

# Add phenotype information
pheno_data_dict={'basic_medical_screening':{
						'Birth/pregnancy complications':'med_cond_birth',
						'Preterm birth':'birth_prem',
						'Microcephaly':'growth_microceph',
						'Macrocephaly':'growth_macroceph',
						'Facial dysmorphology':'birth_def_fac',
						'Seizures':'neuro_sz',
						'Heart defects':'birth_def_thorac_heart',
						'Hearing loss':'visaud_deaf',
						'Vision problems':'med_cond_visaud',
						'Feeding problems':'feeding_dx',
						'Obesity':'growth_obes',
						'ID/DD':'dev_id',
						'Motor delay':'dev_motor',
						'Speech delay':'dev_speech',
						'Language disorder':'dev_lang_dis',
						'Learning disability':'dev_ld',
						'ADHD':'behav_adhd',
						'OCD':'mood_ocd',
						'Schizophrenia':'schiz',
						'BPD':'mood_bipol',
						'Depression':'mood_dep',
						'Anxiety':'mood_anx',
						'Sleep trouble':'sleep_dx'},
				'background_history_child.csv':{
						'Aide in school':'sped_aide'}}
for file in (pheno_data_dict.keys()):
	phenos=list(pheno_data_dict[file].keys())
	cols=[pheno_data_dict[file][key] for key in phenos]
	
	filedf=pd.read_csv(f'{PHENOTYPE_PATH}/{file}')
	filedf=filedf[['subject_sp_id']+cols]
	filedf.columns=['Sample']+phenos
	df=pd.merge(df, filedf, on='Sample', how='left')
	
indiv_reg=pd.read_csv(f'{PHENOTYPE_PATH}/individuals_registration.csv')
indiv_reg['PDD']=0
indiv_reg.loc[indiv_reg.diagnosis.str.contains('Pervasive', na=False), 'PDD']=1
indiv_reg['ASD']=1
indiv_reg.loc[indiv_reg.diagnosis==0, 'ASD']=0
indiv_reg['Sample']=indiv_reg.subject_sp_id.to_list()

df=pd.merge(df, indiv_reg[['Sample', 'PDD', 'ASD']], on='Sample', how='left')

# Save
df.to_csv(OUTPUT_TABLE, index=False)


	
