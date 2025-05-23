import pandas as pd
import numpy as np

# Condense MyCode data into a single table

# Input and output tables
SAMPLE_TABLE="/path/to/MyCode/sample/information.csv" # This file contains sample information, including sample ID, family ID,sex, and available data types

SNVS="/path/to/SNV/variants.csv" # Use the output of script 1_Variant calling/MyCode/1_SNV_annotation/13_loeuf_scores.py
CNVS="/path/to/CNV/variants.csv" # Use the output of script 1_Variant calling/MyCode/2_CNV_calling_annotation/9_annotate_loeuf.py
PRS="/path/to/PRS/scores.csv" # Use the output of script 1_Variant calling/MyCode/3_PRS/7_merge_scores.py

ICD10_CHAPTERS="/path/to/MyCode/ICD10/Chapters.csv" # Use the chapter output of script 1_parse_ICD10.py
ICD10_PHENOTYPES="/path/to/MyCode/ICD10/phenotypes.csv" # Use the phenotype output of script 1_parse_ICD10.py

OUTPUT_TABLE="/path/to/output/MyCode/data/table.csv"

# Load samples
df=pd.read_csv(SAMPLE_TABLE)[['SEQN_ID', 'FID', 'PT_SEX', 'gVCF', 'ICD10']]

# Convert sex annotations
df['Sex']=df.Sex.map({'Male':'M', 'Female':'F'})

# SNVs
# All coding SNVs = LOF + missense (CADD>=25) + splice (CADD>=25)
snvs=pd.read_csv(SNVS)
snvs=snvs[snvs.Sample.isin(df.Sample.to_list())]
snvs=snvs[(snvs.Mut_type=='lof') | (snvs.CADD_PHRED>=25)]

df['All_coding_SNVs']=df.Sample.map(snvs.Sample.value_counts().to_dict())
df['All_coding_SNVs_LF']=df.Sample.map(snvs[snvs.LOEUF<=0.35].Sample.value_counts().to_dict())

df['Missense']=df.Sample.map(snvs[snvs.Mut_type=='missense'].Sample.value_counts().to_dict())
df['Missense_LF']=df.Sample.map(snvs[(snvs.Mut_type=='missense') & (snvs.LOEUF<=0.35)].Sample.value_counts().to_dict())

df['LOF']=df.Sample.map(snvs[snvs.Mut_type=='lof'].Sample.value_counts().to_dict())
df['LOF_LF']=df.Sample.map(snvs[(snvs.Mut_type=='lof') & (snvs.LOEUF<=0.35)].Sample.value_counts().to_dict())

df['Splice']=df.Sample.map(snvs[snvs.Mut_type=='splice'].Sample.value_counts().to_dict())
df['Splice_LF']=df.Sample.map(snvs[(snvs.Mut_type=='splice') & (snvs.LOEUF<=0.35)].Sample.value_counts().to_dict())

# If sample has SNV data, fill all NA with 0
for col in ['All_coding_SNVs', 'Missense', 'LOF', 'Splice', 'All_coding_SNVs_LF', 'Missense_LF', 'LOF_LF', 'Splice_LF']:
	df.loc[(df.gVCF) & (df[col].isnull()), col]=0

# CNVs
# We need the number of genes deleted and duplicated from microarray data
cnv=pd.read_csv(CNVS, sep='\t')
# Remove 16p12.1 samples and calls
cnv=cnv[cnv.Sample.isin(df.Sample.to_list())]
cnv=cnv[~((cnv.Type=='DEL') & (cnv.Pathogenic=='16p12.1'))]

cnv.replace('.', np.nan, inplace=True)
cnv.LOEUF=cnv.LOEUF.astype(float)

df['Genes_del']=df.Sample.map(cnv[cnv.Type=='DEL']['Sample'].value_counts().to_dict())
df['Genes_dup']=df.Sample.map(cnv[cnv.Type=='DUP']['Sample'].value_counts().to_dict())

df['Genes_del_LF']=df.Sample.map(cnv[(cnv.Type=='DEL') & (cnv.LOEUF<=0.35)]['Sample'].value_counts().to_dict())
df['Genes_dup_LF']=df.Sample.map(cnv[(cnv.Type=='DUP') & (cnv.LOEUF<=0.35)]['Sample'].value_counts().to_dict())

# All samples have CNV data, so fill any NAs with 0
df=df.fillna({'Genes_del':0, 'Genes_dup':0, 'Genes_del_LF':0, 'Genes_dup_LF':0})

# PRS
prs=pd.read_csv(PRS)
prs.columns=['FID', 'Sample', 'autism_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'SCZ_PRS']
df=pd.merge(df, prs[['Sample', 'autism_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'SCZ_PRS']], on='Sample', how='left')

# Clean up column annotations
df=df[['Sample', 'Sex', 'ICD10', 
		'All_coding_SNVs', 'Missense', 'LOF', 'Splice', 'Genes_del', 'Genes_dup',
		'All_coding_SNVs_LF', 'Missense_LF', 'LOF_LF', 'Splice_LF', 'Genes_del_LF', 'Genes_dup_LF',
		'intelligence_PRS', 'SCZ_PRS', 'educational_attainment_PRS', 'autism_PRS',]]

df.columns=['Sample', 'Sex', 'ICD10',
			'All coding SNVs', 'Missense', 'LOF', 'Splice', 'Genes del.', 'Genes dup.',
			'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del. (LF)', 'Genes dup. (LF)',
			'Intelligence PRS', 'SCZ PRS', 'Education PRS', 'Autism PRS']

# Add in ICD10 chapter information
icddf=pd.read_csv(ICD10_CHAPTERS)
icddf=icddf[['Sample',
			'Chapter_20', 'Chapter_30', 'Chapter_40', 'Chapter_50', 'Chapter_60', 'Chapter_70', 'Chapter_80', 'Chapter_90',
			'Chapter_100', 'Chapter_110', 'Chapter_120', 'Chapter_130', 'Chapter_140', 'Chapter_150', 'Chapter_170']]
icddf.columns=['Sample',
				'Neoplasms', 'Blood', 'Endocrine/Metabolic', 'Mental/behavioral disorders','Nervous system',
				'Eye', 'Ear', 'Circulatory system', 'Respiratory system', 'Digestive system', 'Skin/subcutaeous tissue',
				'Musc. system/connective tissue', 'Genitourinary system', 'Pregnancy/childbirth', 'Congenital malformations']

df=pd.merge(df, icddf, on='Sample', how='left')

# Add in comparison phenotype information
uniondf=pd.read_csv(ICD10_PHENOTYPES)
uniondf=uniondf[['Sample',
					'sleep_ICD10', 'addiction_ICD10', 'depression_ICD10', 'anxiety_ICD10', 'psychosis_ICD10']]
uniondf.columns=['Sample', 'Sleep trouble (ICD10)', 'Addiction (ICD10)', 'Depression (ICD10)', 'Anxiety (ICD10)', 'Psychosis (ICD10)']

df=pd.merge(df, uniondf, on='Sample', how='left')

# Save
df.to_csv(OUTPUT_TABLE, index=False)

