import pandas as pd
import numpy as np

import random

# Identify controls for burden comparisons

# Input and output files
CNVS="/path/to/CNV/variants.csv" # Use the output of script 1_Variant calling/SPARK/2_CNV_calling_annotation/6_frequency_filter.py
CNV_QC="/path/to/CNV/QC/file.txt" # Use the output of script 1_Variant calling/SPARK/2_CNV_calling_annotation/2_check_QC.py

SNVS="/path/to/SNV/burden/table.csv" # Use the output of script 1_Variant calling/SPARK/1_SNV_annotation/8_burden_table.py

MASTERTABLE="/path/to/SPARK/mastertables/SPARK.iWGS_v1.1.mastertable.2023_03.tsv" # This file can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/

OUTPUT_TABLE="/path/to/output/table.csv"

# Load data
df=pd.read_csv(CNVS)

# Get 16p12.1 deletion samples
del_16p=df[(df.Type=='DEL') & (df.NEJM_Name=='16p12.1')]['Sample'].to_list()

# Identify samples with large (>=500kb), rare (<=0.1%) CNVs or CNVs in known pathogenic regions (to exclude from No CNV group)
rm_samps=list(df[((df.Length>=500000) & (df.Control_freq<=0.001)) | (df.NEJM_Name!='.')]['Sample'].unique())

# Check the QC metrics to identify any samples without CNVs
cnv_qc=pd.read_csv(CNV_QC, sep='\t')
cnv_qc['Case_Control']='No CNV Control'
cnv_qc.loc[cnv_qc.Sample.isin(del_16p), 'Case_Control']='16p12.1 deletion'
df=cnv_qc[['Sample', 'Case_Control']].copy()

# Add family information
fam=pd.read_csv(MASTERTABLE, sep='\t')
fam.index=fam.spid.to_list()

df['Family']=df.Sample.map(fam.sfid.to_dict())

# Add sex and age
df['Sex']=df.Sample.map(fam.sex.to_dict()).map({1:'male', 2:'female'})
df['Age_months']=df.Sample.map(fam.age_m.to_dict())

# Remove samples with no age
df=df[~df.Age_months.isnull()]

# Add in SNV and CNV burden
# SNV
snvs=pd.read_csv(snvs)
df=pd.merge(df, snvs, on='Sample', how='left')
df['SNV']=True
df.loc[df.All_coding_SNVs.isnull(), 'SNV']=False

# Restrict samples to those with SNV data
df=df[df.SNV]

# Subset age and sex matched controls
df['age_sex']=df.Age_months.astype(str)+'.'+df.Sex

case_as=pd.DataFrame(df[df.Case_Control=='16p12.1 deletion']['age_sex'].value_counts())
case_as.columns=['case_num']

case_as['nocnv_num']=case_as.index.map(df[df.Case_Control=='No CNV Control']['age_sex'].value_counts().to_dict())
case_as['nocnv_max_vals']=(case_as.nocnv_num/case_as.case_num).astype(int)

case_as['nocnv_controls_needed']=case_as.case_num*min(case_as.nocnv_max_vals.to_numpy())

cases=df[(df.Case_Control=='16p12.1 deletion') & (df.age_sex.isin(case_as.index.to_list()))]

df=df[df.Case_Control!='16p12.1 deletion']
df=df[df.age_sex.isin(cases.age_sex.to_list())]

# Randomly select controls for each case
df['keep']=False
age_sex=list(cases.age_sex.unique())
age_sex.sort()
for ags in age_sex:
	poss_samps=df[df.age_sex==ags].copy()
	
	needed=case_as.loc[ags, 'nocnv_controls_needed']
	random.seed(205)
	chosen=random.sample(poss_samps[poss_samps.Case_Control=='No CNV Control']['Sample'].to_list(), needed)
	df.loc[df.Sample.isin(chosen), 'keep']=True

df=df[df.keep]
df=pd.concat([cases, df])

df=df[['Sample', 'Case_Control',
		'All_coding_SNVs', 'Missense', 'LOF', 'Splice',
		'All_coding_SNVs_LF', 'Missense_LF', 'LOF_LF', 'Splice_LF']]

df.columns=['Sample', 'Case_Control',
			'All coding SNVs', 'Missense', 'LOF', 'Splice',
			'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)']

# Save
df.to_csv(OUTPUT_TABLE, index=False)
