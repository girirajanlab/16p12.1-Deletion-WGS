import pandas as pd
import numpy as np

# Input and output files
COHORTS="/path/to/cohort/file.csv" # Use the cohort output from script 1_identify_first_hits.py
SNVS="/path/to/SNV/variants.csv" # Use the SNV output of script 1_identify_first_hits.py
CNVS="/path/to/CNV/variants.csv" # Use the CNV output of script 1_identify_first_hits.py
STRS="/path/to/STR/variants.csv" # Use the output of script 1_Variant calling/SSC/3_STR_annotation/10_annotate_loeuf.py
PRS="/path/to/PRS/scores.csv" # Use the output of script 1_Variant calling/SSC/4_PRS/6_merge_scores.py

SANDERS_TABLE_S1="/path/to/Sanders/Table_S1.xlsx"
# This is Table S1 from Sanders et al. Neuron 2015 (PubMed: https://pubmed.ncbi.nlm.nih.gov/26402605/) and can be downloaded here: https://ars.els-cdn.com/content/image/1-s2.0-S0896627315007734-mmc2.xlsx

STR_SAMPLES="/path/to/STR/sample/lists.tsv" # Use the output of script 1_Variant calling/SSC/3_STR_annotation/11_identify_samples.py

# The following phenotype files can all be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
CORE_DESCRIPTIVE="/path/to/SSC/proband/ssc_core_descriptive.csv" 
CBCL_2_5="/path/to/SSC/proband/cbcl_2_5.csv"
CBCL_6_18="/path/to/SSC/proband/cbcl_6_18.csv"
SRS="/path/to/SSC/proband/srs_parent.csv"
RBSR="/path/to/SSC/proband/rbs_r.csv"
DCDQ="/path/to/SSC/proband/dcdq.csv"
BMI="/path/to/SSC/proband/ssc_hwhc.csv"

OUTPUT_PATH="/path/to/output/data/tables/" # There will be four output tables, one for each SSC cohort with a distinct group of primary variants. All output tables will be named for the primary variant group and output to this folder

# Gather burden data for each cohort
# For a given cohort, remove primary variants
cohorts=['DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'No primary variant']
cohort=pd.read_csv(COHORTS)

for co in cohorts:
	df=cohort[cohort[co]][['Sample', co]].copy()
	df['Cohort']=co
	# SNV
	snvs=pd.read_csv(SNVS)
	if co=='DBD Tier 1 SNVs':
		snvs=snvs[~snvs.primary]
	snvs=snvs[snvs.Sample.isin(df.Sample.to_list())]
	
	snvs['low_LOEUF']=snvs.LOEUF<=0.35
	
	snv_burden=snvs[['Sample', 'Mut_type', 'low_LOEUF', 'Chrom']].groupby(['Sample', 'Mut_type', 'low_LOEUF']).count()
	snv_burden.reset_index(inplace=True)
	snv_burden['column']=''
	snv_burden.loc[(snv_burden.Mut_type=='missense') & (~snv_burden.low_LOEUF), 'column']='Missense'
	snv_burden.loc[(snv_burden.Mut_type=='missense') & (snv_burden.low_LOEUF), 'column']='Missense (LF)'
	snv_burden.loc[(snv_burden.Mut_type=='lof') & (~snv_burden.low_LOEUF), 'column']='LOF'
	snv_burden.loc[(snv_burden.Mut_type=='lof') & (snv_burden.low_LOEUF), 'column']='LOF (LF)'
	snv_burden.loc[(snv_burden.Mut_type=='splice') & (~snv_burden.low_LOEUF), 'column']='Splice'
	snv_burden.loc[(snv_burden.Mut_type=='splice') & (snv_burden.low_LOEUF), 'column']='Splice (LF)'
	
	snv_burden=pd.pivot(snv_burden, index='Sample', columns='column', values='Chrom')
	snv_burden.fillna(0, inplace=True)
	
	# Add the LF variant counts into the non-LF variant counts
	for vtype in ['Missense', 'LOF', 'Splice']:
		snv_burden[vtype]=snv_burden[vtype]+snv_burden[vtype+' (LF)']
	
	snv_burden['All coding SNVs']=snv_burden.Missense+snv_burden.LOF+snv_burden.Splice
	snv_burden['All coding SNVs (LF)']=snv_burden['Missense (LF)']+snv_burden['LOF (LF)']+snv_burden['Splice (LF)']
	
	snv_burden=snv_burden.astype(int)
	
	df=pd.merge(df, snv_burden, on='Sample', how='left')
	
	# CNV
	cnvs=pd.read_csv(CNVS)
	if co=='Large rare deletions':
		cnvs=cnvs[(~cnvs.primary) | (cnvs['Del/Dup']!='Del')]
	if co=='Large rare duplications':
		cnvs=cnvs[(~cnvs.primary) | (cnvs['Del/Dup']!='Dup')]
	cnvs=cnvs[cnvs.patientID.isin(df.Sample.to_list())]
	
	cnvs['low_LOEUF']=cnvs.LOEUF<=0.35
	cnvs=cnvs[cnvs['Del/Dup'].isin(['Del', 'Dup'])]
	
	cnv_burden=cnvs[['patientID', 'Del/Dup', 'low_LOEUF', 'Gene_id_']].groupby(['patientID', 'Del/Dup', 'low_LOEUF']).count()
	cnv_burden.reset_index(inplace=True)
	cnv_burden['column']=''
	cnv_burden.loc[(cnv_burden['Del/Dup']=='Del') & (~cnv_burden.low_LOEUF), 'column']='Genes del.'
	cnv_burden.loc[(cnv_burden['Del/Dup']=='Del') & (cnv_burden.low_LOEUF), 'column']='Genes del. (LF)'
	cnv_burden.loc[(cnv_burden['Del/Dup']=='Dup') & (~cnv_burden.low_LOEUF), 'column']='Genes dup.'
	cnv_burden.loc[(cnv_burden['Del/Dup']=='Dup') & (cnv_burden.low_LOEUF), 'column']='Genes dup. (LF)'

	cnv_burden=pd.pivot(cnv_burden, index='patientID', columns='column', values='Gene_id_')
	cnv_burden.fillna(0, inplace=True)
	# Add the LF variant counts into the non-LF variant counts
	for vtype in ['Genes del.', 'Genes dup.']:
		cnv_burden[vtype]=cnv_burden[vtype]+cnv_burden[vtype+' (LF)']
	
	cnv_burden=cnv_burden.astype(int)
	
	df=pd.merge(df, cnv_burden, left_on='Sample', right_on='patientID', how='left')
	
	# If sample has no CNVs, check if microarray was run
	# If so, make their CNVs 0, not NA
	micro_cnv_df = pd.read_excel(SANDERS_TABLE_S1)
	micro_pros=micro_cnv_df[micro_cnv_df['CNV_ThisManuscriptSsc']!='0']['Proband'].to_list()
	for col in ['Genes del.', 'Genes dup.', 'Genes del. (LF)', 'Genes dup. (LF)']:
		df.loc[(df[col].isnull()) & (df.Sample.isin(micro_pros)), col]=0
	
	# STRs
	str=pd.read_csv(STRS)
	str=str[str.Sample.isin(df.Sample.to_list())]
	df['STRs']=df.Sample.map(str.Sample.value_counts().to_dict())
	df['STRs (LF)']=df.Sample.map(str[str.LOEUF<0.35].Sample.value_counts().to_dict())
	
	# If sample has SNV data, they should also have STR
	# Replace any NA STR with 0
	df.loc[(~df.Missense.isnull()) & (df['STRs'].isnull()), 'STRs']=0
	df.loc[(~df.Missense.isnull()) & (df['STRs (LF)'].isnull()), 'STRs (LF)']=0

	# Not all samples were able to be assessed for STRs on every chromosome
	# Use the samples_with_strs.txt to restrict STR analysis to only samples with STRs assessed on every chromosome
	str_samps=pd.read_csv(STR_SAMPLES, sep='\t', header=None, names=['SEQID', 'Sample'])
	df.loc[~df.Sample.isin(str_samps.Sample.to_list()), 'STRs']=np.nan
	df.loc[~df.Sample.isin(str_samps.Sample.to_list()), 'STRs (LF)']=np.nan
	
	# PRS
	prs=pd.read_csv(PRS)
	prs=prs[['IID', 'schizophrenia_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'autism_PRS']]
	prs.columns=['Sample', 'SCZ PRS', 'Intelligence PRS', 'Education PRS', 'Autism PRS']
	df=pd.merge(df, prs, on='Sample', how='left')
	
	# Sex, Age, and FSIQ
	core_df = pd.read_csv(CORE_DESCRIPTIVE)
	core_df=core_df[['individual', 'sex', 'age_at_ados', 'ssc_diagnosis_full_scale_iq']]
	core_df.columns=['Sample', 'Sex', 'Age', 'Full scale IQ']
	df=pd.merge(df, core_df, on='Sample', how='left')
	df['Sex']=df['Sex'].map({'male':'M', 'female':'F'})

	# ABCL/CBCL
	cbcl_2_5_df = pd.read_csv(CBCL_2_5)
	cbcl_2_5_df=cbcl_2_5_df[['individual', 'externalizing_problems_t_score', 'internalizing_problems_t_score']]
	cbcl_6_18_df = pd.read_csv(CBCL_6_18)
	cbcl_6_18_df=cbcl_6_18_df[['individual', 'externalizing_problems_t_score', 'internalizing_problems_t_score']]
	cbcl=pd.concat([cbcl_2_5_df, cbcl_6_18_df])
	cbcl.columns=['Sample', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)']

	# If samples have duplicate scores, take the 6-18 score
	cbcl.drop_duplicates(subset='Sample', keep='last', inplace=True)
	df=pd.merge(df, cbcl, on='Sample', how='left')

	# SRS
	srs_df = pd.read_csv(SRS)
	srs_df=srs_df[['individual', 'total']]
	srs_df.columns=['Sample', 'Social responsiveness (SRS)']
	df=pd.merge(df, srs_df, on='Sample', how='left')

	# RBS_R
	rbsr_df = pd.read_csv(RBSR)
	rbsr_df=rbsr_df[['individual', 'overall_score']]
	rbsr_df.columns=['Sample', 'Repetitive behavior (RBS-R)']
	df=pd.merge(df, rbsr_df, on='Sample', how='left')

	# DCDQ
	dcdq_df = pd.read_csv(DCDQ)
	dcdq_df=dcdq_df[['individual', 'total']]
	dcdq_df.columns=['Sample', 'Coordination disorder (DCDQ)']
	df=pd.merge(df, dcdq_df, on='Sample', how='left')

	# BMI
	hwhc_df = pd.read_csv(BMI)
	hwhc_df=hwhc_df[['individual', 'bmi_z_score']]
	hwhc_df.columns=['Sample', 'BMI z-score']
	df=pd.merge(df, hwhc_df, on='Sample', how='left')
	
	# Organize columns and save
	df=df[['Sample', 'Cohort', 'Sex', 'Age',
			'All coding SNVs', 'Missense', 'LOF', 'Splice', 'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)',
			'Genes del.', 'Genes dup.', 'Genes del. (LF)', 'Genes dup. (LF)',
			'STRs', 'STRs (LF)',
			'SCZ PRS', 'Intelligence PRS', 'Education PRS', 'Autism PRS',
			'Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Social responsiveness (SRS)',
			'Repetitive behavior (RBS-R)', 'Coordination disorder (DCDQ)', 'BMI z-score']]
	
	df.to_csv(f'{OUTPUT_PATH}/{co}.csv', index=False)
	