import pandas as pd

# Input and output files
FAMILIES="/path/to/Searchlight/families.xlsx"

SNVS="/path/to/SNV/variants.csv" # Use the output of script 1_Variant calling/Searchlight/1_SNV_annotation/7_loeuf_scores.py
CNVS="/path/to/CNV/calls.txt" # Use the output from script 1_Variant calling/Searchlight/2_CNV_annotation/8_annotate_loeuf.py
PRS="/path/to/PRS/scores.csv" # Use the output of script 1_Variant calling/Searchlight/3_PRS/6_merge_scores.py

# The following phenotype files can all be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/
DIAGNOSIS_SUMM="/path/to/Simons_Searchlight_Phase1_16p11.2_Dataset_v11.0/diagnosis_summary.csv"
CBCL_2_5="/path/to/Simons_Searchlight_Phase1_16p11.2_Dataset_v11.0/cbcl_2_5.csv"
CBCL_6_18="/path/to/Simons_Searchlight_Phase1_16p11.2_Dataset_v11.0/cbcl_6_18.csv"
ABCL_18_59="/path/to/Simons_Searchlight_Phase1_16p11.2_Dataset_v11.0/cbcl_2_5.csv"
SRS_PARENT="/path/to/Simons_Searchlight_Phase1_16p11.2_Dataset_v11.0/srs_parent.csv"
SRS_ADULT="/path/to/Simons_Searchlight_Phase1_16p11.2_Dataset_v11.0/srs_adult.csv"
BSI="/path/to/Simons_Searchlight_Phase1_16p11.2_Dataset_v11.0/bsi.csv"
BMI_HC="/path/to/Simons_Searchlight_Phase1_16p11.2_Dataset_v11.0/htwhc.csv"

OUTPUT_PATH="/path/to/output/data/tables/" # There will be two output tables, one each for 16p11.2 deletion and duplication cases. Output tables will be named for the CNV and output to this folder

# Create tables for Searchlight
df=pd.read_excel(FAMILIES)
df=df[df.relationship_to_iip=='Initially identified proband']
df['Sample']=df.sfari_id

cohorts=['16p11.2 deletion', '16p11.2 duplication']
family_types=['16p-deletion', '16p-duplication']
for i in range(2):
	subdf=df[df.family_type==family_types[i]]
	
	# SNVs
	snvs=pd.read_csv(SNVS)
	snvs=snvs[snvs.Sample.isin(subdf.Sample.to_list())]
	
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
	
	subdf=pd.merge(subdf, snv_burden, on='Sample', how='left')
	
	# CNVs
	cnvs=pd.read_csv(CNVS)
	cnvs['low_LOEUF']=cnvs.LOEUF<=0.35
	
	cnv_burden=cnvs[['Sample', 'CNV_Type', 'low_LOEUF', 'Gene_id_']].groupby(['Sample', 'CNV_Type', 'low_LOEUF']).count()
	cnv_burden.reset_index(inplace=True)
	cnv_burden['column']=''
	cnv_burden.loc[(cnv_burden.CNV_Type=='Del') & (~cnv_burden.low_LOEUF), 'column']='Genes del.'
	cnv_burden.loc[(cnv_burden.CNV_Type=='Del') & (cnv_burden.low_LOEUF), 'column']='Genes del. (LF)'
	cnv_burden.loc[(cnv_burden.CNV_Type=='Dup') & (~cnv_burden.low_LOEUF), 'column']='Genes dup.'
	cnv_burden.loc[(cnv_burden.CNV_Type=='Dup') & (cnv_burden.low_LOEUF), 'column']='Genes dup. (LF)'

	cnv_burden=pd.pivot(cnv_burden, index='Sample', columns='column', values='Gene_id_')
	cnv_burden.fillna(0, inplace=True)
	# Add the LF variant counts into the non-LF variant counts
	for vtype in ['Genes del.', 'Genes dup.']:
		cnv_burden[vtype]=cnv_burden[vtype]+cnv_burden[vtype+' (LF)']
	
	cnv_burden=cnv_burden.astype(int)
	
	subdf=pd.merge(subdf, cnv_burden, on='Sample', how='left')
	
	# PRS
	prs=pd.read_csv(PRS)
	prs=prs[['IID', 'schizophrenia_PRS', 'intelligence_PRS', 'educational_attainment_PRS', 'autism_PRS']]
	prs.columns=['Sample', 'SCZ PRS', 'Intelligence PRS', 'Education PRS', 'Autism PRS']
	subdf=pd.merge(subdf, prs, on='Sample', how='left')
	
	# Phenotypes
	# FSIQ
	diagsum=pd.read_csv(DIAGNOSIS_SUMM)
	diagsum['Sample']=diagsum.individual.str.replace('-', '.')
	diagsum['Full scale IQ']=diagsum.best_full_scale_iq
	subdf=pd.merge(subdf, diagsum[['Sample', 'Full scale IQ']], on='Sample', how='left')
	
	# ABCL/CBCL
	cbcl_2_5_df = pd.read_csv(CBCL_2_5)
	cbcl_2_5_df=cbcl_2_5_df[['individual', 'cbcl_2_5.externalizing_problems_t_score', 'cbcl_2_5.internalizing_problems_t_score']]
	cbcl_2_5_df.columns=['individual', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)']
	
	cbcl_6_18_df = pd.read_csv(CBCL_6_18)
	cbcl_6_18_df=cbcl_6_18_df[['individual', 'cbcl_6_18.externalizing_problems_t_score', 'cbcl_6_18.internalizing_problems_t_score']]
	cbcl_6_18_df.columns=['individual', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)']
	
	abcl_18_59_df = pd.read_csv(ABCL_18_59)
	abcl_18_59_df=abcl_18_59_df[['individual', 'abcl_18_59.externalizing_t_score', 'abcl_18_59.internalizing_t_score']]
	abcl_18_59_df.columns=['individual', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)']
	
	cbcl=pd.concat([cbcl_2_5_df, cbcl_6_18_df])
	cbcl=pd.concat([cbcl, abcl_18_59_df])
	cbcl['Sample']=cbcl.individual.str.replace('-', '.')

	# If samples have duplicate scores, take the most recent
	cbcl.drop_duplicates(subset='Sample', keep='last', inplace=True)
	subdf=pd.merge(subdf, cbcl[['Sample', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)']], on='Sample', how='left')
	
	# SRS
	srs_par=pd.read_csv(pheno_dir+'srs_parent.csv')
	srs_par=srs_par[['individual', 'srs_parent.total']]
	srs_par.columns=['individual', 'Social responsiveness (SRS)']
	
	srs_adult=pd.read_csv(pheno_dir+'srs_adult.csv')
	srs_adult=srs_adult[['individual', 'srs_adult.total']]
	srs_adult.columns=['individual', 'Social responsiveness (SRS)']
	
	srs=pd.concat([srs_par, srs_adult])
	srs['Sample']=srs.individual.str.replace('-', '.')
	
	subdf=pd.merge(subdf, srs[['Sample', 'Social responsiveness (SRS)']], on='Sample', how='left')
	
	# BSI
	bsi=pd.read_csv(BSI)
	bsi['Sample']=bsi.individual.str.replace('-', '.')
	bsi['Autism behavior (BSI)']=bsi['bsi.sum_ever']
	
	subdf=pd.merge(subdf, bsi[['Sample', 'Autism behavior (BSI)']], on='Sample', how='left')
	
	# BMI and HC
	hwhc=pd.read_csv(BMI_HC)
	hwhc['Sample']=hwhc.individual.str.replace('-', '.')
	hwhc['BMI z-score']=hwhc['htwhc.bmi_z_score']
	hwhc['Head circumference z-score']=hwhc['htwhc.head_circum_z_score']
	
	subdf=pd.merge(subdf, hwhc[['Sample', 'BMI z-score', 'Head circumference z-score']], on='Sample', how='left')
	
	# Clean up annotations
	cohort_map={'16p-duplication':'16p11.2 duplication', '16p-deletion':'16p11.2 deletion'}
	subdf['Cohort']=subdf.family_type.map(cohort_map)

	subdf['Sex']=subdf.sex.map({'male':'M', 'female':'F'})
	
	subdf=subdf[['Sample', 'Cohort', 'Sex',
				'All coding SNVs', 'Missense', 'LOF', 'Splice', 'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)',
				'Genes del.', 'Genes dup.', 'Genes del. (LF)', 'Genes dup. (LF)',
				'SCZ PRS', 'Intelligence PRS', 'Education PRS', 'Autism PRS',
				'Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)',
				'Social responsiveness (SRS)', 'Autism behavior (BSI)', 'BMI z-score', 'Head circumference z-score']]
	
	# Save
	subdf.to_csv(f'{OUTPUT_PATH}/{cohorts[i]}.csv', index=False)
