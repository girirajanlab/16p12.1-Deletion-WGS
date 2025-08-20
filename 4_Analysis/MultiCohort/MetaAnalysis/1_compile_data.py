import pandas as pd

# Compile data for meta analyses

# Input and output files
UKB_QUESTIONNAIRE="/path/to/UKB/questionnaire/stats.csv" # Use the output of script 4_Analysis/UKB/Figure 6/4_phenotype_variant_associations.ipynb
DD_ADULT="/path/to/DD/adult/stats.csv" # Use the output of script 4_Analysis/DD cohort/Figure 6/1_adult_associations.py
UKB_EHR="/path/to/UKB/EHR/stats.csv" # Use the output of script 4_Analysis/UKB/Figure 6/5_icd_variant_associations.ipynb
AoU="/path/to/AoU/EHR/stats.csv"
MYCODE="/path/to/MyCode/EHR/stats.csv" # Use the output of 4_Analysis/MyCode/Figure 6/1_phenotype_variant_associations.py
DD_CHILD="/path/to/DD/child/stats.csv" # Use the output of 4_Analysis/DD cohort/Figure 6/2_child_associations.py
SPARK="/path/to/SPARK/stats.csv" # Use the output of 4_Analysis/SPARK/Figure 6/1_phenotype_variant_associations.py

QUESTIONNAIRE_STATS="/path/to/output/DD/UKB/stats.csv"
EHR_STATS="/path/to/output/UKB/MyCode/AoU/stats.csv"
CHILD_STATS="/path/to/output/DD/SPARK/stats.csv"

# Compile adult questionnaire data
ukb=pd.read_csv(UKB_QUESTIONNAIRE)
ukb['Cohort']='UKB'
dd=pd.read_csv(DD_ADULT)
dd['Cohort']='DD (adults)'

quest_meta=pd.concat([ukb, dd])

# Restrict to phenotype-variant combinations assessed in both cohorts
def filter_res(df, threshold=2):
	# Remove any cases where SD=0 (not given weight in meta-analysis
	df=df[~((df['Case SD']==0) | (df['Control SD']==0))]
	
	# Remove any variant_phenotype associations that could not be assessed in all cohorts
	df['vp']=df.Variant+df.Phenotype
	
	vp_count=df.vp.value_counts()
	vp_count=vp_count[vp_count==threshold]
	df=df[df.vp.isin(vp_count.index.to_list())]
	
	return df

quest_meta=filter_res(quest_meta, 2)

# Save
quest_meta.to_csv(QUESTIONNAIRE_STATS, index=False)

# Compile EHR data
ukb=pd.read_csv(UKB_EHR)
ukb['Cohort']='UKB'
mycode=pd.read_csv(MYCODE)
mycode['Cohort']='MyCode'
aou=pd.read_csv(AoU)
aou['Cohort']='AoU'

icd_meta=pd.concat([ukb, mycode])
icd_meta=pd.concat([icd_meta, aou])
icd_meta=filter_res(icd_meta, 3)

# Save
icd_meta.to_csv(EHR_STATS, index=False)

# Compile child data
dd=pd.read_csv(DD_CHILD)
dd['Cohort']='DD (children)'
spark=pd.read_csv(SPARK)
spark['Cohort']='SPARK'

child_meta=pd.concat([dd, spark])
child_meta=filter_res(child_meta, 2)

child_meta.to_csv(CHILD_STATS, index=False)
