import pandas as pd
import numpy as np

# Input and output files
SAMPLE_DF="/path/to/table/of/sample/ids.csv"
ESTONIA_DF="/path/to/table/of/Estonian/sample/ids.csv"
# These files contain sample IDs and family, relationship, sex, and 16p12.1 deletion carrier status information

SNVS="/path/to/SNV/variants.csv" # Use the output of script 1_Variant calling/DD cohort/2_SNV_annotation/coding_annotations/14_loeuf_scores.py
PROMOTER="/path/to/promoter/SNV/variants.csv" # Use the output of promoter variants from script 1_Variant calling/DD cohort/2_SNV_annotation/noncoding_annotations/promoter_UTR5/8_intracohort_filter.py
UTR5="/path/to/5' UTR/SNV/variants.csv" # Use the output of 5' UTR variants from script 1_Variant calling/DD cohort/2_SNV_annotation/noncoding_annotations/promoter_UTR5/8_intracohort_filter.py
ENHANCER="/path/to/enhancer/SNV/variants.tsv" # Use the output of script 1_Variant calling/DD cohort/2_SNV_annotation/noncoding_annotations/enhancer/8_intracohort_filter.py
CNVS="/path/to/CNV/calls.txt" # Use the output from script 1_Variant calling/DD cohort/3_CNV_calling_annotation/merge_all_cnvs/2_annotate_loeuf.py
STRS="/path/to/STR/calls.csv" # Use the output from script 1_Variant calling/DD cohort/4_STR_calling_annotation/11_annotate_loeuf.py
PRS="/path/to/PRS/scores.csv" # Use the output from script 1_Variant calling/DD cohort/5_PRS/6_merge_scores.py

CHILD_DOMAINS="/path/to/child/domain/scores.xlsx"
DE_VRIES="/path/to/de/vries/scores.xlsx"
BMI_HC="/path/to/BMI/and/Head/Circumference/data.xlsx"
HRS_MAT="/path/to/HRS_MAT/scores.csv"
SRS="/path/to/SRS/scores.csv"
MILESTONES="/path/to/developmental/milestone/data.csv"
ADD_PHENOTYPES="/path/to/additional/family/reported/phenotype/data.xlsx"
ADULT_QUESTIONNAIRE="/path/to/adult/questionnaire/data.csv"

OUTPUT="TableS1A.csv"

# Compile data
df = pd.read_csv(SAMPLE_DF)
df2 = pd.read_csv(ESTONIA_DF)
df2['Estonian']='X'

df=pd.concat([df, df2])

###########
# Genetic #
###########
# SNVs
# Coding
snv=pd.read_csv(SNVS)

for vtype in ['all', 'missense', 'lof', 'splice']:
	label={'all':'All_coding_SNVs', 'missense':'Missense', 'lof':'LOF', 'splice':'Splice'}[vtype]
	vtypes=[vtype]
	if vtype=='all':
		vtypes=['missense', 'lof', 'splice']
	for louef in [False, True]:
		if louef:
			label+='_LF'
			samp_counts=snv[(snv.Mut_type.isin(vtypes)) & (snv.LOEUF<=0.35)].Sample.value_counts().to_dict()
		else:
			samp_counts=snv[(snv.Mut_type.isin(vtypes))].Sample.value_counts().to_dict()
		
		df[label]=df.Sample.map(samp_counts)
		df.loc[(df.Sample.isin(snv.Sample.to_list())) & (df[label].isnull()), label]=0

# Promoter
promoter=pd.read_csv(PROMOTER_UTR5)
samp_counts=promoter.Sample.value_counts().to_dict()
df['Promoter']=df.Sample.map(samp_counts)
df.loc[(df.Sample.isin(snv.Sample.to_list())) & (df.Promoter.isnull()), 'Promoter']=0

# 5' UTR
utr5=pd.read_csv(UTR5, sep='\t')
samp_counts=utr5.Sample.value_counts().to_dict()
df["5' UTR"]=df.Sample.map(samp_counts)
df.loc[(df.Sample.isin(snv.Sample.to_list())) & (df["5' UTR"].isnull()), "5' UTR"]=0

# Enhancer
enhancer=pd.read_csv(ENHANCER, sep='\t')
samp_counts=enhancer.Sample.value_counts().to_dict()
df["Enhancer"]=df.Sample.map(samp_counts)
df.loc[(df.Sample.isin(snv.Sample.to_list())) & (df["Enhancer"].isnull()), 'Enhancer']=0

# CNVs
cnv=pd.read_csv(CNVS, sep='\t')
cnv.loc[cnv.LOEUF=='.', 'LOEUF']=np.nan
cnv.LOEUF=cnv.LOEUF.astype(float)
for vtype in ['DEL', 'DUP']:
	label={'DEL':'Genes_del', 'DUP':'Genes_dup'}[vtype]
	for louef in [False, True]:
		if louef:
			label+='_LF'
			samp_counts=cnv[(cnv.Type==vtype) & (cnv.LOEUF<=0.35)].Sample.value_counts().to_dict()
		else:
			samp_counts=cnv[(cnv.Type==vtype)].Sample.value_counts().to_dict()
		
		df[label]=df.Sample.map(samp_counts)
		df.loc[((df.Sample.isin(snv.Sample.to_list())) | (df.Sample.isin(cnv.Sample.to_list()))) & (df[label].isnull()), label]=0

# STRs
str=pd.read_csv(STRS, index_col=0)
samp_counts=str.Sample.value_counts().to_dict()
df['STRs']=df.Sample.map(samp_counts)
df.loc[(df.Sample.isin(snv.Sample.to_list())) & (df["STRs"].isnull()), 'STRs']=0

str.loc[str.LOEUF=='.', 'LOEUF']=np.nan
samp_counts=str[str.LOEUF<=0.35].Sample.value_counts().to_dict()
df['STRs_LF']=df.Sample.map(samp_counts)
df.loc[(df.Sample.isin(snv.Sample.to_list())) & (df["STRs_LF"].isnull()), 'STRs_LF']=0

# PRS
prs=pd.read_csv(PRS)
df=pd.merge(df, prs, right_on='IID', left_on='Sample', how='left')

##############
# Phenotypic #
##############
# Child domains
child_df = pd.read_excel(CHILD_DOMAINS, header=None)
child_df = child_df.set_index(0)
child_df = child_df.T
child_df=child_df[['Child Code', 'Developmental delay/ motor delay/ speech delay/intellectual disability', 'Behavioral', 'Psychiatric', 'Nervous System', 'Congenital anomalies', 'Growth / \nCraniofacial/Skeletal Abnormalities']]
child_df.columns=['Sample', 'ID_DD', 'Behavioral', 'Psychiatric', 'Nerv_Sys', 'Cong_anom', 'Growth_Skeletal']

df=pd.merge(df, child_df, on='Sample', how='left')

# De Vries
devrie_df = pd.read_excel(DE_VRIES, header=None)
devrie_df = devrie_df.set_index(0)
devrie_df = devrie_df.T
dv_dict=dict(zip(devrie_df['Child Code'].to_list(), devrie_df.Total.to_list()))

df['DeVries_score']=df.Sample.map(dv_dict)

# BMI and HC
pheno_df = pd.read_excel(BMI_HC, header=None)
pheno_df = pheno_df.set_index(0)
pheno_df = pheno_df.T
pheno_df = pheno_df[['Child Code', 'BMI Z score', 'Head Circumference Z score']]
pheno_df.columns=['Sample', 'BMI_Z', 'HC_Z']

df=pd.merge(df, pheno_df, on='Sample', how='left')

# HRS-MAT
hrs=pd.read_csv(HRS_MAT)
hrs_mean=hrs[['Sample', 'iq']].groupby('Sample').agg('mean').iq.to_dict()
hrs.drop_duplicates(subset='Sample', keep='last', inplace=True)
hrs_map=dict(zip(hrs.Sample.to_list(), hrs.iq.to_list()))
df['HRS_MAT']=df.Sample.map(hrs_mean)

# SRS
srs_df = pd.read_csv(SRS)
srs_df.index=srs_df['sample'].to_list()
df['SRS'] = df.Sample.map(srs_df.srs_total_raw.to_dict())

# Add in child developmental milestone data
milestonedf=pd.read_csv(MILESTONES)[['Sample', 'Smile?', 'Laugh?', 'Roll over?', 'Sit up without aid?', 'Crawl?',
										'Pull themselves up to standing?', 'Stand without aid?', 'Take first steps/cruise?',
										'Speak first words?', 'Walk/take steps with no assistance?', 'Speak two-word sentences?',
										'Speak in complete sentences?']]
milestonedf.columns=['Sample', 'Age_smile', 'Age_laugh', 'Age_roll_over', 'Age_sit_up', 'Age_crawl', 'Age_pull_up', 'Age_stand', 'Age_steps', 'Age_words', 'Age_walk', 'Age_short_sentence', 'Age_complete_sentence']
df=pd.merge(df, milestonedf, on='Sample', how='outer')

# Add in additional child phenotype data
giri=pd.read_excel(ADD_PHENOTYPES)
keep_cols=['birth_pregnancy_complications', 'preterm_birth', 'microcephaly', 'macrocephaly', 'seizures', 'heart_defects',
			'hearing_loss', 'strabismus', 'vision_problems', 'feeding_problems', 'obesity', 'id_dd', 'motor_delay', 'speech_delay', 'school_aide',
			'learning_disability', 'language_disorder', 'perv_dev_delay', 'asd', 'ashd', 'ocd', 'schiz', 'bpd', 'depression', 'anxiety', 'sleep_trouble']
giri['Sample']=giri.index
giri=giri[['Sample']+keep_cols]

pheno_map={'birth_pregnancy_complications':'Birth/pregnancy complications', 'preterm_birth':'Preterm birth', 'microcephaly':'Microcephaly', 'macrocephaly':'Macrocephaly', 'facial_birth_defect':'Facial dysmorphology',
			'strabismus':'Strabismus', 'seizures':'Seizures', 'heart_defects':'Heart defects', 'hearing_loss':'Hearing loss', 'vision_problems':'Vision problems', 'feeding_problems':'Feeding problems',
			'obesity':'Obesity', 'id_dd':'ID/DD', 'motor_delay':'Motor delay', 'speech_delay':'Speech delay', 'language_disorder':'Language disorder', 'school_aide':'Aide in school',
			'learning_disability':'Learning disability', 'asd':'ASD', 'ashd':'ADHD', 'ocd':'OCD', 'schiz':'Schizophrenia', 'bpd':'BPD', 'depression':'Depression', 'anxiety':'Anxiety',
			'sleep_trouble':'Sleep trouble', 'perv_dev_delay':'PDD'}
giri.columns=['Sample']+[pheno_map[i] for i in keep_cols]

# Add in adult phenotypes
quest=pd.read_csv(ADULT_QUESTIONNAIRE)
quest=quest[['Sample', 'depressed_mood', 'lost_interest', 'sleep', 'mood',
			'anxious_feeling', 'worried_more', 'drugs_for_anxiety', 'anxiety_interferes_life',
			'drinking_interferes_life', 'unable_stop_drinking', 'drug_addiction',
			'unreal_sounds', 'unreal_visions', 'conspiracy']]
# Group symptoms into broad phenotypes
# Depression
quest['depression']=np.nan
quest.loc[(quest.depressed_mood==0) | (quest.lost_interest==0), 'depression']=0
quest.loc[(quest.depressed_mood==1) | (quest.lost_interest==1), 'depression']=1

# Anxiety
quest['anxiety']=np.nan
quest.loc[(quest.anxious_feeling==0) | (quest.worried_more==0) | (quest.drugs_for_anxiety==0) | (quest.anxiety_interferes_life==0), 'anxiety']=0
quest.loc[(quest.anxious_feeling==1) | (quest.worried_more==1) | (quest.drugs_for_anxiety==1) | (quest.anxiety_interferes_life==1), 'anxiety']=1

# Addiction
quest['addiction']=np.nan
quest.loc[(quest.drinking_interferes_life==0) | (quest.unable_stop_drinking==0) | (quest.drug_addiction==0), 'addiction']=0
quest.loc[(quest.drinking_interferes_life==1) | (quest.unable_stop_drinking==1) | (quest.drug_addiction==1), 'addiction']=1

# Psychosis
quest['psychosis']=np.nan
quest.loc[(quest.unreal_sounds==0) | (quest.unreal_visions==0) | (quest.conspiracy==0), 'psychosis']=0
quest.loc[(quest.unreal_sounds==1) | (quest.unreal_visions==1) | (quest.conspiracy==1), 'psychosis']=1

quest=quest[['Sample', 'depression', 'sleep', 'mood', 'anxiety', 'addiction', 'psychosis']]
quest.columns=['Sample', 'Depression (Questionnaire)', 'Sleep trouble (Questionnaire)', 'Mood lability (Questionnaire)', 'Anxiety (Questionnaire)', 'Addiction (Questionnaire)', 'Psychosis (Questionnaire)']

giri=pd.concat([giri, quest])

# Remove any duplicate samples
giri.drop_duplicates(subset='Sample', inplace=True, keep='first')

df=pd.merge(df, giri, on='Sample', how='left')

# Rename relationships
rel_map={'P':'Proband', 'M':'Mother', 'F':'Father', 'S':'Sibling', 'A':'Aunt', 'U':'Uncle', 'C':'Cousin', 'GF':'Grandfather', 'GM':'Grandmother', 'Child':'Child', 'SF':'Step-Father', 'SM':'Step-Mother', 'R':'Child'}
for rm in list(rel_map.keys()):
	for sfx in ['C', 'NC']:
		rel_map[rm+sfx]=rel_map[rm]
df['Relationship']=df.Relationship.map(rel_map)

# Update 16p12.1 deletion carrier information
df['16p12.1 deletion']=df.Carrier.map({'C':'Carrier', 'NC':'Noncarrier'})

# Add in available data information
df['WGS']=np.nan
df.loc[~df.All_coding_SNVs.isnull(), 'WGS']='X'
df['Microarray']=np.nan
df.loc[~df['Microarray batch'].isnull(), 'Microarray']='X'
df.loc[df['Microarray batch'].isnull(), 'Microarray batch']='.'
df.loc[df['Microarray batch'].str.contains('fail'), 'Microarray']=np.nan
df['PRS']=np.nan
df.loc[~df.schizophrenia_PRS.isnull(), 'PRS']='X'

# Reorder and subset columns
df=df[['Sample', 'Family', 'Relationship', 'Age', 'Sex', 'Estonian', '16p12.1 deletion', 'WGS', 'Microarray', 'PRS',
		'Missense', 'Missense_LF', 'LOF', 'LOF_LF', 'Splice', 'Splice_LF', 'All_coding_SNVs', 'All_coding_SNVs_LF',
		'Genes_del', 'Genes_dup', 'Genes_del_LF', 'Genes_dup_LF', 'STRs', 'STRs_LF', 'Enhancer', "5' UTR", 'Promoter',
		'intelligence_PRS', 'schizophrenia_PRS', 'educational_attainment_PRS', 'autism_PRS',
		'ID_DD', 'Behavioral', 'Psychiatric', 'Nerv_Sys', 'Cong_anom', 'Growth_Skeletal', 'DeVries_score',
		'BMI_Z', 'HC_Z', 'HRS_MAT', 'SRS',
		'Age_smile', 'Age_laugh', 'Age_roll_over', 'Age_sit_up', 'Age_crawl', 'Age_pull_up', 'Age_stand', 'Age_steps', 'Age_words', 'Age_walk', 'Age_short_sentence', 'Age_complete_sentence',
		'Depression (Questionnaire)', 'Anxiety (Questionnaire)', 'Sleep trouble (Questionnaire)', 'Psychosis (Questionnaire)', 'Addiction (Questionnaire)', 'Mood lability (Questionnaire)',
		'Birth/pregnancy complications', 'Preterm birth', 'Microcephaly', 'Macrocephaly', 'Strabismus', 'Depression', 'Anxiety', 'Sleep trouble',
		'Seizures', 'Heart defects', 'Hearing loss', 'Vision problems', 'Feeding problems', 'Obesity', 'ID/DD', 'Motor delay', 'Speech delay', 'Language disorder', 'Aide in school',
		'Learning disability', 'ASD', 'ADHD', 'OCD', 'Schizophrenia', 'BPD', 'PDD']]
	
df.columns=['Sample', 'Family', 'Relationship', 'Age (years)', 'Sex', 'Estonian Biobank Sample', '16p12.1 deletion', 'WGS', 'Microarray', 'PRS',
			'Missense', 'Missense (LF)', 'LOF', 'LOF (LF)', 'Splice', 'Splice (LF)', 'All coding SNVs', 'All coding SNVs (LF)',
			'Genes del.', 'Genes dup.', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs', 'STRs (LF)', 'Enhancer', "5' UTR", 'Promoter',
			'Intelligence PRS', 'SCZ PRS', 'Education PRS', 'Autism PRS',
			'ID/DD (Child domain)', 'Behavioral features (Child domain)', 'Psychiatric features (Child domain)', 'Nervous System Abnormalities (Child domain)', 'Congenital Anomalies (Child domain)', 'Growth/Skeletal Defects (Child domain)',
			'De Vries Score',
			'BMI Z Score', 'Head Circumference Z Score', 'HRS-MAT', 'SRS Raw Score',
			'Age first smiled (months)', 'Age first laughed (months)', 'Age first rolled over (months)', 'Age first sat up without assistance (months)', 'Age first crawled (months)', 'Age able to pull self up to stand (months)',
			'Age able to stand alone (months)', 'Age took first steps (months)', 'Age spoke first words (months)', 'Age first walked alone (months)', 'Age first spoke 2-3 word sentences (months)', 'Age spoke complete sentences (months)',
			'Depression (Questionnaire)', 'Anxiety (Questionnaire)', 'Sleep trouble (Questionnaire)', 'Psychosis (Questionnaire)', 'Addiction (Questionnaire)', 'Mood lability (Questionnaire)',
			'Birth/pregnancy complications', 'Preterm birth', 'Microcephaly', 'Macrocephaly', 'Strabismus', 'Depression', 'Anxiety', 'Sleep trouble',
			'Seizures', 'Heart defects', 'Hearing loss', 'Vision problems', 'Feeding problems', 'Obesity', 'ID/DD', 'Motor delay', 'Speech delay', 'Language disorder', 'Aide in school',
			'Learning disability', 'ASD', 'ADHD', 'OCD', 'Schizophrenia', 'BPD', 'PDD']
			
df=df[['Sample', 'Family', 'Relationship', 'Age (years)', 'Sex', 'Estonian Biobank Sample', '16p12.1 deletion', 'WGS', 'Microarray', 'PRS',
		'All coding SNVs', 'Missense', 'LOF', 'Splice', 'Enhancer', "5' UTR", 'Promoter', 'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)',
		'Genes del.', 'Genes dup.', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs', 'STRs (LF)',
		'Intelligence PRS', 'SCZ PRS', 'Education PRS', 'Autism PRS',
		'ID/DD (Child domain)', 'Behavioral features (Child domain)', 'Psychiatric features (Child domain)', 'Nervous System Abnormalities (Child domain)', 'Congenital Anomalies (Child domain)', 'Growth/Skeletal Defects (Child domain)',
		'De Vries Score',
		'BMI Z Score', 'Head Circumference Z Score', 'HRS-MAT', 'SRS Raw Score',
		'Age first smiled (months)', 'Age first laughed (months)', 'Age first rolled over (months)', 'Age first sat up without assistance (months)', 'Age first crawled (months)', 'Age able to pull self up to stand (months)',
		'Age able to stand alone (months)', 'Age took first steps (months)', 'Age spoke first words (months)', 'Age first walked alone (months)', 'Age first spoke 2-3 word sentences (months)', 'Age spoke complete sentences (months)',
		'Depression (Questionnaire)', 'Anxiety (Questionnaire)', 'Sleep trouble (Questionnaire)', 'Psychosis (Questionnaire)', 'Addiction (Questionnaire)', 'Mood lability (Questionnaire)',
		'Birth/pregnancy complications', 'Preterm birth', 'Microcephaly', 'Macrocephaly', 'Strabismus', 'Depression', 'Anxiety', 'Sleep trouble',
		'Seizures', 'Heart defects', 'Hearing loss', 'Vision problems', 'Feeding problems', 'Obesity', 'ID/DD', 'Motor delay', 'Speech delay', 'Language disorder', 'Aide in school',
		'Learning disability', 'ASD', 'ADHD', 'OCD', 'Schizophrenia', 'BPD', 'PDD']]

# Sort
df.Relationship=pd.Categorical(df.Relationship, ['Proband', 'Sibling', 'Mother', 'Father', 'Aunt', 'Uncle', 'Cousin', 'Grandfather', 'Grandmother', 'Child', 'Step-Father', 'Step-Mother'])
df.sort_values(by=['Family', 'Relationship', 'Sample'], inplace=True)

# Save
df.to_csv(OUTPUT, index=False)