import pandas as pd
import numpy as np

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Input and output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
PATHOGENIC_SNVS="/output/final/pathogenic/SNVs.csv" # Use the output of script 2_update_pathogenic_SNVs.py
SNVS="/path/to/SNV/variants.csv" # Use the output of script 1_Variant calling/DD cohort/2_SNV_annotation/coding_annotations/14_loeuf_scores.py
CNVS="/path/to/CNV/calls.txt" # Use the output from script 1_Variant calling/DD cohort/3_CNV_calling_annotation/merge_all_cnvs/2_annotate_loeuf.py

OUTPUT_FIGURE="/path/to/output.pdf" # This will be the figure presented in Fig. S2B
OUTPUT_STATS="/path/to/output/statistics.csv" # This is the data presented in Table S2C

# Compare the burden in probands with additional pathogenic variants and probands without
df=pd.read_csv(TABS1A)
dddf=pd.read_csv(PATHOGENIC_SNVS)

dddf['DualDiagnosis']=~dddf['None']

df=pd.merge(df, dddf[['Sample', 'ClinVar', 'Tier S SFARI', 'Tier 1 DBD', 'Pathogenic CNV', 'DualDiagnosis']], on='Sample', how='inner')

# Remove the second "first hit" variant from dual diagnosis probands
snv=pd.read_csv(SNVS)
for vtype in ['all', 'missense', 'lof', 'splice']:
	label={'all':'All coding SNVs', 'missense':'Missense', 'lof':'LOF', 'splice':'Splice'}[vtype]
	vtypes=[vtype]
	if vtype=='all':
		vtypes=['missense', 'lof', 'splice']
	for louef in [False, True]:
		if louef:
			label+=' (LF)'
			samp_counts=snv[(snv.Mut_type.isin(vtypes)) & (snv.LOEUF<=0.35)].Sample.value_counts().to_dict()
		else:
			samp_counts=snv[(snv.Mut_type.isin(vtypes))].Sample.value_counts().to_dict()
		
		df['tmp']=df.Sample.map(samp_counts)
		df.loc[df.tmp.isnull(), 'tmp']=0
		df[label]=df[label]-df.tmp

cnv=pd.read_csv(CNVS, sep='\t')
cnv=cnv[(cnv.Sample.isin(df.Sample.to_list())) & (cnv.NEJM!='.')]
cnv.loc[cnv.LOEUF=='.', 'LOEUF']=np.nan
cnv.LOEUF=cnv.LOEUF.astype(float)
for vtype in ['DEL', 'DUP']:
	label={'DEL':'Genes del.', 'DUP':'Genes dup.'}[vtype]
	for louef in [False, True]:
		if louef:
			label+=' (LF)'
			samp_counts=cnv[(cnv.Type==vtype) & (cnv.LOEUF<=0.35)].Sample.value_counts().to_dict()
		else:
			samp_counts=cnv[(cnv.Type==vtype)].Sample.value_counts().to_dict()
		
		df['tmp']=df.Samplee.map(samp_counts)
		df.loc[df.tmp.isnull(), 'tmp']=0
		df[label]=df[label]-df.tmp

vars=['All coding SNVs', 'Missense', 'LOF', 'Splice', 'Enhancer', "5' UTR", 'Promoter', 'Genes del.', 'Genes dup.', 'STRs',
		'Intelligence PRS', 'SCZ PRS', 'Education PRS', 'Autism PRS',
		'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)']

def cohens_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)

# Run t-tests
stat_lst=[]
for v in vars:
	subdf=df[~df[v].isnull()]
	
	x=subdf[subdf.DualDiagnosis][v].to_numpy()
	y=subdf[~subdf.DualDiagnosis][v].to_numpy()
	
	res=stats.ttest_ind(x, y, alternative='two-sided')
	
	# Save stats
	stat_lst.append([v, 'Two tailed t test', np.mean(x), np.std(x), len(x), np.mean(y), np.std(y), len(y), cohens_d(x, y), res.statistic, res.pvalue])

# Save statistics
statdf=pd.DataFrame(stat_lst, columns=['Variant', 'Test', 'DDP_Mean', 'DDP_SD', 'DDP_N', 'NDDP_Mean', 'NDDP_SD', 'NDDP_N', "Cohen's D", 'statistic', 'p value'])
statdf.to_csv(OUTPUT_STATS, index=False)

# Make heatmap of results
plotdf=statdf.pivot(index='Test', columns='Variant', values="Cohen's D")
statdf['star']=''
statdf.loc[statdf['p value']<=0.05, 'star']='*'
stardf=statdf.pivot(index='Test', columns='Variant', values="star")

plotdf=plotdf[vars]
stardf=stardf[vars]

biggest=max(abs(statdf["Cohen's D"].to_numpy()))

colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plotdf, cmap=cmap, vmin=-biggest, vmax=biggest, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf)
plt.tight_layout()
plt.savefig(OUTPUT_FIGURE)
plt.close()
