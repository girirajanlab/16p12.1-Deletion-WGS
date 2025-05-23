import pandas as pd

import scipy.stats as stats

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Compare the frequencies of phenotypes across cohorts

# Input and output files
INPUT_FILES=['/path/to/output/DD/cohort/prevalence.csv', # Use the output of script 4_Analysis\DD cohort\Figure 6\3_phenotype_prevalence.py
				'/path/to/output/SPARK/cohort/prevalence.csv', # Use the output of script 4_Analysis\SPARK\Figure 6\2_phenotype_prevalence.py
				'/path/to/output/MyCode/cohort/prevalence.csv', # Use the output of script 4_Analysis\MyCode\Figure 6\2_phenotype_prevalence.py
				'/path/to/output/UKB/cohort/prevalence.csv',
				'/path/to/output/AoU/cohort/prevalence.csv']

OUTPUT_STATS="/path/to/output/statistics.csv" # Output statistics presented in Table S5B
OUTPUT_FIG="/path/to/output/figure.pdf" # Output figures presenteed in Figs 6A, S5C-D

# Load data
df=pd.DataFrame()
for f in INPUT_FILES:
	fdf=pd.read_csv(f)
	# For AoU, skip any phenotypes present in less than 20 people
	if fdf.Cohort.to_list()[0]=='All of Us':
		fdf=fdf[fdf.Notes.isnull()]
	df=pd.concat([df, fdf])

df=df[df.Total>0]
df['Phenotype Prevalence (percent)']=100*df.Count/df.Total
df.reset_index(inplace=True, drop=True)

# Perform Fisher's Exact tests between specific groups
phenos=list(df.Phenotype.unique())

stat_lst=[]
for p in phenos:
	for groups in [['DD (adults)', 'UKB (Questionnaire)'], ['MyCode', 'UKB (ICD10)'], ['MyCode', 'All of Us'], ['UKB (ICD10)', 'All of Us'], ['DD (children)', 'SPARK']]:
		subdf=df[(df.Cohort.isin(groups)) & (df.Phenotype==p)]
		if subdf.shape[0]!=2:
			continue
		
		subdf.index=subdf.Cohort.to_list()
		fishdf=subdf[['Count', 'Without']].astype(int)
		
		res=stats.fisher_exact(fishdf)
		or_res=stats.contingency.odds_ratio(fishdf)
		ci=or_res.confidence_interval()
		
		annos={}
		for g in groups:
			g_count=subdf.loc[g, 'Count']
			g_total=subdf.loc[g, 'Total']
			g_perc=100*g_count/g_total
			if g_perc>=10:
				g_anno='%i%% (%i/%i)' % (g_perc, g_count, g_total)
			else:
				g_anno='%.1f%% (%i/%i)' % (g_perc, g_count, g_total)
			annos[g]=g_anno
		
		stat_lst.append([p, groups[0], groups[1], annos[groups[0]], annos[groups[1]], res.statistic, ci[0], ci[1], res.pvalue])
statdf=pd.DataFrame(stat_lst, columns=['Phenotype', 'Cohort 1', 'Cohort 2', 'Cohort 1 prevalence', 'Cohort 2 prevalence', "Fisher's Exact statistic", '95% C.I. lower', '95% C.I. upper', 'p value'])

# Save
statdf.to_csv(OUTPUT_STATS, index=False)

# Plot prevalence of specific phenotypes in 3 batches
# 1. Mixed children adult: ID/DD, ASD, Anxiety, Depression, Psychosis
# 2. Adult only: Sleep trouble, Addiction
# 3. Child only: All child only phenotypes with at least 20% prevalence in one cohort

pdf=PdfPages(OUTPUT_FIG)

palette={'DD (children)':'#254061', 'SPARK':'#b12422', 'blank2':'white',
			'DD (adults)':'#B9CDE5', 'UKB (Questionnaire)':'#b084cc', 'blank':'white', 'MyCode':'#87ba44', 'UKB (ICD10)':'#45235F', 'All of Us':'#E3A857'}

# Mixed phenotypes
fig, ax=plt.subplots(figsize=(8, 4))
pheno_order=['ID/DD', 'ASD', 'Anxiety', 'Depression', 'Psychosis']
cohort_order=['DD (children)', 'SPARK', 'blank2', 'DD (adults)', 'UKB (Questionnaire)', 'blank', 'MyCode', 'UKB (ICD10)', 'All of Us']
sns.barplot(data=df[(df.Phenotype!='Depression') | (~df.Cohort.isin(['DD (children)', 'SPARK']))], x='Phenotype', y='Phenotype Prevalence (percent)', hue='Cohort', hue_order=cohort_order, palette=palette, order=pheno_order)
ax.tick_params(axis='x', which='major', labelsize=8)
# Add significance
statdf['star']='n.s.'
statdf.loc[statdf['p value']<=0.05, 'star']='*'

cohort_shift={'DD (children)':-(4*7/80), 'SPARK':-(3*7/80),
				'DD (adults)':-(7/80), 'UKB (Questionnaire)':0, 
				'MyCode':(2*7/80), 'UKB (ICD10)':(3*7/80), 'All of Us':(4*7/80)}
pheno_centers={'ID/DD':0, 'ASD':1, 'Anxiety':2, 'Depression':3, 'Psychosis':4}

for idx, row in statdf.iterrows():
	pheno=row.Phenotype
	if pheno not in pheno_order:
		continue
	c1_shift=cohort_shift[row['Cohort 1']]
	c2_shift=cohort_shift[row['Cohort 2']]
	center=pheno_centers[pheno]
	max_prev=row[['Cohort 1 prevalence', 'Cohort 2 prevalence']].str.split('%', expand=True)[0].astype(float).max()
	plot_h=max_prev+2
	if (row['Cohort 1']=='MyCode') and (row['Cohort 2']=='All of Us'):
		plot_h=max_prev+7
	
	plt.plot([center+c1_shift, center+c2_shift], [plot_h, plot_h], color='k', lw=1)
	plt.text((center+c1_shift+center+c2_shift)/2, plot_h+1, row.star, ha='center')

plt.legend(bbox_to_anchor=(1, 1))
plt.tight_layout()
pdf.savefig()
plt.close()

# Adult only phenotypes
fig, ax=plt.subplots(figsize=(6, 4))
pheno_order=['Sleep trouble', 'Addiction']
cohort_order=['DD (adults)', 'UKB (Questionnaire)', 'blank', 'MyCode', 'UKB (ICD10)']
sns.barplot(data=df, x='Phenotype', y='Phenotype Prevalence (percent)', hue='Cohort', hue_order=cohort_order, palette=palette, order=pheno_order)
ax.tick_params(axis='x', which='major', labelsize=8)
# Add significance
statdf['star']='n.s.'
statdf.loc[statdf['p value']<=0.05, 'star']='*'

cohort_shift={'DD (adults)':-0.3, 'UKB (Questionnaire)':-0.15, 
				'MyCode':0.15, 'UKB (ICD10)':0.3}
pheno_centers={'Sleep trouble':0, 'Addiction':1}

for idx, row in statdf.iterrows():
	pheno=row.Phenotype
	if pheno not in pheno_order:
		continue
	if row['Cohort 1'] not in cohort_order:
		continue
	c1_shift=cohort_shift[row['Cohort 1']]
	c2_shift=cohort_shift[row['Cohort 2']]
	center=pheno_centers[pheno]
	max_prev=row[['Cohort 1 prevalence', 'Cohort 2 prevalence']].str.split('%', expand=True)[0].astype(float).max()
	plot_h=max_prev+2
	if (row['Cohort 1']=='MyCode') and (row['Cohort 2']=='All of Us'):
		plot_h=max_prev+7
	
	plt.plot([center+c1_shift, center+c2_shift], [plot_h, plot_h], color='k', lw=1)
	plt.text((center+c1_shift+center+c2_shift)/2, plot_h+1, row.star, ha='center')

plt.legend(bbox_to_anchor=(1, 1))
plt.tight_layout()
pdf.savefig()
plt.close()

# Child only phenotypes
phenos=list(df[(df.Cohort.isin(['DD (children)', 'SPARK'])) & (df['Phenotype Prevalence (percent)']>=20)].Phenotype.unique())
orderdf=df[(df.Cohort=='DD (children)') & (df.Phenotype.isin(phenos)) & (df.Phenotype.isin(statdf.Phenotype.to_list()))  & (~df.Phenotype.isin(['ID/DD', 'ASD', 'Anxiety']))].copy()
orderdf.sort_values(by='Phenotype Prevalence (percent)', ascending=False, inplace=True)
pheno_order=orderdf.Phenotype.to_list()

fig, ax=plt.subplots(figsize=(8, 4))
cohort_order=['DD (children)', 'SPARK']
sns.barplot(data=df, x='Phenotype', y='Phenotype Prevalence (percent)', hue='Cohort', hue_order=cohort_order, palette=palette, order=pheno_order)
ax.tick_params(axis='x', which='major', labelsize=8)
# Add significance
statdf['star']='n.s.'
statdf.loc[statdf['p value']<=0.05, 'star']='*'

cohort_shift={'DD (children)':-0.25, 'SPARK':0.25}
pheno_centers={}
for i, pheno in enumerate(pheno_order):
	pheno_centers[pheno]=i

for idx, row in statdf.iterrows():
	pheno=row.Phenotype
	if pheno not in pheno_order:
		continue
	if row['Cohort 1'] not in cohort_order:
		continue
	c1_shift=cohort_shift[row['Cohort 1']]
	c2_shift=cohort_shift[row['Cohort 2']]
	center=pheno_centers[pheno]
	max_prev=row[['Cohort 1 prevalence', 'Cohort 2 prevalence']].str.split('%', expand=True)[0].astype(float).max()
	plot_h=max_prev+2
	
	plt.plot([center+c1_shift, center+c2_shift], [plot_h, plot_h], color='k', lw=1)
	plt.text((center+c1_shift+center+c2_shift)/2, plot_h+1, row.star, ha='center')

plt.legend(bbox_to_anchor=(1, 1))
plt.tight_layout()
pdf.savefig()
plt.close()

pdf.close()