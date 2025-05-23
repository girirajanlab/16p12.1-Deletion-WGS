import pandas as pd

import pandas as pd
import numpy as np

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap

matplotlib.rcParams['pdf.fonttype'] = 42

# Compare burden between probands and carrier parents and 16p12.1 del carriers in the Estonian Biobank

# Input and output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
OUTPUT_HEATMAP="/path/to/output/heatmap/figure.pdf" # This plot will be the heatmap in Fig S5F
OUTPUT_STATS="/path/to/output/statistics/file.csv" # These are the statistics presented in Table S5K

# Load data
df=pd.read_csv(TABS1A)
df=df[(df.Relationship.isin(['Proband', 'Mother', 'Father'])) | (~df['Estonian Biobank Sample'].isnull())]
df=df[df['16p12.1 deletion']=='Carrier']

vars=['All coding SNVs', 'Missense', 'LOF', 'Splice', 'Enhancer', "5' UTR", 'Promoter', 'Genes del.', 'Genes dup.', 'STRs',
		'Intelligence PRS', 'SCZ PRS', 'Education PRS', 'Autism PRS',
		'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)']
df.loc[df.Relationship!='Proband', 'Relationship']='Carrier parent'
df.loc[~df['Estonian Biobank Sample'].isnull(), 'Relationship']='Estonian Biobank Sample'

df=df[['Sample', 'Relationship']+vars]

# Run t-tests
def cohens_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)

statlst=[]
for dd in ['Proband', 'Carrier parent']:
	for v in vars:
		subdf=df[(df.Relationship.isin([dd, 'Estonian Biobank Sample'])) & (~df[v].isnull())][['Relationship', v]].copy()
		
		x=subdf[subdf.Relationship==dd][v].to_list()
		y=subdf[subdf.Relationship=='Estonian Biobank Sample'][v].to_list()
		
		alt='less'
		test_type='One tailed t-test'
		if 'PRS' in v:
			alt='two-sided'
			test_type='Two tailed t-test'
		res=stats.ttest_ind(y, x, alternative=alt)
		
		statlst.append([dd, v, test_type, len(x), len(y), np.mean(x), np.mean(y), cohens_d(y, x), res.statistic, res.pvalue])
statdf=pd.DataFrame(statlst, columns=['DD cohort samples', 'Variant', 'Test', 'DD cohort N', 'Estonian N', 'DD cohort mean', 'Estonian mean', "Cohen's D", 'statistic', 'p value'])

# FDR correction
statdf['BH FDR']=np.nan
statdf.loc[~statdf.Variant.str.contains('LF'), 'BH FDR']=stats.false_discovery_control(statdf[~statdf.Variant.str.contains('LF')]['p value'].to_numpy(), method='bh')
statdf.loc[statdf.Variant.str.contains('LF'), 'BH FDR']=stats.false_discovery_control(statdf[statdf.Variant.str.contains('LF')]['p value'].to_numpy(), method='bh')

# Save
statdf.to_csv(OUTPUT_STATS, index=False)

# Make heatmaps

# Add significance annotations
statdf['star']=''
statdf.loc[statdf['p value']<=0.05, 'star']='*'
statdf.loc[statdf['BH FDR']<=0.05, 'star']='**'

plotdf=statdf.pivot(index='Variant', columns='DD cohort samples', values="Cohen's D")
stardf=statdf.pivot(index='Variant', columns='DD cohort samples', values="star")

plotdf=plotdf.loc[vars, ['Proband', 'Carrier parent']]
stardf=stardf.loc[vars, ['Proband', 'Carrier parent']]

biggest=max(abs(statdf["Cohen's D"].to_numpy()))

fig, ax=plt.subplots(figsize=(10, 10))
colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plotdf, cmap=cmap, vmin=-biggest, vmax=biggest, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf)
plt.tight_layout()
plt.savefig(OUTPUT_HEATMAP)
plt.close()