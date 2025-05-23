import pandas as pd
import numpy as np

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap

matplotlib.rcParams['pdf.fonttype'] = 42

# Compare the burden in probands with and without specific phenotypic domains

# Input and output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
OUTPUT_HEATMAP="/path/to/output/heatmap/figure.pdf" # This plot will be the heatmap in Fig S4C
OUTPUT_STATS="/path/to/output/statistics/file.csv" # These are the statistics presented in Table S4D

# Load data
df=pd.read_csv(TABS1A)
df=df[(df.Relationship=='Proband') & (df['Estonian Biobank Sample'].isnull())]

child_domains=['Behavioral features (Child domain)', 'Psychiatric features (Child domain)', 'Nervous System Abnormalities (Child domain)', 'Congenital Anomalies (Child domain)', 'Growth/Skeletal Defects (Child domain)']
vars=['All coding SNVs', 'Missense', 'LOF', 'Splice', 'Enhancer', "5' UTR", 'Promoter', 'Genes del.', 'Genes dup.', 'STRs',
		'Intelligence PRS', 'SCZ PRS', 'Education PRS', 'Autism PRS',
		'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)']
		
# Split phenotypic domains at a value to create approximately even groups
splits={'Behavioral features (Child domain)':1, 'Psychiatric features (Child domain)':0, 'Nervous System Abnormalities (Child domain)':0, 'Congenital Anomalies (Child domain)':1, 'Growth/Skeletal Defects (Child domain)':2}
for key in splits.keys():
	df.loc[df[key]<=splits[key], key]=0
	df.loc[df[key]>splits[key], key]=1

# Run t-tests
def cohens_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)

statlst=[]
for cd in child_domains:
	for v in vars:
		subdf=df[(~df[cd].isnull()) & (~df[v].isnull())]
		
		x=subdf[subdf[cd]==1][v].to_list()
		y=subdf[subdf[cd]==0][v].to_list()
		
		res=stats.ttest_ind(x, y)
		
		statlst.append([cd, v, len(x), len(y), np.mean(x), np.mean(y), cohens_d(x, y), res.statistic, res.pvalue])
statdf=pd.DataFrame(statlst, columns=['Phenotype', 'Variant', 'N with', 'N without', 'Mean with', 'Mean without', "Cohen's D", 'statistic', 'p value'])

statdf['BH FDR']=np.nan
statdf.loc[~statdf.Variant.str.contains('LF'), 'BH FDR']=stats.false_discovery_control(statdf[~statdf.Variant.str.contains('LF')]['p value'].to_numpy(), method='bh')
statdf.loc[statdf.Variant.str.contains('LF'), 'BH FDR']=stats.false_discovery_control(statdf[statdf.Variant.str.contains('LF')]['p value'].to_numpy(), method='bh')

# Clean up annotations
statdf['Phenotype']=statdf.Phenotype.str.split(" \(", expand=True)[0]

# Save
statdf.to_csv(OUTPUT_STATS, index=False)

# Make heatmaps

# Add significance annotations
statdf['star']=''
statdf.loc[statdf['p value']<=0.05, 'star']='*'
statdf.loc[statdf['BH FDR']<=0.05, 'star']='**'

plotdf=statdf.pivot(index='Variant', columns='Phenotype', values="Cohen's D")
stardf=statdf.pivot(index='Variant', columns='Phenotype', values="star")

phenolst=[i.split(' (')[0] for i in child_domains]
plotdf=plotdf.loc[vars, phenolst]
stardf=stardf.loc[vars, phenolst]

biggest=max(abs(statdf["Cohen's D"].to_numpy()))

fig, ax=plt.subplots(figsize=(10, 10))
colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plotdf, cmap=cmap, vmin=-biggest, vmax=biggest, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf)
plt.tight_layout()
plt.savefig(OUTPUT_HEATMAP)
plt.close()