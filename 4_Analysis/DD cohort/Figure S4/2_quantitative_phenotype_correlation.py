import pandas as pd
import numpy as np

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap

matplotlib.rcParams['pdf.fonttype'] = 42

# Correlate burden and quantiative phenotypes in probands

# Input and output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
OUTPUT_HEATMAP="/path/to/output/heatmap/figure.pdf" # This plot will be the heatmap in Fig S4D
OUTPUT_STATS="/path/to/output/statistics/file.csv" # These are the statistics presented in Table S4E

# Load data
df=pd.read_csv(TABS1A)
df=df[(df.Relationship=='Proband') & (df['Estonian Biobank Sample'].isnull())]

quant_phenos=['BMI Z Score', 'Head Circumference Z Score', 'HRS-MAT', 'SRS Raw Score']
vars=['All coding SNVs', 'Missense', 'LOF', 'Splice', 'Enhancer', "5' UTR", 'Promoter', 'Genes del.', 'Genes dup.', 'STRs',
		'Intelligence PRS', 'SCZ PRS', 'Education PRS', 'Autism PRS',
		'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)']

# Run Pearon correlations
statlst=[]
for qp in quant_phenos:
	for v in vars:
		subdf=df[(~df[v].isnull()) & (~df[qp].isnull())]
		
		x=subdf[qp].to_numpy()
		y=subdf[v].to_numpy()
		
		res=stats.pearsonr(x, y)
		
		statlst.append([qp, v, len(x), res.statistic, res.pvalue])
statdf=pd.DataFrame(statlst, columns=['Phenotype', 'Variant', 'N', 'Pearson R', 'p value'])

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

plotdf=statdf.pivot(index='Variant', columns='Phenotype', values="Pearson R")
stardf=statdf.pivot(index='Variant', columns='Phenotype', values="star")

plotdf=plotdf.loc[vars, quant_phenos]
stardf=stardf.loc[vars, quant_phenos]

biggest=max(abs(statdf["Pearson R"].to_numpy()))

fig, ax=plt.subplots(figsize=(10, 10))
colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plotdf, cmap=cmap, vmin=-biggest, vmax=biggest, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf)
plt.tight_layout()
plt.savefig(OUTPUT_HEATMAP)
plt.close()