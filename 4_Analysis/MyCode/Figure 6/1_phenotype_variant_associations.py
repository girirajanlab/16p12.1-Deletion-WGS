import pandas as pd
import numpy as np

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Perform t-tests for the burden of specific variant classes between MyCode samples with and without specific phenotypes

# Input and output files
MYCODE_DATA="/path/to/MyCode/data/files.csv" # Use the output of script 3_Data preparation\MyCode\2_condense_data.py

OUTPUT_STATS="/path/to/output/statistics/file.csv" # These are the statistics presented in Table S5F
OUTPUT_HEATMAP="/path/to/output/heatmap/figure.pdf" # This plot will be the heatmap in Fig 6E middle

# Load data
df=pd.read_csv(MYCODE_DATA)

vars=['All coding SNVs', 'Missense', 'LOF', 'Genes del.', 'Genes dup.',
		'Intelligence PRS', 'SCZ PRS', 'Education PRS', 'Autism PRS',
		'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Genes del. (LF)', 'Genes dup. (LF)']
phenos=['Neoplasms', 'Blood', 'Endocrine/Metabolic', 'Mental/behavioral disorders', 'Nervous system',
		'Eye', 'Ear', 'Circulatory system', 'Respiratory system', 'Digestive system', 'Skin/subcutaeous tissue',
		'Musc. system/connective tissue', 'Genitourinary system', 'Pregnancy/childbirth', 'Congenital malformations']

# Remove splice variants from All coding SNVs
df['All coding SNVs']=df['All coding SNVs']-df.Splice
df['All coding SNVs (LF)']=df['All coding SNVs (LF)']-df['Splice (LF)']

def cohens_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)

stat_lst=[]
for v in vars:
	for p in phenos:
		subdf=df[(~df[v].isnull()) & (~df[p].isnull())]
		
		# Restrict analyses to phenotype/variant combinations where at least 5 or 10% of the cohort are available as cases and controls
		cutoff=subdf.shape[0]*0.1
		if cutoff<5:
			cutoff=5
		
		x=subdf[subdf[p]>0][v].to_numpy()
		y=subdf[subdf[p]==0][v].to_numpy()
		
		if len(x)<cutoff or len(y)<cutoff:
			print(p)
			continue
		
		res=stats.ttest_ind(x, y, alternative='two-sided')
		
		# Save stats
		stat_lst.append([p, v, 'Two tailed t-test', np.mean(x), np.std(x), len(x), np.mean(y), np.std(y), len(y), cohens_d(x, y), res.statistic, res.pvalue])
statdf=pd.DataFrame(stat_lst, columns=['Phenotype', 'Variant', 'Test', 'Case mean', 'Case SD', 'Case N', 'Control mean', 'Control SD', 'Control N', "Cohen's D", 'statistic', 'p value'])

# Multiple test correct
# Correct separately over LF and non-LF categories
statdf['BH FDR']=np.nan
statdf.loc[(statdf.Variant.str.contains('LF')) & (~statdf['p value'].isnull()), 'BH FDR']=stats.false_discovery_control(statdf[(statdf.Variant.str.contains('LF')) & (~statdf['p value'].isnull())]['p value'].to_numpy(), method='bh')
statdf.loc[(~statdf.Variant.str.contains('LF')) & (~statdf['p value'].isnull()), 'BH FDR']=stats.false_discovery_control(statdf[(~statdf.Variant.str.contains('LF')) & (~statdf['p value'].isnull())]['p value'].to_numpy(), method='bh')

# Save
statdf.to_csv(OUTPUT_STATS, index=False)

# Make a heatmap of specific results
plotdf=statdf.pivot(index='Variant', columns='Phenotype', values="Cohen's D")
statdf['star']=''
statdf.loc[statdf['p value']<=0.05, 'star']='*'
statdf.loc[statdf['BH FDR']<=0.05, 'star']='**'
stardf=statdf.pivot(index='Variant', columns='Phenotype', values="star")

pheno_order=['Neoplasms', 'Blood', 'Mental/behavioral disorders', 'Nervous system',
             'Circulatory system', 'Respiratory system', 'Skin/subcutaeous tissue',
             'Genitourinary system']

plotdf=plotdf.loc[vars, pheno_order]
stardf=stardf.loc[vars, pheno_order]

biggest=0.57

colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
fig, axs = plt.subplots(figsize=(10, 10))
sns.heatmap(data=plotdf, cmap=cmap, vmin=-biggest, vmax=biggest, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf)
plt.tight_layout()
plt.savefig(OUTPUT_HEATMAP)
plt.close()