import pandas as pd
import numpy as np

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Perform paired t-tests between proband and parent burden for all variant classes

# Input and Output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
OUTPUT_BOXPLOT="/path/to/output/boxplot/figure.pdf" # This plot will be the boxplots in the bottom of Fig 3A and Fig S2E
OUTPUT_HEATMAP="/path/to/output/heatmap/figure.pdf" # This plot will be the heatmap in Fig 3A
OUTPUT_STATS="/path/to/output/statistics/file.csv" # These are the statistics presented in Table S2E

# Load data
df=pd.read_csv(TABS1A)

# Subset to needed data
df=df[df['Estonian Biobank Sample']!='X']
df=df[df.Relationship.isin(['Proband', 'Mother', 'Father'])]
df=df[df['16p12.1 deletion'].isin(['Carrier', 'Noncarrier'])]

df['Group']='Proband'
df.loc[(df.Relationship.isin(['Mother', 'Father'])) & (df['16p12.1 deletion']=='Carrier'), 'Group']='Carrier parent'
df.loc[(df.Relationship.isin(['Mother', 'Father'])) & (df['16p12.1 deletion']=='Noncarrier'), 'Group']='Noncarrier parent'

df=df[['Sample', 'Family', 'Group']+vars]

vars=['All coding SNVs', 'Missense', 'LOF', 'Splice', 'Enhancer', "5' UTR", 'Promoter', 'Genes del.', 'Genes dup.', 'STRs',
		'Intelligence PRS', 'SCZ PRS', 'Education PRS', 'Autism PRS',
		'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)']
df=df[~df[vars].isnull().all(axis=1)]

# Remove families with multiple Noncarrier parents (do novo cases)
dn_fams=df[(df.Group=='Noncarrier parent')].Family.value_counts()
dn_fams=sorted(dn_fams[dn_fams>1].index.to_list())
print(dn_fams)
df=df[~df.Family.isin(dn_fams)]

# Group by family for paired t-tests
famdf=df.melt(id_vars=['Sample', 'Family', 'Group'], value_vars=vars)
famdf=famdf.pivot(index='Family', columns=['Group', 'variable'], values='value')
print(famdf)

# Run t tests and make boxplots comparing probands and their parents
pdf=PdfPages(OUTPUT_BOXPLOT)
def cohens_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)

stat_lst=[]
for v in vars:
	for par in ['Carrier parent', 'Noncarrier parent']:
		pro_col=famdf['Proband'][v]
		par_col=famdf[par][v]
		
		subdf=pd.merge(pro_col, par_col, right_index=True, left_index=True, suffixes=['_Proband', '_'+par])
		subdf.columns=[i.split('_')[1] for i in subdf.columns.to_list()]
		subdf=subdf[(~subdf.Proband.isnull()) & (~subdf[par].isnull())]
		
		x=subdf.Proband.to_numpy()
		y=subdf[par].to_numpy()
		
		alt='greater'
		test='One tailed t-test'
		if 'PRS' in v:
			alt='two-sided'
			test='Two tailed t-test'
		
		res=stats.ttest_rel(x, y, alternative=alt)
		
		# Save stats
		stat_lst.append([v, par, test, len(x), cohens_d(x, y), res.statistic, res.pvalue])
		
		# Make boxplot
		plotdf=pd.DataFrame({v:subdf.Proband.to_list()+subdf[par].to_list(), 'Group':['Proband']*subdf.shape[0]+[par]*subdf.shape[0]})
		sns.boxplot(plotdf, x='Group', y=v, order=['Proband', par], color='white', width=0.2, fliersize=0)
		sns.stripplot(plotdf, x='Group', y=v, order=['Proband', par], marker='$\circ$', color='k')
		
		y=plotdf[v].max()
		plt.plot([0, 1], [y*1.03, y*1.03], color='k')
		plt.text(0.5, y*1.07, 'p=%.3f' % res.pvalue, color='k', ha='center', va='top')
		plt.title(v)
		plt.ylabel('Variant burden')
		plt.xlabel('')
		plt.xlim(-0.5, 1.5)
		plt.tight_layout()
		pdf.savefig()
		plt.close()
		

pdf.close()
statdf=pd.DataFrame(stat_lst, columns=['Variant', 'Parent', 'Test', 'Sample size (pairs)', "Cohen's D", 'statistic', 'p value'])

# Multiple test correct
# Correct separately over LF and non-LF categories
statdf['BH FDR']=np.nan
statdf.loc[statdf.Variant.str.contains('LF'), 'BH FDR']=stats.false_discovery_control(statdf[statdf.Variant.str.contains('LF')]['p value'].to_numpy(), method='bh')
statdf.loc[~statdf.Variant.str.contains('LF'), 'BH FDR']=stats.false_discovery_control(statdf[~statdf.Variant.str.contains('LF')]['p value'].to_numpy(), method='bh')

# Save statistics
statdf.to_csv(OUTPUT_STATS, index=False)

# Make heatmap of results
plotdf=statdf.pivot(index='Parent', columns='Variant', values="Cohen's D")
statdf['star']=''
statdf.loc[statdf['p value']<=0.05, 'star']='*'
statdf.loc[statdf['BH FDR']<=0.05, 'star']='**'
stardf=statdf.pivot(index='Parent', columns='Variant', values="star")

plotdf=plotdf[vars]
stardf=stardf[vars]

biggest=max(abs(statdf["Cohen's D"].to_numpy()))

colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plotdf, cmap=cmap, vmin=-biggest, vmax=biggest, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf)
plt.tight_layout()
plt.savefig(OUTPUT_HEATMAP)
plt.close()