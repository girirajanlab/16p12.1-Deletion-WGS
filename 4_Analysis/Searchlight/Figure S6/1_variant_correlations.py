import pandas as pd
import numpy as np

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
from matplotlib.colors import LinearSegmentedColormap

matplotlib.rcParams['pdf.fonttype'] = 42

# Correlate variant burden and quantitative phenotypes in SSC probands

# Input and output files
DATA_DIR="/path/to/Searchlight/data/directory/" # Use output directory of 3_Data preparation\Searchlight\1_compile_data.py

OUTPUT_STATS="/path/to/output/statistics.csv" # Statistics presented in Table S6E
OUTPUT_HEATMAP="/pathh/to/output/heatmaps.pdf" # Forest plots presented in Fig S6B

# Define variables
quant_phenos=['Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Social responsiveness (SRS)', 'Autism behavior (BSI)', 'BMI z-score', 'Head circumference z-score']
cohorts=['16p11.2 deletion', '16p11.2 duplication']
vars=['All coding SNVs', 'Missense', 'LOF', 'Splice', 'Genes del.', 'Genes dup.',
		'Intelligence PRS', 'SCZ PRS', 'Education PRS', 'Autism PRS',
		'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del. (LF)', 'Genes dup. (LF)']

statdf=pd.DataFrame()
for cohort in cohorts:
	df=pd.read_csv(f'{DATA_DIR}/{cohort}.csv')

	statlst=[]
	for v in vars:
		for qp in quant_phenos:
			subdf=df[(~df[v].isnull()) & (~df[qp].isnull())]
			
			# Pearson R
			x=subdf[v].to_numpy()
			y=subdf[qp].to_numpy()
			
			res=stats.pearsonr(x, y)
			
			statlst.append([cohort, v, qp, subdf.shape[0], res.statistic, res.pvalue])
	cohort_df=pd.DataFrame(statlst, columns=['Cohort', 'Variant', 'Phenotype', 'Sample size', 'Pearson R', 'p value'])
	
	# Multiple testing correction
	cohort_df['BH FDR']=np.nan
	cohort_df.loc[~cohort_df.Variant.str.contains('LF'), 'BH FDR']=stats.false_discovery_control(cohort_df[~cohort_df.Variant.str.contains('LF')]['p value'].to_numpy(), method='bh')
	cohort_df.loc[cohort_df.Variant.str.contains('LF'), 'BH FDR']=stats.false_discovery_control(cohort_df[cohort_df.Variant.str.contains('LF')]['p value'].to_numpy(), method='bh')
	
	statdf=pd.concat([statdf, cohort_df])

# Save
statdf.to_csv(OUTPUT_STATS, index=False)

# Make heatmaps

# Add significance annotations
statdf['star']=''
statdf.loc[statdf['p value']<=0.05, 'star']='*'
statdf.loc[statdf['BH FDR']<=0.05, 'star']='**'

pdf=PdfPages(OUTPUT_HEATMAP)
for cohort in cohorts:
	subdf=statdf[statdf.Cohort==cohort]
	
	plotdf=subdf.pivot(index='Variant', columns='Phenotype', values="Pearson R")
	stardf=subdf.pivot(index='Variant', columns='Phenotype', values="star")

	plotdf=plotdf.loc[vars, quant_phenos]
	stardf=stardf.loc[vars, quant_phenos]

	biggest=max(abs(statdf["Pearson R"].to_numpy()))

	fig, ax=plt.subplots(figsize=(10, 10))
	colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
	cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
	cmap.set_bad('#CCCCCC')
	sns.heatmap(data=plotdf, cmap=cmap, vmin=-biggest, vmax=biggest, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf)
	plt.tight_layout()
	pdf.savefig()
	plt.close()
pdf.close()
