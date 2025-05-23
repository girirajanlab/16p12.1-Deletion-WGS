import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype']=42

import scipy.stats as stats

# Compare the distributions of BMI and HC in 16p12.1 del probands for deviation from a mean 0

# Input and Output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
OUTPUT_FIG="/path/to/output/figure.pdf"
OUTPUT_STATS="/paht/to/output/statistics/file.csv"

# Load data
df=pd.read_csv(TABS1A)

# Restrict to probands with available data
df=df[df['Estonian Biobank Sample']!='X']
df=df[['Sample', 'Relationship', 'BMI Z Score', 'Head Circumference Z Score']]
df=df[(~df['BMI Z Score'].isnull()) | (~df['Head Circumference Z Score'].isnull())]
df=df[df.Relationship=='Proband']

# Plot histograms/kdes for each group for each phenotype and calculate deviation from a mean 0
pdf=PdfPages(OUTPUT_FIG)
stat_lst=[]
for pheno in ['BMI Z Score', 'Head Circumference Z Score']:
	subdf=df[~df[pheno].isnull()]
	sns.kdeplot(data=subdf, x=pheno, color='#006C67', legend=False, cut=0, common_norm=False)
	
	# Calculate SD shift and one tailed T-test p values from 0
	sd=subdf[pheno].std()

	if pheno=='Head Circumference Z Score':
		buffer=0.27
	else:
		buffer=0.29
	
	alt='greater'
	if pheno=='Head Circumference Z Score':
		alt='less'
	
	x=subdf[pheno].to_numpy()
	
	res = stats.ttest_1samp(x, popmean=0, alternative=alt)
	
	shift=abs(np.mean(x)/sd)

	# Annotate on the graph
	plt.plot([np.mean(x), 0], [buffer, buffer], color='k')
	anno='%.2f\np=%.3f' % (shift, res.pvalue)
	plt.text(np.mean(x)/2, buffer, anno, color='k', ha='center', va='center')

	# Plot vertical line at 0
	lo, hi=plt.ylim()
	plt.plot([0, 0], [lo, hi], color='red', ls='--', zorder=-1)
	plt.ylim(lo, hi)

	# Save statistics
	stat_lst.append([pheno, subdf.shape[0], res.statistic, res.pvalue, np.mean(x), shift])
	
	plt.xlabel(pheno)
	plt.title('Proband '+pheno)

	pdf.savefig()
	plt.close()
pdf.close()

# Save statistics to file
statdf=pd.DataFrame(stat_lst, columns=['Phenotype', 'Sample size', 'Statistic', 'p value', 'Proband mean', 'Shift'])
statdf.to_csv(OUTPUT_STATS, index=False)