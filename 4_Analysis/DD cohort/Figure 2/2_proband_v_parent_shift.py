import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype']=42

import scipy.stats as stats

# Compare the SRS and HRS-MAT distribution between 16p12.1 deletion probands and their parents

# Input and Output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
OUTPUT_FIG="/path/to/output/figure.pdf"
OUTPUT_STATS="/paht/to/output/statistics/file.csv"

# Subset needed data
df=pd.read_csv('../../../1_Variant_preparation/Table_S1A.csv')
df=df[df['Estonian Biobank Sample']!='X']
df=df[df.Relationship.isin(['Proband', 'Mother', 'Father'])]
df=df[df['16p12.1 deletion'].isin(['Carrier', 'Noncarrier'])]

# Clean up annotations
df['Group']='Proband'
df.loc[(df.Relationship.isin(['Mother', 'Father'])) & (df['16p12.1 deletion']=='Carrier'), 'Group']='Carrier parent'
df.loc[(df.Relationship.isin(['Mother', 'Father'])) & (df['16p12.1 deletion']=='Noncarrier'), 'Group']='Noncarrier parent'

df=df[['Sample', 'Relationship', '16p12.1 deletion', 'Group', 'SRS Raw Score', 'HRS-MAT']]
df=df[(~df['SRS Raw Score'].isnull()) | (~df['HRS-MAT'].isnull())]

# Plot histograms/kdes for each group for each phenotype and calculate Mann Whitney statistics
pdf=PdfPages(OUTPUT_FIG)
stat_lst=[]
for pheno in ['SRS Raw Score', 'HRS-MAT']:
	subdf=df[~df[pheno].isnull()]
	
	sns.kdeplot(data=subdf, x=pheno, hue='Group', hue_order=['Proband', 'Carrier parent', 'Noncarrier parent'], palette=['#006C67', '#DA5334', '#3C8CCF'], legend=False, cut=0, common_norm=False)

	# Calculate SD shift and one tailed Mann Whitney p values
	sd=subdf[subdf.Group=='Proband'][pheno].std()

	if pheno=='SRS Raw Score':
		buffer=0.015
	else:
		buffer=0.034
		
	for par in ['Carrier parent', 'Noncarrier parent']:
		x=subdf[subdf.Group=='Proband'][pheno].to_numpy()
		y=subdf[subdf.Group==par][pheno].to_numpy()
		
		alt='greater'
		if pheno=='HRS-MAT':
			alt='less'
		
		res = stats.mannwhitneyu(x, y, alternative=alt)
		
		shift=abs((np.mean(x)-np.mean(y))/sd)

		# Annotate on the graph
		plt.plot([np.mean(x), np.mean(y)], [buffer, buffer], color={'Carrier parent':'#DA5334', 'Noncarrier parent':'#3C8CCF'}[par])
		anno='%.2f\np=%.2E' % (shift, res.pvalue)
		plt.text((np.mean(x)-np.mean(y))/2+np.mean(y), buffer, anno, color='k', ha='center', va='center')

		# Save statistics
		stat_lst.append([pheno, par, len(x), len(y), res.statistic, res.pvalue, np.mean(x), np.mean(y), shift])
		
		if pheno=='SRS Raw Score':
			buffer+=0.002
		else:
			buffer+=0.005

	plt.xlabel(pheno)
	plt.title('Proband vs. Parent '+pheno.upper())

	pdf.savefig()
	plt.close()
pdf.close()

# Save statistics to file
statdf=pd.DataFrame(stat_lst, columns=['Phenotype', 'Parent', 'Proband n', 'Parent n', 'Statistic', 'p value', 'Proband mean', 'Parent mean', 'Shift'])
statdf.to_csv(OUTPUT_STATS, index=False)