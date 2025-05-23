import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype']=42

import scipy.stats as stats

import random

# Compare the SRS distribution between 16p12.1 deletion probands and SSC probands and 16p11.2 CNV probands from Simons Searchlight

# Input and Output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
SSC_DATA="/path/to/SSC/cohort/data.csv" # Use the output of script 3_Data preparation/SSC/3_merge_data.py
SEARCHLIGHT_DATA="/path/to/Searchlight/SRS/data.csv" # Use the output of script 3_Data preparation/Searchlight/3_gather_SRS.py

OUTPUT_FIG="/path/to/output/figure.pdf"
OUTPUT_SRS="/path/to/output/SRS/statistics.csv" # This file will contain the output statistics of SRS comparisons
OUTPUT_HRS_MAT="/path/to/output/HRS-MAT/statistics.csv" # This file will contain the output statistics of comparing the proband HRS-MAT distribution to the SSC cohort mean reported in Hansen 2019 J Autism Dev Disord (https://pubmed.ncbi.nlm.nih.gov/27738852/)

# 16p12.1 deletion proband data
df=pd.read_csv(TABS1A)
df=df[df['Estonian Biobank Sample']!='X']
df=df[['Sample', 'Relationship', 'SRS Raw Score']]
df=df[(~df['SRS Raw Score'].isnull())]
df=df[df.Relationship=='Proband']
df['Cohort']='16p12.1 del. probands'

# SSC proband data
ssc=pd.read_csv(SSC_DATA)
ssc=ssc[~ssc['Social responsiveness (SRS)'].isnull()]
ssc['SRS Raw Score']=ssc['Social responsiveness (SRS)']
ssc['Cohort']='Autism probands'

# 16p11.2 proband data
searchlight=pd.read_csv(SEARCHLIGHT)
searchlight['Cohort']='16p11.2 CNV probands'

df=pd.concat([df, ssc, searchlight])

# Plot kdes for each group for SRS and calculate Mann Whitney statistics
pdf=PdfPages(OUTPUT_FIG)
stat_lst=[]

sns.kdeplot(data=df, x='SRS Raw Score', hue='Cohort', hue_order=['16p12.1 del. probands', 'Autism probands', '16p11.2 CNV probands'], palette=['#006C67', '#F6921E', '#5B12B3'], legend=False, cut=0, common_norm=False)

# Calculate SD shift and one tailed Mann Whitney p values
sd=df[df.Cohort=='16p12.1 del. probands']['SRS Raw Score'].std()

buffer=0.015
for par in ['Autism probands', '16p11.2 CNV probands']:
	x=df[df.Cohort=='16p12.1 del. probands']['SRS Raw Score'].to_numpy()
	y=df[df.Cohort==par]['SRS Raw Score'].to_numpy()
	
	res = stats.mannwhitneyu(x, y, alternative='two-sided')
	
	shift=abs((np.mean(x)-np.mean(y))/sd)

	# Annotate on the graph
	plt.plot([np.mean(x), np.mean(y)], [buffer, buffer], color={'Autism probands':'#F6921E', '16p11.2 CNV probands':'#5B12B3'}[par])
	anno='%.2f\np=%.2E' % (shift, res.pvalue)
	plt.text((np.mean(x)-np.mean(y))/2+np.mean(y), buffer, anno, color='k', ha='center', va='center')
	
	buffer+=0.0015

	# Save statistics
	stat_lst.append(['SRS Raw Score', par, len(x), len(y), res.statistic, res.pvalue, np.mean(x), np.mean(y), shift])

plt.xlabel('SRS Raw Score')
plt.title('16p12.1 del. probands vs. SSC and SVIP SRS Raw Score')

pdf.savefig()
plt.close()

# Save statistics to file
statdf=pd.DataFrame(stat_lst, columns=['Phenotype', 'Autsim cohort', 'Proband n', 'Cohort n', 'Statistic', 'p value', '16p12.1 proband mean', 'Autism proband mean', 'Shift'])
statdf.to_csv(OUTPUT_SRS, index=False)


# Also compare 16p12.1 deletion proband HRS-MAT scores to mean SSC proband scores (mu=99) from Hansen 2019 J Autism Dev Disord (https://pubmed.ncbi.nlm.nih.gov/27738852/)
df=pd.read_csv(TABS1A)
df=df[['Sample', 'Relationship', 'HRS-MAT']]
df=df[(~df['HRS-MAT'].isnull())]
df=df[df.Relationship=='Proband']
df['Cohort']='16p12.1 del. probands'

sns.kdeplot(data=df, x='HRS-MAT', color='#006C67', legend=False, cut=0)

# One sample t-test
x=df['HRS-MAT'].to_numpy()
sd=df['HRS-MAT'].std()
res = stats.ttest_1samp(x, popmean=99, alternative='less')

# Calculate SD shift
shift=abs((np.mean(x)-99)/sd)

statdf=pd.DataFrame([['HRS-MAT', len(x), res.statistic, res.pvalue, np.mean(x), shift]], columns=['Phenotype', 'Sample size', 'Statistic', 'p value', 'Proband mean', 'Shift'])
print(statdf)
statdf.to_csv(OUTPUT_HRS_MAT, index=False)

# Annotate on the graph
buffer=0.027
plt.plot([np.mean(x), 99], [buffer, buffer], color='k')
anno='%.2f\np=%.3f' % (shift, res.pvalue)
plt.text((np.mean(x)-99)/2+99, buffer, anno, color='k', ha='center', va='center')

# Plot vertical line at SSC mean (99)
lo, hi=plt.ylim()
plt.plot([99, 99], [lo, hi], color='#F6921E', ls='--', zorder=-1)
plt.ylim(lo, hi)

plt.xlabel('HRS-MAT')
plt.title('16p12.1 del. probands vs. Mean SSC HRS-MAT')

pdf.savefig()
plt.close()

pdf.close()
