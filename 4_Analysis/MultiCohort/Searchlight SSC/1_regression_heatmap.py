import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Make a heatmap of the results from SSC and Searchlight

# Input and output files
SSC_RES="/path/to/SSC/regression/results.csv" # Use the output of script 4_Analysis\SSC\Figure 7\1_regressions.py
SEARCHLIGHT_RES="/path/to/Searchlight/regression/results.csv" # Use the output of script 4_Analysis\Searchlight\Figure 7\1_regressions.py

OUTPUT_FIG="/path/to/output/heatmap.pdf" # Figure presented in Fig 7A

# Set variables
quant_phenos=['Full scale IQ', 'Social responsiveness (SRS)', 'Coordination disorder (DCDQ)', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Repetitive behavior (RBS-R)',
				'Autism behavior (BSI)', 'BMI z-score', 'Head circumference z-score']
cohorts=['16p11.2 deletion', '16p11.2 duplication', 'DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'No primary variant']
vars=['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs',
		'All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)', 'SCZ PRS']

# Load data
ssc=pd.read_csv(SSC_RES)
searchlight=pd.read_csv(SEARCHLIGHT_RES)

df=pd.concat([ssc, searchlight])

mod_map={'Phenotype ~ Sex + All coding SNVs + Genes del. + Genes dup. + STRs + SCZ PRS':'B',
			'Phenotype ~ Sex + All coding SNVs (LF) + Genes del. (LF) + Genes dup. (LF) + STRs (LF) + SCZ PRS':'C',
			'Phenotype ~ Sex + All coding SNVs + Genes del. + Genes dup. + SCZ PRS':'B',
			'Phenotype ~ Sex + All coding SNVs (LF) + Genes del. (LF) + Genes dup. (LF) + SCZ PRS':'C'}

# Update model names for consistency
df['Model']=df.Model.map(mod_map)

# Remove covariates and intercepts
df=df[~df.Variable.isin(['Sex', 'Intercept'])]

# Restrict to variable/model/phenotype combinations with significant results
df['var_mod_pheno']=df.Variable+'.'+df.Model+'.'+df.Phenotype
sig_vmp=df[df['p value']<=0.05]['var_mod_pheno'].to_list()

df=df[df.var_mod_pheno.isin(sig_vmp)]

# Add significance annotations
df['star']=''
df.loc[df['p value']<=0.05, 'star']='*'

# Sort
df.Phenotype=pd.Categorical(df.Phenotype, quant_phenos)
df.Variable=pd.Categorical(df.Variable, vars)

df.sort_values(by=['Phenotype', 'Model', 'Variable'], inplace=True)

# Pivot for plotting
plotdf=pd.pivot(df, index=['Phenotype', 'Model', 'Variable'], columns='Cohort', values='Estimate')
stardf=pd.pivot(df, index=['Phenotype', 'Model', 'Variable'], columns='Cohort', values='star')

plotdf.reset_index(inplace=True)
stardf.reset_index(inplace=True)

plotdf.sort_values(by=['Phenotype', 'Model', 'Variable'], inplace=True)
stardf.sort_values(by=['Phenotype', 'Model', 'Variable'], inplace=True)

plotdf.index=plotdf[['Phenotype', 'Model', 'Variable']]
stardf.index=stardf[['Phenotype', 'Model', 'Variable']]

plotdf=plotdf[cohorts]
stardf=stardf[cohorts]

biggest=df.Estimate.abs().max()

colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plotdf, cmap=cmap, vmin=-biggest, vmax=biggest, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf)
plt.tight_layout()
plt.savefig(OUTPUT_FIG)
plt.close()
