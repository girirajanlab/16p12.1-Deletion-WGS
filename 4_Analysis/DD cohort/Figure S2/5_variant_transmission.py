import pandas as pd
import numpy as np

import scipy.stats as stats

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype']=42

# Look at tranmission of variants from each parent of probands

# Input and output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
SNVS_INHERITANCE="/path/to/SNVs/annotated/with/inheritance.txt" # Use the output of script 1_Variant calling/DD cohort/2_SNV_annotation/coding_annotations/15_inheritance_annotations.py
CNVS_INHERITANCE="/path/to/CNVs/annotated/with/inheritance.txt" # Use the output of script 1_Variant calling/DD cohort/3_CNV_calling_annotation/merge_all_cnvs/2_annotate_loeuf.py
STRS_INHERITANCE="/path/to/STRs/annotated/with/inheritance.txt" # Use the output of script 1_Variant calling/DD cohort/4_STR_calling_annotation/12_inheritance_annotation.py

OUTPUT_STATS="/path/to/output/statistics.csv" # These are the statistics presented in Table S2D
OUTPUT_FIGURE="/path/to/output/figure.pdf" # This will be the figure presented in Fig. S2F

# Load files
df=pd.read_csv(TABS1A)
df=df[df.Relationship=='Proband']

# Parse SNVs
snvs=pd.read_csv(SNVS_INHERITANCE, sep='\t')
snvs=snvs[(snvs.inheritance!='.') & (snvs.inheritance!='both')]
snvs=snvs[snvs.Sample.isin(df.Sample.to_list())]

# Parse CNVs
cnvs=pd.read_csv(CNVS_INHERITANCE, sep='\t')
cnvs=cnvs[(cnvs.inheritance!='.') & (cnvs.inheritance!='both')]
cnvs=cnvs[cnvs.Sample.isin(df.Sample.to_list())]

cnvs['Mut_type']=cnvs['Type']
cnvs['Gene_symbol']=cnvs.Gene_Symbol

# Parse STRs
strs=pd.read_csv(STRS_INHERITANCE)
strs=strs[(strs.inheritance!='.') & (strs.Sample.isin(df.Sample.to_list()))]
strs['Mut_type']='STRs'

# Combine variants
vdf=pd.concat([cnvs[['Sample', 'Mut_type', 'inheritance', 'LOEUF', 'Gene_symbol']], snvs[['Sample', 'Mut_type', 'inheritance', 'LOEUF', 'Gene_symbol']]])
vdf=pd.concat([vdf, strs[['Sample', 'Mut_type', 'inheritance', 'LOEUF', 'Gene_symbol']]])
vdf.loc[vdf.LOEUF=='.', 'LOEUF']=np.nan
vdf.LOEUF=vdf.LOEUF.astype(float)
vdf['low_LOEUF']=vdf.LOEUF<=0.35
vdf['Variant']=vdf.Mut_type.map({'missense':'Missense', 'lof':'LOF', 'splice':'Splice', 'DEL':'Genes del.', 'DUP':'Genes dup.', 'STRs':'STRs'}) + vdf.low_LOEUF.map({False:'', True:' (LF)'})

inhdf=vdf[['Sample', 'Variant', 'low_LOEUF', 'inheritance', 'Gene_symbol']].groupby(['Sample', 'Variant', 'low_LOEUF', 'inheritance']).agg('count')
inhdf.columns=['Burden']
inhdf.reset_index(inplace=True)

# Annotate mother/father
inhdf['MF']=inhdf.inheritance.str[0]
inhdf.loc[inhdf.inheritance=='de novo', 'MF']='de novo'

# Annotate carrier/noncarrier
inhdf['CNC']=inhdf.inheritance.map({'MC':'C', 'FC':'C', 'MNC':'NC', 'FNC':'NC', 'de novo':'de novo'})

# Male a list of de novo samples
dn_samps=inhdf[(inhdf.CNC=='NC') & (inhdf.Variant=='Missense') & (~inhdf.low_LOEUF)].Sample.value_counts()
dn_samps=dn_samps[dn_samps>1].index.to_list()

vars=['Missense', 'LOF', 'Splice', 'Genes del.', 'Genes dup.', 'STRs']
cols=['MF', 'CNC']
parents=[['M', 'F'], ['C', 'NC']]

par_map={'M':'Mother', 'F':'Father', 'C':'Carrier', 'NC':'Noncarrier'}

outlst=[]
for i in range(2):
	col=cols[i]
	par=parents[i]
	
	par1=par_map[par[0]]
	par2=par_map[par[1]]
	
	# Make a combined variant-parent column
	subdf=inhdf.copy()
	subdf['VP']=subdf.Variant+subdf[col]
	
	# If the analysis is carrier/noncarrier, remove de novo samples
	subdf=subdf[~subdf.Sample.isin(dn_samps)]
	
	# Pivot frame for analysis
	subdf=subdf.pivot(index='Sample', columns='VP', values='Burden')
	
	# Add in any missing columns
	for v in vars:
		for p in par+['de novo']:
			if (v+p) not in subdf.columns.to_list():
				subdf[(v+p)]=np.nan
			if (v+' (LF)'+p) not in subdf.columns.to_list():
				subdf[(v+' (LF)'+p)]=np.nan
	
	# Remove any samples missing all of one parent
	for p in par:
		subdf=subdf[~(subdf[[i+p for i in vars]+[i+' (LF)'+p for i in vars]].isnull().all(axis=1))]
	
	# Fill NA with 0
	subdf.fillna(0, inplace=True)
	subdf=subdf.astype(int)
	
	# Add the LF burden to the non-LF
	for v in vars:
		for p in par+['de novo']:
			if (v+p) not in subdf.columns.to_list():
				subdf[(v+p)]=0
			else:
				subdf[(v+p)]=subdf[v+p]+subdf[v+' (LF)'+p]
	
	# Create All coding SNVs columns
	for p in par+['de novo']:
		for lf in ['', ' (LF)']:
			subdf['All coding SNVs'+lf+p]=subdf[[v+lf+p for v in ['Missense', 'LOF', 'Splice']]].sum(axis=1)
	
	# Do a paired t-test for the number of variants inherited from either parent
	for lf in ['', ' (LF)']:
		for v in ['All coding SNVs']+vars:
			x=subdf[v+lf+par[0]].to_numpy()
			y=subdf[v+lf+par[1]].to_numpy()
		
			res=stats.ttest_rel(x, y)
			ci=res.confidence_interval()
			outlst.append([v+lf, 'Paired t-test', par1, par2, np.mean(x), np.std(x), np.mean(y), np.std(y), len(y), res.statistic, ci.low, ci.high, res.pvalue])		
	
	# Make boxplot
	pdf=PdfPages(f'Figures/1_{col}_boxplots.pdf')
	for lf in ['', ' (LF)']:
		for v in ['All coding SNVs']+vars:
			plotdf=subdf[[v+lf+par[0], v+lf+par[1]]].melt()
			plotdf['Parent']=plotdf.VP.str.replace(v, '')
			
			sns.boxplot(data=plotdf, x='Parent', y='value', width=0.5, color='white', fliersize=0)
			sns.stripplot(data=plotdf, x='Parent', y='value', jitter=0.2, color='k')
			
			plt.title(v+lf)
			
			pdf.savefig()
			plt.close()
	pdf.close()
				

outdf=pd.DataFrame(outlst, columns=['Variant', 'Test', 'Parent1', 'Parent2', 'Parent1_Mean', 'Parent1_SD', 'Parent2_Mean', 'Parent2_SD', 'N', 'Statistic', '95% CI lower', '95% CI upper', 'p value'])
# Save
outdf.to_csv(OUTPUT_STATS, index=False)


# Plot confidence intervals from paired t-tests
pdf=PdfPages(OUTPUT_FIGURE)
for pars in [['Mother', 'Father'], ['Carrier', 'Noncarrier']]:
	# Subset out relevant comparisons
	p1=pars[0]
	subdf=outdf[outdf.Parent1==p1].copy()
	subdf.reset_index(inplace=True)
	
	# Annotate order
	idx=subdf.index.to_list()
	idx.reverse()
	subdf['idx']=idx
	
	# Annotate difference between parents
	subdf['pardiff']=subdf.Parent1_Mean-subdf.Parent2_Mean
	
	# Plot
	fig, ax = plt.subplots()
	sns.scatterplot(data=subdf, x='pardiff', y='idx', color='k', marker='d')
	for index, row in subdf.iterrows():
		plt.plot([row['95% CI lower'], row['95% CI upper']], [row.idx, row.idx], color='k')
	
	# Add a line at 0
	ax.axvline(0, color='grey', ls='--', zorder=0)
	
	# Annotate rows
	plt.yticks(subdf.idx.to_list(), subdf.Variant.to_list())
	plt.ylabel('')
	
	plt.xlabel('Difference in parental burden transmission')
	
	# Add title
	plt.title('/'.join(pars))
	
	plt.tight_layout()
	
	pdf.savefig()
	plt.close()
pdf.close()