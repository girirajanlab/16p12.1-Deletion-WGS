import pandas as pd
import numpy as np

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Check the enrichment of genes with specific variant types in neuronal subtypes using Fisher's Exact tests

# Input and Output files
CELLTYPE_EXPR="/path/to/combined/preferential/expression/data.csv" # Use the combined output of script 2_Analysis preparation\Brain_Celltype_Expression\1_parse_Allen_data.py
GENE_LIST_DIR="/path/to/proband/gene/list/directory" # Use the output of script 3_Data preparation\DD cohort\2_make_genelists.py
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py

OUTPUT_FIG="/path/to/output/figure.pdf" # This plot will be the forest plot presented in Fig 3C
OUTPUT_STATS="/paht/to/output/statistics/file.csv" # These are the statistics presented in Table S2F

# Laod data
expdf=pd.read_csv(CELLTYPE_EXPR)

vars=['All_variants', 'All coding SNVs', 'Missense', 'LOF', 'Splice', 'Genes del', 'Genes dup', 'STRs', 'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del (LF)', 'Genes dup (LF)', 'STRs (LF)']
celltypes=[i for i in expdf.columns.to_list() if i not in ['gene_symbol', 'gene_id']]

stat_lst=[]
for v in vars:
	vdf=pd.read_csv(f'{GENE_LIST_DIR}/proband_genes/{v}.csv', header=None, names=['Gene'])
	for ct in celltypes:
		ctdf=expdf[['gene_id', ct]].copy()
		ctdf[v]=0
		ctdf.loc[ctdf.gene_id.isin(vdf.Gene.to_list()), v]=1
		
		count_df=ctdf[[ct, v]].groupby([ct, v]).size().to_frame()
		count_df.reset_index(inplace=True)
		count_df=count_df.pivot(index=ct, columns=v, values=0)
		count_df.fillna(0, inplace=True)
		count_df=count_df.astype(int)
		
		res=stats.fisher_exact(count_df)
		or_res=stats.contingency.odds_ratio(count_df)
		ci=or_res.confidence_interval()
		stat_lst.append([v, ct, res.statistic, ci[0], ci[1], res.pvalue])

statdf=pd.DataFrame(stat_lst, columns=['Variant', 'Cell type', 'Odds ratio', '95% C.I. lower', '95% C.I. upper', 'p value'])

# Multiple testing correction
statdf['BH FDR']=np.nan
statdf.loc[statdf.Variant.str.contains('LF'), 'BH FDR']=stats.false_discovery_control(statdf[statdf.Variant.str.contains('LF')]['p value'].to_numpy(), method='bh')
statdf.loc[~statdf.Variant.str.contains('LF'), 'BH FDR']=stats.false_discovery_control(statdf[~statdf.Variant.str.contains('LF')]['p value'].to_numpy(), method='bh')

# Save
statdf.to_csv(OUTPUT_STATS, index=False)

# Make forest plots for significant results, separated by variant type
statdf=statdf[statdf['BH FDR']<=0.05]
sig_vars=[i for i in vars if i in statdf.Variant.to_list()]

# Annotate results with broader class of cell type
statdf['class']=statdf['Cell type'].map({'Vip':'Inhibitory', 'L5/6 NP':'Excitatory', 'Sst':'Inhibitory', 'L2/3 IT':'Excitatory', 'Exc':'Excitatory',
														'Inh':'Inhibitory', 'Pvalb':'Inhibitory', 'Endo_subclass':'Endothelial', 'Endo_superclass':'Endothelial',
														'Oligo_subclass':'Oligodendrocyte', 'Oligo_superclass':'Oligodendrocyte', 'L5 IT':'Excitatory'})
statdf['class_type']=statdf['Cell type'].map({'Vip':'Sub', 'L5/6 NP':'Sub', 'Sst':'Sub', 'L2/3 IT':'Sub', 'Exc':'Super',
														'Inh':'Super', 'Pvalb':'Sub', 'Endo_subclass':'Sub', 'Endo_superclass':'Super',
														'Oligo_subclass':'Sub', 'Oligo_superclass':'Super', 'L5 IT':'Sub'})
statdf.sort_values(by=['Variant', 'class', 'class_type'], inplace=True, ascending=[True, True, False])

palette={'Inhibitory':'#40449B', 'Excitatory':'#A54198', 'Endothelial':'#77966D' , 'Oligodendrocyte':'#FF8484'}

pdf=PdfPages(OUTPUT_FIG)
for v in sig_vars:
	fig, ax = plt.subplots(figsize=(3, 6))
	subdf=statdf[statdf.Variant==v].copy()
	subdf.reset_index(inplace=True)
	
	# Points
	sns.scatterplot(data=subdf, x='Odds ratio', y='Cell type', hue='class', palette=palette, legend=False, edgecolor='None')
	
	# Confidence intervals
	for idx, row in subdf.iterrows():
		color=palette[row['class']]
		lo=row['95% C.I. lower']
		hi=row['95% C.I. upper']
		
		plt.plot([lo, hi], [idx, idx], color=color)
	
		# Also plot stars
		plt.text(hi*1.05, idx, '**', color='k', ha='center', va='center')
	
	# Add a dotted line at 1
	lo, hi = plt.ylim()
	plt.plot([1, 1], [lo, hi], color='k', ls=':', zorder=-1)
	plt.ylim(lo, hi)
	
	lo, hi = plt.xlim()
	plt.xlim(0, hi*1.05)
	
	plt.title(v)
	plt.tight_layout()
	
	pdf.savefig()
	plt.close()

pdf.close()