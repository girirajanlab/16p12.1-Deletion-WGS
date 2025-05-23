import pandas as pd
import numpy as np

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Check the enrichment of genes with specific variant types in 16p12.1 deletion NPC co-expression modules using Fisher's Exact tests

# Input and output files
MODULES="/path/to/WGCNA/module/genes.csv" # Use the gene output from script 2_Analysis preparation\NPC_RNAseq\2_WGCNA.R
GENE_ANNO="/path/to/gene/annotations.csv" # Use the output from script 2_Analysis preparation\Gene_Annotations\5_add_loeuf.py
GENE_LIST_DIR="/path/to/proband/gene/list/directory" # Use the output of script 3_Data preparation\DD cohort\2_make_genelists.py

OUTPUT_FIG="/path/to/output/figure.pdf" # This plot will be the heatmap presented in Fig 3D
OUTPUT_STATS="/path/to/output/statistics/file.csv" # These are the statistics presented in Table S2E

# Load in WGCNA modules
modules=pd.read_csv(MODULES)

# Identify modules with 16p12.1 deletion genes
genes_16p12={'UQCRC2':'ENSG00000140740',
				'PDZD9':'ENSG00000155714',
				'MOSMO':'ENSG00000185716',
				'VWA3A':'ENSG00000175267',
				'SDR42E2':'ENSG00000183921',
				'EEF2K':'ENSG00000103319',
				'POLR3E':'ENSG00000058600',
				'CDR2':'ENSG00000140743'}
gene_ids_16p12=[genes_16p12[i] for i in genes_16p12.keys()]
modules['16p12.1 gene']=modules.Gene.isin(gene_ids_16p12)

modules_16p12=list(modules[modules['16p12.1 gene']].Module.unique())

# Remove the grey module
modules_16p12.remove('grey')

modules['placeholder']=1
df=pd.pivot(modules, index='Gene', columns='Module', values='placeholder')

# Add in all possible genes
genedf=pd.read_csv(GENE_ANNO)
add_genes=genedf[~genedf.gene_id.isin(df.index.to_list())]['gene_id'].to_list()
add_df=pd.DataFrame(index=add_genes)

df=pd.concat([df, add_df])
df.fillna(0, inplace=True)
df=df.astype(int)

# Run Fisher's Exact tests
vars=['All_variants', 'All coding SNVs', 'Missense', 'LOF', 'Splice',
       'Genes del', 'Genes dup', 'STRs',
       'All coding SNVs (LF)']

stat_lst=[]
for v in vars:
	vdf=pd.read_csv(f'{GENE_LIST_DIR}/proband_genes/{v}.csv', header=None, names=['Gene_id'])
	for m in modules_16p12:
		mdf=df[[m]].copy()
		mdf[v]=0
		mdf.loc[vdf.Gene_id.to_list(), v]=1
		
		count_df=mdf[[m, v]].groupby([m, v]).size().to_frame()
		count_df.reset_index(inplace=True)
		count_df=count_df.pivot(index=m, columns=v, values=0)
		count_df.fillna(0, inplace=True)
		count_df=count_df.astype(int)
		
		res=stats.fisher_exact(count_df)
		or_res=stats.contingency.odds_ratio(count_df)
		ci=or_res.confidence_interval()
		stat_lst.append([m, v, res.statistic, ci[0], ci[1], res.pvalue])
		
statdf=pd.DataFrame(stat_lst, columns=['Module', 'Variant', 'Odds ratio', '95% C.I. lower', '95% C.I. upper', 'p value'])

# Multiple testing correction
statdf['BH FDR']=stats.false_discovery_control(statdf['p value'].to_numpy(), method='bh')

# Save results
statdf.to_csv(OUTPUT_STATS, index=False)

# Create a heatmap of LCL enrichment results
statdf['star']=''
statdf.loc[statdf['p value']<=0.05, 'star']='*'
statdf.loc[statdf['BH FDR']<=0.05, 'star']='**'

plotdf=statdf.pivot(index='Module', columns='Variant', values="Odds ratio")
stardf=statdf.pivot(index='Module', columns='Variant', values="star")

plotdf=plotdf.loc[modules_16p12, vars]
stardf=stardf.loc[modules_16p12, vars]

colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plotdf, cmap=cmap, vmin=0, vmax=2.9, center=1, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf, annot_kws={'fontsize':5})
plt.tight_layout()
plt.savefig(OUTPUT_FIG)
plt.close()

df_16p=modules[modules.Gene.isin(gene_ids_16p12)].copy()
reversed_dict = {v: k for k, v in genes_16p12.items()}
df_16p['gene_symbol']=df_16p.Gene.map(reversed_dict)

df_16p.sort_values(by=['Module', 'gene_symbol'], inplace=True)

print(df_16p)