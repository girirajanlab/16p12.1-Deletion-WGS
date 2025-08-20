import pandas as pd

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap

matplotlib.rcParams['pdf.fonttype'] = 42

# Calculate enrichment of secondary variants in curated disease gene sets

# INput and output files
GENE_ANNO="/path/to/gene/annotations.csv" # Use the output of script 2_Analysis preparation\Gene_Annotations\5_add_loeuf.py
GENE_LIST_DIR="/path/to/proband/gene/list/directory" # Use the output of script 3_Data preparation\DD cohort\2_make_genelists.py

OUTPUT_STATS="/path/to/output/statistics/file.csv" # These are the statistics presented in Table S2F
OUTPUT_FIG="/path/to/output/heatmap/figure.pdf" # Heatmap presented in Fig S3A

# Load in gene sets
genesets=['ASD_risk_genes_TADA_FDR0.3', 'ASD_coexpression_networks_Willsey2013', 'BrainExpressed_Kang2011',
              'Constrained_LOEUF', 'PSD_Genes2Cognition', 'Developmental_delay_DDD', 'CHD8_targets_Cotney2015_Sugathan2014',
              'FMRP_targets_Darnell2011', 'Geisinger_DBD_Tier', 'DD_G2P', 'SFARI_gene_score', 'SZDB_schizophrenia', 'Wang_Epilepsy']
vars=['Missense', 'LOF', 'Splice',
       'Genes del', 'Genes dup', 'STRs']
	   
gdf=pd.read_csv(GENE_ANNO)
for gs in genesets:
	gdf.loc[(gdf[gs]!='0') & (gdf[gs]!=0), gs]=1
	gdf.loc[(gdf[gs]=='0') | (gdf[gs]==0), gs]=0
	gdf[gs]=gdf[gs].astype(int)

stat_lst=[]
for v in vars:
	vdf=pd.read_csv(f'{GENE_LIST_DIR}/{v}.csv', header=None, names=['gene_id'])
	for gs in genesets:
		mdf=gdf[['gene_id', gs]].copy()
		mdf[v]=0
		mdf.loc[mdf.gene_id.isin(vdf.gene_id.to_list()), v]=1
		
		count_df=mdf[[gs, v]].groupby([gs, v]).size().to_frame()
		count_df.reset_index(inplace=True)
		count_df=count_df.pivot(index=gs, columns=v, values=0)
		count_df.fillna(0, inplace=True)
		count_df=count_df.astype(int)
		
		res=stats.fisher_exact(count_df)
		or_res=stats.contingency.odds_ratio(count_df)
		ci=or_res.confidence_interval()
		stat_lst.append([gs, v, res.statistic, ci[0], ci[1], res.pvalue])

statdf=pd.DataFrame(stat_lst, columns=['Gene set', 'Variant', 'Odds ratio', '95% C.I. lower', '95% C.I. upper', 'p value'])

# Multiple testing correction
statdf['BH FDR']=stats.false_discovery_control(statdf['p value'].to_numpy(), method='bh')

# Clean up annotations
setmap={'ASD_risk_genes_TADA_FDR0.3':'ASD candidate genes (TADA)', 'ASD_coexpression_networks_Willsey2013':'ASD midfetal co-expression', 'BrainExpressed_Kang2011':'Brain-expressed genes',
		'Constrained_LOEUF':'Constrained genes (LF)', 'PSD_Genes2Cognition':'Post-synaptic density genes', 'Developmental_delay_DDD':'DD candidate genes (DDD)', 'CHD8_targets_Cotney2015_Sugathan2014':'CHD8 targets',
		'FMRP_targets_Darnell2011':'FMRP targets', 'Geisinger_DBD_Tier':'Developmental brain disorder', 'DD_G2P':'DD candidate genes (G2P)', 'SFARI_gene_score':'ASD candidate genes (SFARI)',
		'SZDB_schizophrenia':'Schizophrenia candidate genes', 'Wang_Epilepsy':'Epilepsy candidate genes'}
statdf['Gene set']=statdf['Gene set'].map(setmap)

# Plot a heatmap of disease gene set enrichment by variant type
statdf.to_csv(OUTPUT_STATS, index=False)

# Add significance annotations
statdf['star']=''
statdf.loc[statdf['p value']<=0.05, 'star']='*'
statdf.loc[statdf['BH FDR']<=0.05, 'star']='**'

plotdf=statdf.pivot(index='Gene set', columns='Variant', values="Odds ratio")
stardf=statdf.pivot(index='Gene set', columns='Variant', values="star")

setlst=['ASD candidate genes (TADA)', 'ASD midfetal co-expression', 'Brain-expressed genes', 'Constrained genes (LF)','Post-synaptic density genes', 'DD candidate genes (DDD)', 'CHD8 targets',
		'FMRP targets', 'Developmental brain disorder', 'DD candidate genes (G2P)', 'ASD candidate genes (SFARI)', 'Schizophrenia candidate genes', 'Epilepsy candidate genes']
plotdf=plotdf.loc[setlst, vars]
stardf=stardf.loc[setlst, vars]

colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plotdf, cmap=cmap, vmin=0, vmax=3.25, center=1, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf, annot_kws={'fontsize':5})
plt.tight_layout()
plt.savefig(OUTPUT_FIG)
plt.close()