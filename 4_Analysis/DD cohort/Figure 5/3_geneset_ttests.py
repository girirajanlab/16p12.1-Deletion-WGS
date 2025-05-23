import pandas as pd
import numpy as np

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap

matplotlib.rcParams['pdf.fonttype'] = 42

# Compare the number of genes in each gene set present in probands with and without specific phenotypic domains

# Input and output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
GENE_ANNO="/path/to/gene/annotations.csv" # Use the output of script 2_Analysis preparation\Gene_Annotations\5_add_loeuf.py
GENE_LIST_DIR="/path/to/proband/gene/list/directory" # Use the output of script 3_Data preparation\DD cohort\3_genelist_by_proband.py

OUTPUT_STATS="/path/to/output/statistics/file.csv" # These are the statistics presented in Table SC
OUTPUT_HEATMAP="/path/to/output/heatmap/figure.pdf" # Heatmap presented in Fig 5D

# Load data
df=pd.read_csv(TABS1A)

child_domains=['Behavioral features (Child domain)', 'Psychiatric features (Child domain)', 'Nervous System Abnormalities (Child domain)', 'Congenital Anomalies (Child domain)', 'Growth/Skeletal Defects (Child domain)']
genesets=['ASD_risk_genes_TADA_FDR0.3', 'ASD_coexpression_networks_Willsey2013', 'SFARI_gene_score',
              'Developmental_delay_DDD', 'Geisinger_DBD_Tier', 'DD_G2P',
			  'CHD8_targets_Cotney2015_Sugathan2014', 'FMRP_targets_Darnell2011',
			  'Constrained_LOEUF', 'BrainExpressed_Kang2011', 'PSD_Genes2Cognition',
              'Wang_Epilepsy', 'SZDB_schizophrenia']

df=df[(df.Relationship=='Proband') & (df.WGS=='X')][['Sample']+child_domains]
df=df[~df[child_domains].isnull().all(axis=1)]

# Load in gene sets
gdf=pd.read_csv(GENE_ANNO)
for gs in genesets:
	gdf.loc[(gdf[gs]!='0') & (gdf[gs]!=0), gs]=1
	gdf.loc[(gdf[gs]=='0') | (gdf[gs]==0), gs]=0
	gdf[gs]=gdf[gs].astype(int)

# Get the number of genes each proband has in each geneset
pros=df.Sample.to_list()
gsdf=pd.DataFrame(index=pros, columns=genesets)
for pro in pros:
	p_genes=pd.read_csv(f'{GENE_LIST_DIR}/{pro}/All_variants.csv', header=None, names=['gene_id'])['gene_id'].to_list()
	pdf=gdf[gdf.gene_id.isin(p_genes)][genesets].sum(axis=0)
	gsdf.loc[pro]=pdf
	
gsdf['Sample']=gsdf.index.to_list()
df=pd.merge(df, gsdf, on='Sample')

# T-tests of the number of genes in each domain for probands with and without specific domains
def cohens_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)

statlst=[]
for cd in child_domains:
	for gs in genesets:
		subdf=df[(~df[cd].isnull())][[cd, gs]]
		
		x=subdf[subdf[cd]>0][gs].to_list()
		y=subdf[subdf[cd]==0][gs].to_list()
		
		res=stats.ttest_ind(x, y, alternative='greater')
		
		statlst.append([cd, gs, len(x), len(y), np.mean(x), np.mean(y), cohens_d(x, y), res.statistic, res.pvalue])
statdf=pd.DataFrame(statlst, columns=['Phenotype', 'Gene set', 'N with', 'N without', 'Mean with', 'Mean without', "Cohen's D", 'statistic', 'p value'])

statdf['BH FDR']=stats.false_discovery_control(statdf['p value'].to_numpy(), method='bh')

setmap={'ASD_risk_genes_TADA_FDR0.3':'ASD candidate genes (TADA)', 'ASD_coexpression_networks_Willsey2013':'ASD midfetal co-expression', 'BrainExpressed_Kang2011':'Brain-expressed genes',
		'Constrained_LOEUF':'Constrained genes (LF)', 'PSD_Genes2Cognition':'Post-synaptic density genes', 'Developmental_delay_DDD':'DD candidate genes (DDD)', 'CHD8_targets_Cotney2015_Sugathan2014':'CHD8 targets',
		'FMRP_targets_Darnell2011':'FMRP targets', 'Geisinger_DBD_Tier':'Developmental brain disorder', 'DD_G2P':'DD candidate genes (G2P)', 'SFARI_gene_score':'ASD candidate genes (SFARI)',
		'SZDB_schizophrenia':'Schizophrenia candidate genes', 'Wang_Epilepsy':'Epilepsy candidate genes'}

# Clean up annotations
statdf['Phenotype']=statdf.Phenotype.str.split(" \(", expand=True)[0]
statdf['Gene set']=statdf['Gene set'].map(setmap)

# Save
statdf.to_csv(OUTPUT_STATS, index=False)

# Make heatmaps

# Add significance annotations
statdf['star']=''
statdf.loc[statdf['p value']<=0.05, 'star']='*'
statdf.loc[statdf['BH FDR']<=0.05, 'star']='**'

plotdf=statdf.pivot(index='Gene set', columns='Phenotype', values="Cohen's D")
stardf=statdf.pivot(index='Gene set', columns='Phenotype', values="star")

phenolst=[i.split(' (')[0] for i in child_domains]
setlst=[setmap[i] for i in genesets]
plotdf=plotdf.loc[setlst, phenolst]
stardf=stardf.loc[setlst, phenolst]

biggest=max(abs(statdf["Cohen's D"].to_numpy()))

colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plotdf, cmap=cmap, vmin=-biggest, vmax=biggest, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf)
plt.tight_layout()
plt.savefig(OUTPUT_HEATMAP)
plt.close()