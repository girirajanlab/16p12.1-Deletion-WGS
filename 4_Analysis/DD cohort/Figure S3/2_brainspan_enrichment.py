import pandas as pd
import numpy as np

import scipy.stats as stats

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Calculate the enrichment of second hit genes compared to genes expressed in specific brain regions

# Input and output files
BRAINSPAN_EXP="/path/to/brainspan_preferential_tissue_expression_minor_epoch.tsv" # File was generated for a previous publication and available here: https://github.com/girirajanlab/16p12_RNAseq_project/blob/main/datasets/brainspan_preferential_tissue_expression_minor_epoch.tsv
GENE_ANNO="/path/to/gene/annotations.csv" # Use the output of script 2_Analysis preparation\Gene_Annotations\5_add_loeuf.py
GENE_LIST="/path/to/proband/gene/list/directory/All_variants.csv" # Use the output of script 3_Data preparation\DD cohort\2_make_genelists.py

OUTPUT_STATS="/path/to/output/statistics/file.csv" # These are the statistics presented in Table S2G
OUTPUT_FIG="/path/to/output/figure.pdf" # Figure presented in Fig S3B

# Gather Brainspan tissue expression data
bsdf=pd.read_csv(BRAINSPAN_EXP, sep = '\t')
df=pd.read_csv(GENE_ANNO)

bsdf=bsdf[bsdf.ensembl_gene_id.isin(df.gene_id.to_list())]
gene_ids=bsdf.ensembl_gene_id.to_list()
tissues=[i for i in bsdf.columns.to_list() if i!='ensembl_gene_id']

bsdf=bsdf[tissues]
bsdf=bsdf[bsdf.columns[(bsdf.sum(axis=0)>0).values].to_list()]
bsdf['ensembl_gene_id']=gene_ids
tissues=[i for i in bsdf.columns.to_list() if i!='ensembl_gene_id']

vdf=pd.read_csv(GENE_LIST, header=None, names=['gene_id'])

# Fisher's Exact tests
stat_lst=[]
for st in tissues:
	mdf=bsdf[['ensembl_gene_id', st]].copy()
	mdf['variant']=0
	mdf.loc[mdf.ensembl_gene_id.isin(vdf.gene_id.to_list()), 'variant']=1
	
	count_df=mdf[[st, 'variant']].groupby([st, 'variant']).size().to_frame()
	count_df.reset_index(inplace=True)
	count_df=count_df.pivot(index=st, columns='variant', values=0)
	count_df.fillna(0, inplace=True)
	count_df=count_df.astype(int)
	
	res=stats.fisher_exact(count_df)
	or_res=stats.contingency.odds_ratio(count_df)
	ci=or_res.confidence_interval()
	stat_lst.append([st, res.statistic, ci[0], ci[1], res.pvalue])

statdf=pd.DataFrame(stat_lst, columns=['Tissue', 'Odds ratio', '95% C.I. lower', '95% C.I. upper', 'p value'])

# FDR correction
statdf['BH FDR']=stats.false_discovery_control(statdf['p value'].to_numpy(), method='bh')

# Update annotations
# Separate tissue/timepoint annotations
regions = ['Cerebellum', 'Hippocampus', 'Amygdala', 'Striatum', 'Thalamus', 'Medial frontal cortex', 'Orbitofrontal cortex',
           'Dorsolateral frontal cortex', 'Ventrolateral frontal cortex', 'Primary motor cortex',
           'Primary somatosensory cortex', 'Inferior parietal cortex', 'Primary auditory cortex',
           'Superior temporal cortex', 'Inferior temporal cortex', 'Primary visual cortex']
times = ['Early fetal', 'Early mid-fetal', 'Late mid-fetal', 'Late fetal', 'Early infancy', 'Late infancy',
         'Early childhood', 'Middle and late childhood', 'Adolescence', 'Young adulthood', 'Middle adulthood']

statdf['Region']='.'
statdf['Time']='.'
for r in regions:
	statdf.loc[statdf.Tissue.str.contains(r.lower()), 'Region']=r
for t in times:
	statdf.loc[statdf.Tissue.str.contains(t.lower()), 'Time']=t

statdf.Region=pd.Categorical(statdf.Region, regions)
statdf.Time=pd.Categorical(statdf.Time, times)

statdf.sort_values(by=['Region', 'Time'], inplace=True)

statdf=statdf[['Region', 'Time', 'Odds ratio', '95% C.I. lower', '95% C.I. upper', 'p value', 'BH FDR']]

# Save
statdf.to_csv(OUTPUT_STATS, index=False)

statdf['rt']=statdf.Region.astype(str)+'.'+statdf.Time.astype(str)

# Log transform odds ratios
trans_cols=['Odds ratio', '95% C.I. lower', '95% C.I. upper']
for tc in trans_cols:
	statdf[tc]=np.log2(statdf[tc].to_numpy())

# If a region/time does not exist in dataset, add NA
app=[]
for r in regions:
	for t in times:
		rt=r+'.'+t
		if statdf[statdf.rt==rt].shape[0]==0:
			app.append([np.nan, np.nan, np.nan, np.nan, np.nan, rt, r, t])
appdf=pd.DataFrame(app, columns=['Odds ratio', '95% C.I. lower', '95% C.I. upper', 'p value', 'BH FDR', 'rt', 'Region', 'Time'])
statdf=pd.concat([statdf, appdf])

statdf.Region=pd.Categorical(statdf.Region, regions)
statdf.Time=pd.Categorical(statdf.Time, times)
statdf.sort_values(by=['Region', 'Time'], inplace=True)

# Function to plot a region
def plot_lines(region, color, buffer=0, ci = True):
    subdf=statdf[(statdf.Region==region)].copy()
	# Get relevant lines
    estimate = subdf['Odds ratio'].to_list()
    upper = subdf['95% C.I. lower'].to_list()
    lower = subdf['95% C.I. upper'].to_list()
    
    xs = [i for i in range(len(estimate))]

    if ci:
        plt.fill_between(xs, lower, upper, color = color, alpha = 0.08)
    plt.plot(xs, estimate, color = color)
    
    # Add * for p-values
    fdrs = subdf['BH FDR'].to_list()
    ps = subdf['p value'].to_list()

    for i, fdr in enumerate(fdrs):
        if fdr <= 0.05:
            plt.text(xs[i], buffer, '**', color = color, ha='center')
        elif ps[i] <= 0.05:
            plt.text(xs[i], buffer, '*', color = color, ha='center')

def custom_legend(colors):
    lines = []
    for i in colors:
        lines.append(Line2D([0], [0], color=i, lw=2))
    return(lines)

fig, ax=plt.subplots(figsize=(10, 8))
cmap='RdYlBu_r'
plt.plot([0, len(times)-1], [0, 0], linestyle = ':', color = 'k')
map = cm.get_cmap(cmap, 17)
buffer=statdf['95% C.I. upper'].max()
for i, region in enumerate(regions):
	plot_lines(region, map(i), buffer=buffer)
	buffer+=0.1

lo, hi=plt.ylim()
plt.ylim(lo, buffer+0.1)

# Add ticks
ax = plt.gca()
ax.set_xticks([i for i in range(len(times))])
ax.set_xticklabels(times, rotation = 90)

# Legend
custom_lines = custom_legend([map(i) for i in range(16)])
plt.legend(custom_lines, regions, bbox_to_anchor = (1, 1))

plt.tight_layout()

# Save
plt.savefig(OUTPUT_FIG)