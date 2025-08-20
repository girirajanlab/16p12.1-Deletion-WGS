import pandas as pd
import numpy as np

from collections import Counter

import random

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams['pdf.fonttype'] = 42

random.seed(205)

# Define connectivity quartiles for all genes in the parsed STRING DB

# Input and output files
STRING="/path/to/annotated/string/data.csv.gz" # Use the output of script 2_Analysis preparation\STRING_Network\1_parse_string.py
BACKGROUND="/path/to/string/background/genes.csv" # Use the output of script 2_Analysis preparation\STRING_Network\2_string_connectivity.py
ALL_VARS="/path/to/all/proband/variant/genes.csv" # Use the output of script 3_Data preparation/DD cohort/2_make_genelists.py with a single gene list for all proband variants

OUTPUT_FIGURE="/path/to/output/figure.pdf" # Figure presented in Fig. S4D
OUTPUT_STATISTICS="/path/to/output/statistics.csv" # Stats presented in Table S3D

# Load STRING relationships
string=pd.read_csv(STRING, compression='gzip')
string_genes=string.Gene1.to_list()+string.Gene2.to_list()

# Load all genes
background=pd.read_csv(BACKGROUND, header=None, names=['GeneID'])

# Filter the background gene list to include only those present in STRING
background=background[background.GeneID.isin(string_genes)]

# Filter STRING to highest confidence interactions
string.new_combined_score=string.new_combined_score/1000
string=string[string.new_combined_score>=0.9]

# Get a counter of the total number of connections a gene is involved in
connections=dict(Counter(string.Gene1.to_list()+string.Gene2.to_list()))
background['Connections']=background.GeneID.map(connections)
background.fillna(0, inplace=True)
background['Connections']=background.Connections.astype(int)

# Assign connectivity quartiles
background['Quantile']=0
quants=[0.25, 0.5, 0.75, 1]
for q in reversed(quants):
	quantile=background.Connections.quantile(q=q)
	background.loc[background.Connections<=quantile, 'Quantile']=q

hitdf=pd.read_csv(ALL_VARS, header=None, names=['GeneID'])
# Make sure secondary variant is in network
hitdf=hitdf[hitdf.GeneID.isin(string_genes)]
secondary_variants=hitdf.GeneID.to_list()

# Run simulation
df=pd.DataFrame(background[background.GeneID.isin(secondary_variants)].Quantile.value_counts())
df.columns=['Real']
df=df.transpose()

df=df[[i for i in quants if i in df.columns.to_list()]]

# Create random simulations
num_genes=len(lst)
for i in range(1000):
	gene_list=background.GeneID.to_list()
	random.shuffle(gene_list)
	sim_genes=gene_list[0:num_genes]
	sim_vals=pd.DataFrame(background[background.GeneID.isin(sim_genes)].Quantile.value_counts())
	sim_vals.columns=[str(i)]
	sim_vals=sim_vals.transpose()
	df=pd.concat([df, sim_vals])

df.fillna(0, inplace=True)

# Get empirical p values for real values
stat_lst=[]
for q in quants:
	real=df.loc['Real', q]
	sim=df[df.index!='Real'][q].to_numpy()
	
	if real<=np.mean(sim):
		emp_p=sum(sim<=real)
	else:
		emp_p=sum(sim>=real)
	
	stat_lst.append([q, real, np.mean(sim), emp_p/1000])
statdf=pd.DataFrame(stat_lst, columns=['Quantile', 'Actual value', 'Simulation mean', 'Empirical p value'])

# Save
statdf.to_csv(OUTPUT_STATISTICS, index=False)

# Make violin plots comparing real and simulated values
df['idx']=df.index.to_list()
longdf=df.melt(id_vars='idx')
plotdf=longdf[longdf.idx!='Real']

# Violins
sns.violinplot(data=plotdf, x='Quantile', y='value', order=quants, hue='Quantile', palette=['#B4C7E6', '#FCE799', '#C6E0B4', '#F9CCAD'], inner='quartile', cut=0, width=0.5, legend=False)

# Lines for real values
for i in range(4):
	lo, hi=plt.ylim()
	real=df.loc['Real', quants[i]]
	plt.plot([i-0.3, i+0.3], [real, real], color='k', lw=3)

	if real>hi:
		plt.ylim(lo, real*1.1)
		hi=real*1.1

	# p values
	pval=statdf[statdf.Quantile==quants[i]]['Empirical p value'].to_list()[0]
	plt.text(i, hi, 'p=%.3f' % pval, color='k', ha='center', va='bottom')

plt.ylim(lo, hi*1.05)

plt.ylabel('Number of genes in bin')
plt.title(name)

plt.savefig(OUTPUT_FIGURE)
