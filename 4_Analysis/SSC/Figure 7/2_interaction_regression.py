import pandas as pd
import numpy as np

import scipy.stats as stats
import statsmodels.formula.api as smf

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cmx
from matplotlib.backends.backend_pdf import PdfPages

import seaborn as sns

matplotlib.rcParams['pdf.fonttype']=42

import re

# Test for interactions between primary and seconday variants using linear regression

# Input and output files
SSC_DATA="/path/to/SSC/while/cohort/data.csv" # Use the output of script 3_Data preparation\SSC\3_merge_data.py

OUTPUT_STATS="/path/to/output/regression/statistics.csv" # Data presented in Table S6B
OUTPUT_FOREST="/path/to/output/forestplot.pdf" # Figure presented in Fig 7B
OUTPUT_FAN="/path/to/output/fanplot.pdf" # Figure presented in Fig S7E

# Define variables
quant_phenos=['Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Social responsiveness (SRS)', 'Repetitive behavior (RBS-R)', 'Coordination disorder (DCDQ)', 'BMI z-score']
primary_variants=['DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'Any variant']
vars=['All coding SNVs', 'Missense', 'LOF', 'Splice', 'Genes del.', 'Genes dup.', 'STRs',
		'Intelligence PRS', 'SCZ PRS', 'Education PRS',
		'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)']

replace_string='[ .()/-]'

# Helper function
def run_model(moddf, output_col, primary_variant, secondary_variant):
	# Replace problem characters in variable names
	form='%s ~ Age + Sex + %s + %s + %s*%s' % (re.sub(replace_string, '_', output_col), re.sub(replace_string, '_', primary_variant), re.sub(replace_string, '_', secondary_variant),
												re.sub(replace_string, '_', primary_variant), re.sub(replace_string, '_', secondary_variant))
	mod=smf.ols(formula=form, data=moddf)
    
	res=mod.fit()
	
	# Parse model
	num_vars=6
	ci=res.conf_int(alpha=0.05)
	ci.columns=['95% C.I. lower', '95% C.I. upper']
	
	r2=res.rsquared
	
	mod_name='%s ~ Age + Sex + %s + %s + %s*%s' % (output_col, primary_variant, secondary_variant, primary_variant, secondary_variant)
	
	res_dict={'Primary variant':[primary_variant]*num_vars, 'Secondary variant':[secondary_variant]*num_vars, 'Phenotype':[output_col]*num_vars, 'Model':[mod_name]*num_vars,
				'Variable':['Intercept', 'Age', 'Sex', primary_variant, secondary_variant, 'Interaction'], 'R2':[r2]*num_vars}
	mod_res=pd.DataFrame(res_dict)
	mod_res.index=['Intercept', 'Age', 'Sex', re.sub(replace_string, '_', primary_variant), re.sub(replace_string, '_', secondary_variant), re.sub(replace_string, '_', primary_variant)+':'+re.sub(replace_string, '_', secondary_variant)]
	
	mod_res['Estimate']=res.params
	mod_res['Error']=res.bse
	mod_res['p value']=res.pvalues
	mod_res=pd.merge(mod_res, ci, left_index=True, right_index=True)
	
	mod_res=mod_res[['Primary variant', 'Secondary variant', 'Phenotype', 'Model', 'Variable', 'Estimate', 'Error', '95% C.I. lower', '95% C.I. upper', 'p value', 'R2']]
	return mod_res

# Load data
df=pd.read_csv(SSC_DATA)
df=df[(~df[quant_phenos].isnull().all(axis=1)) & (~df.Sex.isnull()) & (~df.Age.isnull())]

# Center and scale numeric variables and convert others to binary
num_vars=vars+quant_phenos+['Age']
for nv in num_vars:
	df[nv]=(df[nv]-df[nv].mean(skipna=True))/df[nv].std(skipna=True)

bin_vars=['Sex']+primary_variants
df['Sex']=df.Sex.map({'M':1, 'F':0})
df[primary_variants]=df[primary_variants].astype(int)

# Run regressions
outdf=pd.DataFrame()
for v in vars:
	for pheno in quant_phenos:
		for pv in primary_variants:
			mtype='linear'
			input_cols=['Sex', 'Age', v, pv]
			moddf=df[~df[input_cols+[pheno]].isnull().any(axis=1)].copy()

			# Scale input cols
			for ic in input_cols+[pheno]:
				if ic!='Sex':
					moddf[ic]=(moddf[ic]-moddf[ic].mean())/moddf[ic].std()
			
			moddf=moddf[[pheno]+input_cols]
			# Replace problem characters in variable names
			cols=moddf.columns.to_list()
			moddf.columns=[re.sub(replace_string, '_', i) for i in cols]
			
			out=run_model(moddf, pheno, pv, v)
			out['Sample size']=moddf.shape[0]
			out.reset_index(inplace=True, drop=True)
			outdf=pd.concat([outdf, out])
outdf=outdf[['Primary variant', 'Secondary variant', 'Phenotype', 'Model', 'Sample size', 'Variable', 'Estimate', 'Error', '95% C.I. lower', '95% C.I. upper', 'p value', 'R2']]

# FDR correction
outdf['BH FDR']=stats.false_discovery_control(outdf['p value'], method='bh')

# Save results
outdf.to_csv(OUTPUT_STATS, index=False)

# Save a copy for forest plots
plotdf=outdf.copy()

# Plot results
outdf=outdf[outdf.Variable=='Interaction']

outdf['Secondary variant']=pd.Categorical(outdf["Secondary variant"], vars)
outdf['Primary variant']=pd.Categorical(outdf['Primary variant'], primary_variants)
outdf['Phenotype']=pd.Categorical(outdf['Phenotype'], quant_phenos)

outdf.sort_values(by=['Phenotype', 'Secondary variant', 'Primary variant'], inplace=True)
outdf.reset_index(drop=True, inplace=True)

# Add some variables for plotting
var_idx_map={}
for v in vars:
	var_idx_map[v]=vars.index(v)
outdf['variable_idx']=outdf['Secondary variant'].map(var_idx_map).astype(int)
pheno_idx_map={}
rev_pheno=quant_phenos.copy()
rev_pheno.reverse()
for p in quant_phenos:
	pheno_idx_map[p]=rev_pheno.index(p)
outdf['phenotype_idx']=outdf.Phenotype.map(pheno_idx_map).astype(int)

colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
biggest_est=outdf.Estimate.abs().max()
norm=matplotlib.colors.Normalize(vmin=-biggest_est, vmax=biggest_est)
outdf['estimate_color']=outdf.Estimate.apply(lambda x: cmap(norm(x)))

outdf['significance_color']='k'
outdf.loc[outdf['BH FDR']<=0.05, 'significance_color']='red'

outdf['neglogp']=-np.log10(outdf['p value'])
outdf['radius']=0.5*(outdf.neglogp)/(outdf.neglogp.max())

outdf['theta1']=outdf['Primary variant'].map({'DBD Tier 1 SNVs':90, 'Large rare deletions':180, 'Large rare duplications':0, 'Any variant':270}).astype(int)
outdf['theta2']=outdf.theta1+90

# Make interaction plot
fig, axs = plt.subplots(ncols=2, figsize=(20, 5), sharex=True, sharey=True)
for idx, row in outdf.iterrows():
	axs[0].add_artist(mpatches.Wedge((row['variable_idx'], row['phenotype_idx']), row['radius'], row['theta1'], row['theta2'], fc=row['estimate_color'], ec=row['significance_color'], lw=1))
axs[0].set_xlim(-1, 18)
axs[0].set_ylim(-1, 7)
axs[0].set_aspect('equal', adjustable='box')

axs[0].set_yticks([i for i in range(0, 7)], rev_pheno)
axs[0].set_xticks([i for i in range(0, 17)], vars, rotation=90)

# Make legend panel
legdf=pd.DataFrame({'pvalue':[0.1, 0.01, 0.001, 0.0001, 0.0001, 0.001, 0.001, 0.001, 0.001], 'theta1':[90, 90, 90, 90, 90, 90, 180, 0, 270], 'y':[3, 2, 1, 2, 1]+[4]*4, 'x':[2, 2, 2, 5, 5, 2, 2, 9, 9],
					'fillcolor':['#CCCCCC']*9, 'edgecolor':['k']*4+['red']+['k']*4,
					'text':['0.1', '0.01', '0.001', '0.0001', 'FDR<0.05', 'DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'Any variant']})
legdf['neglogp']=-np.log10(legdf.pvalue)
legdf['radius']=0.5*(legdf.neglogp)/(outdf.neglogp.max())
legdf['theta2']=legdf.theta1+90

# Adjustments
legdf.loc[legdf.theta1.isin([0, 270]), 'x']=legdf[legdf.theta1.isin([0, 270])].x-legdf[legdf.theta1.isin([0, 270])].radius

legdf['text_x']=legdf.x+0.1
legdf['text_y']=legdf.y+(0.5*legdf.radius)

# Text adjustments
legdf.loc[legdf.theta1.isin([180, 270]), 'text_y']=legdf[legdf.theta1.isin([180, 270])].y-(0.5*legdf[legdf.theta1.isin([180, 270])].radius)
legdf.loc[legdf.theta1.isin([0, 270]), 'text_x']=legdf[legdf.theta1.isin([0, 270])].x+0.1+legdf[legdf.theta1.isin([0, 270])].radius

for idx, row in legdf.iterrows():
	axs[1].add_artist(mpatches.Wedge((row['x'], row['y']), row['radius'], row['theta1'], row['theta2'], fc=row['fillcolor'], ec=row['edgecolor'], lw=1))
	# Add text
	axs[1].text(row['text_x'], row['text_y'], row['text'], color='k', va='center', ha='left')
axs[1].set_aspect('equal', adjustable='box')

# Colorbar
scalarMap=cmx.ScalarMappable(norm=norm, cmap=cmap)
scalarMap.set_array([])
steps=[-biggest_est]
for i in range(4):
	steps.append(steps[i]+(biggest_est/2.5))
steps.append(biggest_est)
plt.colorbar(scalarMap, ax=axs[1], ticks=steps, orientation='vertical')

plt.tight_layout()
plt.savefig(OUTPUT_FAN)
plt.close()

# Create forest plots for each comparison
pdf=PdfPages(OUTPUT_FOREST)
for pheno in quant_phenos:
	for v in vars:
		subdf=plotdf[(plotdf['Secondary variant']==v) & (plotdf.Phenotype==pheno) & ~(plotdf.Variable.isin(['Age', 'Sex', 'Intercept']))].copy()
		
		# Annotate stars
		subdf['star']=''
		subdf.loc[subdf['p value']<=0.05, 'star']='*'
		subdf.loc[subdf['BH FDR']<=0.05, 'star']='**'
		
		# Rename variables for easy plotting
		subdf['plotvar']=subdf.Variable
		subdf.loc[subdf.Variable==subdf['Primary variant'], 'plotvar']='Primary variant'
		subdf.loc[subdf.Variable==subdf['Secondary variant'], 'plotvar']=v
		
		# Add in spacer columns
		subdf=pd.concat([subdf, pd.DataFrame({'Primary variant':['Z']*3, 'plotvar':['Primary variant', v, 'Interaction']})])
		
		# Order values
		subdf['Primary variant']=pd.Categorical(subdf['Primary variant'], primary_variants+['Z'])
		subdf['plotvar']=pd.Categorical(subdf.plotvar, ['Primary variant', v, 'Interaction'])
		subdf.sort_values(by=['plotvar', 'Primary variant'], inplace=True, ascending=False)
		
		# Annotate order
		subdf.reset_index(inplace=True, drop=True)
		idx=subdf.index.to_list()
		idx.reverse()
		subdf['idx']=idx
		subdf=subdf[subdf['Primary variant']!='Z']
		subdf.reset_index(inplace=True, drop=True)
		
		fig, ax = plt.subplots(figsize=(6, 3))
		
		palette={'DBD Tier 1 SNVs':'#f69320', 'Large rare deletions':'#f4e46c', 'Large rare duplications':'#f37a7a', 'Any variant':'k'}
		
		# Points
		sns.scatterplot(data=subdf, x='idx', y='Estimate', hue='Primary variant', hue_order=primary_variants, palette=palette, edgecolor='None', legend=False)
		
		# Confidence intervals
		for index, row in subdf.iterrows():
			color=palette[row['Primary variant']]
			lo=row['95% C.I. lower']
			hi=row['95% C.I. upper']
			idx=row['idx']
			
			plt.plot([idx, idx], [lo, hi], color=color)
		
			# Also plot stars
			if hi>0:
				plt.text(idx, hi*1.2, row['star'], color='k', ha='center', va='center_baseline')
			else:
				plt.text(idx, lo*1.2, row['star'], color='k', ha='center', va='center_baseline')
		
		# Add a dotted line at 0
		ax.axhline(0, color='grey', ls='--')
		
		lo, hi = plt.ylim()
		plt.ylim(lo*1.25, hi*1.25)
		
		# Update x labels
		order=['Primary variant', v, 'Interaction']
		xs=[]
		for o in order:
			xs.append(subdf[subdf.plotvar==o].idx.mean())
		
		plt.xticks(xs, order)
		plt.xlabel('')
		
		plt.title(v+' '+pheno)
		plt.tight_layout()
		
		pdf.savefig()
		plt.close()

pdf.close()
