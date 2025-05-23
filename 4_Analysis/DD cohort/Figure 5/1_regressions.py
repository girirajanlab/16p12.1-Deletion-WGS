import pandas as pd
import numpy as np

import scipy.stats as stats
import statsmodels.api as sm

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Run logistic and linear regressions on 16p12.1 deletion proband phenotypes with three equations:
#	(A) Phenotype ~ Rare coding variants + SCZ PRS + Sex
#	(B) Phenotype ~ All coding SNVs + Genes del. + Genes dup. + STRs + Sex + SCZ PRS
#	(C) Phenotype ~ All coding SNVs (LF) + Genes del. (LF) + Genes dup. (LF) + STRs (LF) + Sex + SCZ PRS

# Input and output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
OUTPUT_STATS="/path/to/output/statistics/file.csv" # These are the statistics presented in Table S4A
OUTPUT_FOREST="/path/to/output/forest/plot/figures.pdf" # Forest plots presented in Fig 4A-B and S4A
OUTPUT_HEATMAP="/path/to/output/heatmap/figure.pdf" # Heatmap presented in Fig S4B

# Define phenotypes and variants
child_domains=['Behavioral features (Child domain)', 'Psychiatric features (Child domain)', 'Nervous System Abnormalities (Child domain)', 'Congenital Anomalies (Child domain)', 'Growth/Skeletal Defects (Child domain)']
quant_phenos=['De Vries Score', 'BMI Z Score', 'Head Circumference Z Score', 'HRS-MAT', 'SRS Raw Score']
vars=['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs',
		'SCZ PRS', 
		'All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)']

# Load data
df=pd.read_csv(TABS1A)
df=df[df['Estonian Biobank Sample']!='X']
df=df[df.Relationship=='Proband']
df=df[((~df[child_domains].isnull().all(axis=1)) | (~df[quant_phenos].isnull().all(axis=1))) & (~df.Sex.isnull()) & (~df['SCZ PRS'].isnull()) & (~df['All coding SNVs'].isnull()) & (~df['STRs'].isnull()) & (~df['Genes del.'].isnull())]
df=df[['Sample', 'Sex']+vars+child_domains+quant_phenos]

# Center and scale numeric variables and convert others to binary
df['Rare coding variants']=df[['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs']].sum(axis=1)

num_vars=['Rare coding variants']+vars+quant_phenos
for nv in num_vars:
	df[nv]=(df[nv]-df[nv].mean(skipna=True))/df[nv].std(skipna=True)

bin_vars=['Sex']+child_domains
df['Sex']=df.Sex.map({'M':1, 'F':0})
# For the phenotypic domains, split at a value to create approximately even groups
splits={'Behavioral features (Child domain)':1, 'Psychiatric features (Child domain)':0, 'Nervous System Abnormalities (Child domain)':0, 'Congenital Anomalies (Child domain)':1, 'Growth/Skeletal Defects (Child domain)':2}
for key in splits.keys():
	df.loc[df[key]<=splits[key], key]=0
	df.loc[df[key]>splits[key], key]=1

def run_model(moddf, input_vars, output_col, model_name, type='logistic'):
	X=sm.add_constant(moddf[input_vars].to_numpy())
	if type=='logistic':
		mod=sm.Logit(moddf[output_col].to_numpy(), X)
	else:
		mod=sm.OLS(moddf[output_col].to_numpy(), X)
    
	res=mod.fit()
	
	# Parse model
	num_vars=len(input_vars)+1
	ci=res.conf_int(alpha=0.05)
	
	if type=='logistic':
		r2=res.prsquared
	else:
		r2=res.rsquared
	
	res_dict={'Phenotype':[output_col]*num_vars, 'Variable':['Intercept']+input_vars, 'Model':[model_name]*num_vars, 'Test':[type]*num_vars, 'N':[moddf.shape[0]]*num_vars,
				'Estimate':res.params, 'Error':res.bse, '95% C.I. lower':[i[0] for i in ci], '95% C.I. upper':[i[1] for i in ci], 'p value':res.pvalues, 'R2':[r2]*num_vars}
	mod_res=pd.DataFrame(res_dict)
	
	return mod_res

model_inputs={'Phenotype ~ Sex + Rare coding variants + SCZ PRS':['Rare coding variants', 'SCZ PRS', 'Sex'],
				'Phenotype ~ Sex + SCZ PRS + All coding SNVs. + Genes del. + Genes dup. + STRs':['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs', 'SCZ PRS', 'Sex'],
				'Phenotype ~ Sex + SCZ PRS + All coding SNVs. (LF) + Genes del. (LF) + Genes dup. (LF) + STRs (LF)':['All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)', 'SCZ PRS', 'Sex']}
outdf=pd.DataFrame(columns=['Phenotype', 'Variable', 'Model', 'Test', 'N', 'Estimate', 'Error', '95% C.I. lower', '95% C.I. upper', 'p value', 'R2'])
for m in model_inputs.keys():
	for pheno in child_domains+quant_phenos:
		# Get model inputs
		mtype='linear'
		if pheno in child_domains:
			mtype='logistic'
		input_cols=model_inputs[m]
		moddf=df[~df[input_cols+[pheno]].isnull().any(axis=1)].copy()
		
		# Scale input cols
		for ic in input_cols:
			if ic!='Sex':
				moddf[ic]=(moddf[ic]-moddf[ic].mean())/moddf[ic].std()
		if mtype=='linear':
			moddf[pheno]=(moddf[pheno]-moddf[pheno].mean())/moddf[pheno].std()
		
		out=run_model(moddf, input_cols, pheno, m, type=mtype)
		outdf=pd.concat([outdf, out])

# Perform FDR corrections on the joint variant models
outdf['BH FDR']=np.nan
for pg in [child_domains, quant_phenos]:
	outdf.loc[outdf.Phenotype.isin(pg), 'BH FDR']=stats.false_discovery_control(outdf[outdf.Phenotype.isin(pg)]['p value'], method='bh')

for v in ['Rare coding variants']+vars:
	for pheno in child_domains+quant_phenos:
		# Get model inputs
		mtype='linear'
		if pheno in child_domains:
			mtype='logistic'
		input_cols=[v, 'Sex']
		moddf=df[~df[input_cols+[pheno]].isnull().any(axis=1)].copy()
		
		modname='Phenotype ~ Sex + '+v
		
		# Scale input cols
		for ic in input_cols:
			if ic!='Sex':
				moddf[ic]=(moddf[ic]-moddf[ic].mean())/moddf[ic].std()
		if mtype=='linear':
			moddf[pheno]=(moddf[pheno]-moddf[pheno].mean())/moddf[pheno].std()
		
		out=run_model(moddf, input_cols, pheno, modname, type=mtype)
		outdf=pd.concat([outdf, out])

# Save results
outdf.Phenotype=outdf.Phenotype.str.split(" \(", expand=True)[0]
outdf.Test=outdf.Test.map({'logistic':'Logistic regression', 'linear':'Linear regression'})
outdf.to_csv(OUTPUT_STATS, index=False)

# Plot results as forest plots
outdf.reset_index(inplace=True, drop=True)

model_inputs={'Phenotype ~ Sex + Rare coding variants + SCZ PRS':['Rare coding variants', 'SCZ PRS', 'Sex'],
				'Phenotype ~ Sex + SCZ PRS + All coding SNVs. + Genes del. + Genes dup. + STRs':['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs', 'SCZ PRS', 'Sex'],
				'Phenotype ~ Sex + SCZ PRS + All coding SNVs. (LF) + Genes del. (LF) + Genes dup. (LF) + STRs (LF)':['All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)', 'SCZ PRS', 'Sex']}

# Add some dummay variables for spacing
mods=['Phenotype ~ Sex + Rare coding variants + SCZ PRS', 'Phenotype ~ Sex + SCZ PRS + All coding SNVs. + Genes del. + Genes dup. + STRs', 'Phenotype ~ Sex + SCZ PRS + All coding SNVs. (LF) + Genes del. (LF) + Genes dup. (LF) + STRs (LF)']
for m in mods:
	outdf=pd.concat([outdf, pd.DataFrame({'Variable':model_inputs[m]*2, 'Model':[m]*(2*len(model_inputs[m])), 'Phenotype':['Z']*(2*len(model_inputs[m]))})])

# Remove variables we don't want to plot
vardf=outdf[outdf.Variable=='Intercept']
outdf=outdf[(~outdf.Variable.isin(['Sex', 'Intercept'])) & (outdf.Model.isin(mods))]
outdf=outdf[(outdf.Variable!='SCZ PRS') | (outdf.Model=='Phenotype ~ Sex + Rare coding variants + SCZ PRS')]

# Organize data
outdf['Variable']=pd.Categorical(outdf.Variable, ['Rare coding variants', 'SCZ PRS', 'All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs', 'All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)'])
outdf['Phenotype']=pd.Categorical(outdf.Phenotype, [i.split(" (")[0] for i in child_domains]+quant_phenos+['Z'])
outdf['star']=''
outdf.loc[outdf['p value']<=0.05, 'star']='*'

outdf.sort_values(by=['Model', 'Variable', 'Phenotype'], inplace=True)
print(outdf)

pdf=PdfPages(OUTPUT_FOREST)
for mod in mods:
	for pgroup in [[i.split(" (")[0] for i in child_domains], quant_phenos, ['De Vries Score', 'BMI Z Score', 'Head Circumference Z Score']]:
		if mod=='Phenotype ~ Sex + Rare coding variants + SCZ PRS':
			fig, ax = plt.subplots(figsize=(4, 3))
		else:
			fig, ax = plt.subplots(figsize=(6, 3))
		subdf=outdf[((outdf.Phenotype.isin(pgroup)) | (outdf.Phenotype=='Z')) & (outdf.Model==mod)].copy()
		subdf.reset_index(inplace=True, drop=True)
		
		subdf['idx']=subdf.index.to_list()
		subdf=subdf[subdf.Phenotype!='Z']
		subdf.reset_index(inplace=True, drop=True)
		
		if 'Behavioral features' in pgroup:
			palette={'Behavioral features':'#9385BF', 'Psychiatric features':'#F9B163', 'Nervous System Abnormalities':'#81B1D2',
						'Congenital Anomalies':'#FFCE00', 'Growth/Skeletal Defects':'#EF7F71'}
		else:
			palette={'De Vries Score':'#283540', 'BMI Z Score':'#9BC4BC', 'Head Circumference Z Score':'#30618C', 'HRS-MAT':'#D3FFE9', 'SRS Raw Score':'#8DDBE0'}
		
		# Points
		sns.scatterplot(data=subdf, x='idx', y='Estimate', hue='Phenotype', hue_order=pgroup, palette=palette, edgecolor='None', legend=False)
		
		# Confidence intervals
		for index, row in subdf.iterrows():
			color=palette[row['Phenotype']]
			lo=row['95% C.I. lower']
			hi=row['95% C.I. upper']
			idx=row['idx']
			
			plt.plot([idx, idx], [lo, hi], color=color)
		
			# Also plot stars
			plt.text(idx, hi*1.1, row['star'], color='k', ha='center', va='center')
		
		# Add a dotted line at 0
		lo, hi = plt.xlim()
		plt.plot([lo-2, hi+2], [0, 0], color='k', ls=':', zorder=-1)
		plt.xlim(lo, hi)
		
		lo, hi = plt.ylim()
		plt.ylim(lo, hi*1.2)
		
		# Change the x labels
		num_phenos=len(pgroup)
		if mod=='Phenotype ~ Sex + Rare coding variants + SCZ PRS':
			if num_phenos==5:
				locs=[2, 9]
			elif num_phenos==3:
				locs=[1, 6]
			plt.xticks(locs, ['Rare coding variants', 'SCZ PRS'])
		if mod=='Phenotype ~ Sex + SCZ PRS + All coding SNVs. + Genes del. + Genes dup. + STRs':
			if num_phenos==5:
				locs=[2, 9, 16, 23]
			elif num_phenos==3:
				locs=[1, 6, 11, 16]
			plt.xticks(locs, ['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs'])
		if mod=='Phenotype ~ Sex + SCZ PRS + All coding SNVs. (LF) + Genes del. (LF) + Genes dup. (LF) + STRs (LF)':
			if num_phenos==5:
				locs=[2, 9, 16, 23]
			elif num_phenos==3:
				locs=[1, 6, 11, 16]
			plt.xticks(locs, ['All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)'])
		
		plt.xlabel('')
		
		plt.title(mod)
		plt.tight_layout()
		
		pdf.savefig()
		plt.close()

pdf.close()

# Plot variance heatmap
# Remove phenotypes we do not want to plot
pheno_order=[i.split(" (")[0] for i in child_domains]
vardf=vardf[vardf.Phenotype.isin(pheno_order)]

# Organize data
mod_order=[mods[0], 'Phenotype ~ Sex + Rare coding variants', 'Phenotype ~ Sex + SCZ PRS',
			mods[1], 'Phenotype ~ Sex + All coding SNVs', 'Phenotype ~ Sex + Genes del.', 'Phenotype ~ Sex + Genes dup.', 'Phenotype ~ Sex + STRs',
			mods[2], 'Phenotype ~ Sex + All coding SNVs (LF)', 'Phenotype ~ Sex + Genes del. (LF)', 'Phenotype ~ Sex + Genes dup. (LF)', 'Phenotype ~ Sex + STRs (LF)']
vardf['Model']=pd.Categorical(vardf.Model, mod_order)
vardf['Phenotype']=pd.Categorical(vardf.Phenotype, pheno_order)

vardf.sort_values(by=['Model', 'Phenotype'], inplace=True)

# Convert to wide form
plotdf=pd.pivot(vardf, index='Model', columns='Phenotype', values='R2')
plotdf=plotdf.loc[mod_order, pheno_order]

# Draw heatmap
colors=["#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('Reds', colors, N=20)
sns.heatmap(data=plotdf, cmap=cmap, vmin=0, vmax=0.14, square=True, fmt='', linecolor='k', linewidths=0.75)

plt.tight_layout()
plt.savefig(OUTPUT_HEATMAP)
	