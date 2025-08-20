import pandas as pd

import scipy.stats as stats
import statsmodels.api as sm

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Run linear regressions on SSC proband phenotypes with two equations:
#	(B) Phenotype ~ All coding SNVs + Genes del. + Genes dup. + STRs + Sex + SCZ PRS
#	(C) Phenotype ~ All coding SNVs (LF) + Genes del. (LF) + Genes dup. (LF) + STRs (LF) + Sex + SCZ PRS

# Input and output files
DATA_DIR="/paht/to/SSC/data/directory/" # Use the output directory of script 3_Data preparation\SSC\2_compile_data.py

OUTPUT_STATS="/path/to/output/regression/statistics.csv" # Data presented in Table S6A
OUTPUT_FIG="/path/to/output/forestplot.pdf" # Figure presented in Fig S7C

# Identify variables
quant_phenos=['Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Social responsiveness (SRS)', 'Repetitive behavior (RBS-R)', 'Coordination disorder (DCDQ)', 'BMI z-score']
cohorts=['DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'No primary variant']
vars=['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs',
		'SCZ PRS', 
		'All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)']

# Helper function
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

# Run regressions
outdf=pd.DataFrame(columns=['Cohort', 'Phenotype', 'Variable', 'Model', 'Test', 'N', 'Estimate', 'Error', '95% C.I. lower', '95% C.I. upper', 'p value', 'R2'])
for cohort in cohorts:
	df=pd.read_csv(f'{DATA_DIR}/{cohort}.csv')
	df=df[(~df[quant_phenos].isnull().all(axis=1)) & (~df.Sex.isnull()) & (~df['SCZ PRS'].isnull()) & (~df['All coding SNVs'].isnull()) & (~df['STRs'].isnull()) & (~df['Genes del.'].isnull())]
	df=df[['Sample', 'Sex']+vars+quant_phenos]

	# Center and scale numeric variables and convert others to binary
	num_vars=vars+quant_phenos
	for nv in num_vars:
		df[nv]=(df[nv]-df[nv].mean(skipna=True))/df[nv].std(skipna=True)

	bin_vars=['Sex']
	df['Sex']=df.Sex.map({'M':1, 'F':0})

	model_inputs={'B':['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs', 'SCZ PRS', 'Sex'],
					'C':['All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)', 'SCZ PRS', 'Sex']}
	for m in model_inputs.keys():
		for pheno in quant_phenos:
			# Get model inputs
			mtype='linear'
			input_cols=model_inputs[m]
			moddf=df[~df[input_cols+[pheno]].isnull().any(axis=1)].copy()
			
			# Scale input cols
			for ic in input_cols:
				if ic!='Sex':
					moddf[ic]=(moddf[ic]-moddf[ic].mean())/moddf[ic].std()
			if mtype=='linear':
				moddf[pheno]=(moddf[pheno]-moddf[pheno].mean())/moddf[pheno].std()
			
			out=run_model(moddf, input_cols, pheno, m, type=mtype)
			out['Cohort']=cohort
			outdf=pd.concat([outdf, out])

# Perform FDR corrections on the joint variant models
outdf['BH FDR']=stats.false_discovery_control(outdf['p value'], method='bh')

# Update model annotations
outdf2=outdf.copy()
outdf2['Model']=outdf2.Model.map({'B':'Phenotype ~ Sex + All coding SNVs + Genes del. + Genes dup. + STRs + SCZ PRS',
									'C':'Phenotype ~ Sex + All coding SNVs (LF) + Genes del. (LF) + Genes dup. (LF) + STRs (LF) + SCZ PRS'})
outdf2['Test']=outdf2.Test.map({'logistic':'Logistic regression', 'linear':'Linear regression'})

# Save results
outdf2.to_csv(OUTPUT_STATS, index=False)

# Plot results as forest plots
outdf.reset_index(inplace=True, drop=True)

# Add some dummay variables for spacing
for m in ['B', 'C']:
	outdf=pd.concat([outdf, pd.DataFrame({'Variable':model_inputs[m]*2, 'Model':[m]*(2*len(model_inputs[m])), 'Cohort':['Z']*(2*len(model_inputs[m]))})])

# Remove variables we don't want to plot
outdf=outdf[(~outdf.Variable.isin(['Sex', 'Intercept']))]

# Organize data
outdf['Variable']=pd.Categorical(outdf.Variable, ['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs', 'All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)', 'SCZ PRS'])
outdf['Cohort']=pd.Categorical(outdf.Cohort, cohorts+['Z'])
outdf['star']=''
outdf.loc[outdf['p value']<=0.05, 'star']='*'

outdf.sort_values(by=['Model', 'Variable', 'Cohort'], inplace=True, ascending=[True, False, True])

pdf=PdfPages(OUTPUT_FIG)
for mod in ['B', 'C']:
	for pheno in quant_phenos:
		fig, ax = plt.subplots(figsize=(4, 6))
		subdf=outdf[((outdf.Phenotype==pheno) | (outdf.Phenotype.isnull())) & (outdf.Model==mod)].copy()
		subdf.reset_index(inplace=True, drop=True)
		
		subdf['idx']=subdf.index.to_list()
		subdf=subdf[subdf.Cohort!='Z']
		subdf.reset_index(inplace=True, drop=True)
		
		palette={'DBD Tier 1 SNVs':'#f69320', 'Large rare deletions':'#f4e46c', 'Large rare duplications':'#f37a7a', 'No primary variant':'#942911'}
		
		# Points
		sns.scatterplot(data=subdf, y='idx', x='Estimate', hue='Cohort', hue_order=cohorts, palette=palette, edgecolor='None', legend=False)
		
		# Confidence intervals
		for index, row in subdf.iterrows():
			color=palette[row['Cohort']]
			lo=row['95% C.I. lower']
			hi=row['95% C.I. upper']
			idx=row['idx']
			
			plt.plot([lo, hi], [idx, idx], color=color)
		
			# Also plot stars
			if hi>0:
				plt.text(hi*1.2, idx, row['star'], color='k', ha='left', va='center_baseline')
			else:
				plt.text(lo*1.2, idx, row['star'], color='k', ha='left', va='center_baseline')
		
		# Add a dotted line at 0
		lo, hi = plt.ylim()
		plt.plot([0, 0], [lo-2, hi+2], color='k', ls=':', zorder=-1)
		plt.ylim(lo, hi)
		
		lo, hi = plt.xlim()
		plt.xlim(lo*1.25, hi*1.25)
		
		# Change the x labels
		locs=[1.5, 7.5, 13.5, 19.5, 25.5]
		if mod=='B':
			plt.yticks(locs, reversed(['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs', 'SCZ PRS']))
		if mod=='C':
			plt.yticks(locs, reversed(['All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)', 'SCZ PRS']))
		
		plt.ylabel('')
		
		plt.title(mod+' '+pheno)
		plt.tight_layout()
		
		pdf.savefig()
		plt.close()

pdf.close()