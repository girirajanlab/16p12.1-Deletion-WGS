import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Plot the significant results from meta analysis

# Input and output files
DD_UKB='/path/to/DD/UKB/meta.csv' # Use the outputs of script 2_meta_analysis.R
UKB_MYCODE_AOU='/path/to/UKB/MyCode/AoU/meta.csv'
DD_SPARK='/path/to/DD/SPARK/meta.csv'

OUTPUT_FIG="/path/to/output/figure.pdf" # Output figures presented in Figs 6G, S6H-I

# Parse results
vars=['All coding SNVs', 'Missense', 'LOF', 'Splice',
		'Genes del.', 'Genes dup.',
		'SCZ PRS', 'Education PRS', 'Autism PRS', 'Intelligence PRS',
		'All coding SNVs (LF)', 'Missense (LF)', 'LOF (LF)', 'Splice (LF)',
		'Genes del. (LF)', 'Genes dup. (LF)']

pdf=PdfPages(OUTPUT_FIG)
for f in [DD_UKB, UKB_MYCODE_AOU, DD_SPARK]:
	df=pd.read_csv(f)
	df=df[df.p_value<=0.05]

	plotdf=df.melt(id_vars=['Phenotype', 'Variant', 'p_value', 'FDR'])
	plotdf['Cohort']=plotdf.variable.str.split('_', expand=True)[0]
	plotdf['Metric']=[a.replace(b+'_', '').strip() for a, b in zip(plotdf['variable'], plotdf['Cohort'])]

	plotdf=pd.pivot(plotdf, index=['Phenotype', 'Variant', 'p_value', 'FDR', 'Cohort'], columns='Metric', values='value')
	plotdf.reset_index(inplace=True)

	order=['DD (children)', 'SPARK', 'DD (adults)', 'UKB', 'MyCode', 'AoU', 'Random']
	plotdf['Cohort']=pd.Categorical(plotdf.Cohort, order)
	plotdf.Variant=pd.Categorical(plotdf.Variant, vars)

	plotdf.sort_values(by=['Phenotype', 'Variant', 'Cohort'], inplace=True)
	plotdf.reset_index(inplace=True, drop=True)

	plotdf=plotdf[['Phenotype', 'Variant', 'Cohort', 'effect', 'CI_low', 'CI_upp', 'p_value', 'FDR']]
	ys=plotdf.index.to_list()
	plotdf['y']=ys

	plotdf['pv']=plotdf.Phenotype.astype(str)+'.'+plotdf.Variant.astype(str)
	pvs=list(plotdf.pv.unique())
	plotdf['y']=plotdf.y+plotdf.pv.apply(lambda x: pvs.index(x)+1)
	ys=plotdf.y.to_list()
	ys.reverse()
	plotdf['y']=ys

	palette={'DD (children)':'#254061', 'SPARK':'#b12422', 'Random':'k',
				'DD (adults)':'#B9CDE5', 'MyCode':'#87ba44', 'UKB':'#45235F', 'AoU':'#E3A857'}
	if f==DD_UKB:
		palette['UKB']='#b084cc'
	fig, axs=plt.subplots(ncols=3, gridspec_kw={'width_ratios': [1, 1, 1.5]}, sharey=True)
	ax=axs[2]
	sns.scatterplot(data=plotdf, x='effect', y='y', hue='Cohort', palette=palette, ax=ax, legend=False)
	for cohort in list(plotdf.Cohort.unique()):
		subdf=plotdf[plotdf.Cohort==cohort]
		ax.plot([subdf.CI_low, subdf.CI_upp], [subdf.y, subdf.y], color=palette[cohort])

	# Formatting
	ax.axvline(0, ls=':', color='darkgrey', zorder=0)
	ax.tick_params(left=False)
	ax.set_xlabel("Cohen's D")

	# Add row annotations
	pvdf=plotdf[['Phenotype', 'Variant', 'y']].groupby(['Phenotype', 'Variant']).agg('mean')
	pvdf=pvdf[~pvdf.y.isnull()]
	pvdf.reset_index(inplace=True)

	for idx, row in pvdf.iterrows():
		axs[1].text(0, row['y'], row['Variant'], color='k', ha='left', va='center')
	axs[1].set_xlim([-0.1, 3])
	axs[1].axis('off')

	phdf=plotdf[['Phenotype', 'y']].groupby(['Phenotype', ]).agg('mean')
	phdf=phdf[~phdf.y.isnull()]
	phdf.reset_index(inplace=True)

	for idx, row in phdf.iterrows():
		axs[0].text(0, row['y'], row['Phenotype'], color='k', ha='left', va='center')
	axs[0].set_xlim([-0.1, 3])
	axs[0].axis('off')

	phenos=sorted(list(set(plotdf.Phenotype.to_list())))
	for p1 in phenos:
		for p2 in phenos:
			if phenos.index(p1)!=phenos.index(p2)+1:
				continue
			pmean=(pvdf[pvdf.Phenotype==p1].y.max()+pvdf[pvdf.Phenotype==p2].y.min())/2
			for i in range(3):
				axs[i].axhline(pmean, color='darkgrey', zorder=0)

	# Save
	pdf.savefig()
	plt.close()
pdf.close()