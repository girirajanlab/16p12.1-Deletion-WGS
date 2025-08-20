import pandas as pd

from pedigree_functions import *

import seaborn as sns

# Plot pedigrees and burden for all multigenerational pedigrees

# Input and output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
FAM="/path/to/sample/FAM/file.fam"

OUTPUT_FIGURE="/path/to/output/figures.pdf"

# Load data
df=pd.read_csv(TABS1A)

gen3=df[df.Relationship.isin(['Grandmother', 'Grandfather'])].Family.to_list()
df=df[df.Family.isin(gen3)]

vars=['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs']
child_domains=['ID/DD (Child domain)', 'Behavioral features (Child domain)', 'Psychiatric features (Child domain)', 'Nervous System Abnormalities (Child domain)', 'Congenital Anomalies (Child domain)', 'Growth/Skeletal Defects (Child domain)']
adult_domains=['Anxiety (Questionnaire)', 'Sleep trouble (Questionnaire)', 'Psychosis (Questionnaire)', 'Addiction (Questionnaire)', 'Mood lability (Questionnaire)']

domain_colors={'ID/DD (Child domain)':'#6A8A82', 'Behavioral features (Child domain)':'#9385BF', 'Psychiatric features (Child domain)':'#F9B163',
				'Nervous System Abnormalities (Child domain)':'#81B1D2', 'Congenital Anomalies (Child domain)':'#FFCE00', 'Growth/Skeletal Defects (Child domain)':'#EF7F71',
				'Anxiety (Questionnaire)':'#50723C', 'Sleep trouble (Questionnaire)': '#3D7EA6', 'Psychosis (Questionnaire)':'#AB6C82', 'Addiction (Questionnaire)':'#C34A36', 'Mood lability (Questionnaire)':'#D88C1C'}

df=df[['Sample', 'Family', 'Relationship', 'Sex', '16p12.1 deletion']+vars+child_domains+adult_domains]

# Load in fam file to determine relationships
fam=pd.read_csv(FAM, sep='\t')
fam['Sample']=fam.IID
fam['Family']=fam.FID

df=pd.merge(df, fam[['Sample', 'Family', 'Mother', 'Father']], on=['Sample', 'Family'], how='outer')

# Get all possible parents and their partners
def plot_pedigrees(subdf, pdf):
	subdf['generations']=generations(subdf)
	subdf['y']=y_position(subdf.generations)
	subdf['x']=x_position(subdf, skip_symbol='UNRECRUITED')

	fig=plt.figure(figsize=(10,6))
	ax1=plt.subplot2grid((3, 3), (1, 0), rowspan=2, colspan=2)
	# Make pedigree
	for idx, row in subdf.iterrows():
		if 'UNRECRUITED' in row.Sample:
			continue
		ec='k'
		lw=1
		if row['16p12.1 deletion']=='Carrier':
			ec='red'
			lw=2
		plot_shape(ax1, row.Sex, row.x, row.y, ec=ec, lw=lw)
		if 'UNRECRUITED' not in row.Mother or 'UNRECRUITED' not in row.Father:
			add_lines_up(ax1, row.x, row.y)
		
		# Add patches for phenotypes
		if row.Relationship in ['Sibling', 'Proband', 'Cousin']:
			phenos=row[child_domains].index[row[child_domains]>0].to_list()
			print(row.Sample, row.Sample, phenos, row[child_domains])
			if row[child_domains].isnull().all():
				plot_shape(ax1, row.Sex, row.x, row.y, fc='lightgrey', ec='None', zo=-1)
		else:
			phenos=row[adult_domains].index[row[adult_domains]>0].to_list()
			if row[adult_domains].isnull().all():
				plot_shape(ax1, row.Sex, row.x, row.y, fc='lightgrey', ec='None', zo=-1)

		if len(phenos)>0:
				plot_phenos(ax1, row.Sex, row.x, row.y, phenos, domain_colors)	
	add_connecting_lines(ax1, subdf, skip_symbol='UNRECRUITED')
	label_individuals(ax1, subdf, iid_col='Sample')
	ax1.axis('off')

	ax1.set_xlim(subdf.x.min()-0.2, subdf.x.max()+0.2)
	ax1.set_ylim(subdf.y.min()-0.2, subdf.y.max()+0.2)
	ax1.set_aspect('equal', adjustable='box')

	# Add custom legend to pedigrees
	ax2=plt.subplot2grid((3, 3), (0, 0), colspan=2)
	cl_dict={}
	for p in child_domains+adult_domains:
		p_short=p.split(' (')[0]
		cl_dict[p_short]=domain_colors[p]
	custom_pheno_legend(ax2, cl_dict)
	ax2.axis('off')

	# Plot change in burden over generations
	ax3=plt.subplot2grid((2, 3), (1, 2))
	burddf=subdf[(subdf.Relationship.isin(['Proband', 'Mother', 'Father', 'Grandmother', 'Grandfather'])) & (subdf['16p12.1 deletion']=='Carrier')][['Sample', 'generations']+vars].copy()
	burddf=burddf.melt(id_vars=['Sample', 'generations'], var_name='Variant', value_name='Variant burden')

	var_palette={'All coding SNVs':'k', 'STRs':'#94B8A3', 'Genes del.':'#D52C28', 'Genes dup.':'#9268AD'}

	sns.lineplot(data=burddf, x='generations', y='Variant burden', hue='Variant', palette=var_palette, ax=ax3, legend=False)
	sns.scatterplot(data=burddf, x='generations', y='Variant burden', hue='Variant', palette=var_palette, ax=ax3, legend=False)
	ax3.set_xticks([2, 1, 0], ['Proband', 'Carrier\nparent', 'Carrier\nGrandparent'])
	ax3.set_xlabel('')

	# Add a custom legend to the linegraph
	ax4=plt.subplot2grid((2, 3), (0, 2))
	custom_burden_legend(ax4, var_palette)
	ax4.axis('off')

	pdf.savefig()
	plt.close()

pdf=PdfPages(OUTPUT_FIGURE)
fams=sorted(list(df.Family.unique()))
for f in fams:
	plot_pedigrees(df[df.Family==f].copy(), pdf)
pdf.close()
