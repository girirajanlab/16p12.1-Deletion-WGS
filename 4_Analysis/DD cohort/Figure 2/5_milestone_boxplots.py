import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import seaborn as sns

import scipy.stats as stats

matplotlib.rcParams['pdf.fonttype']=42

# Plot developmental milestones for probands, carrier siblings, and non-carrier siblings

# Input and Output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
OUTPUT_FIG="/path/to/output/figure.pdf"
OUTPUT_STATS="/paht/to/output/statistics/file.csv"

# Load data
dat=pd.read_csv(TABS1A)
dat=dat[dat['Estonian Biobank Sample']!='X']

# Subset to needed data
milestones=['Age first smiled (months)', 'Age first laughed (months)', 'Age first rolled over (months)', 'Age first sat up without assistance (months)', 'Age first crawled (months)', 'Age able to pull self up to stand (months)',
				'Age able to stand alone (months)', 'Age took first steps (months)', 'Age spoke first words (months)', 'Age first walked alone (months)', 'Age first spoke 2-3 word sentences (months)', 'Age spoke complete sentences (months)']
dat=dat[['Sample', 'Relationship', '16p12.1 deletion']+milestones]
dat=dat[dat.Relationship.isin(['Proband', 'Sibling', 'Cousin'])]
dat=dat[dat['16p12.1 deletion'].isin(['Carrier', 'Noncarrier'])]

# Assign samples into Groups
dat['Group']='Proband'
dat.loc[(dat.Relationship.isin(['Sibling', 'Cousin'])) & (dat['16p12.1 deletion']=='Carrier'), 'Group']='Carrier sibling/cousin'
dat.loc[(dat.Relationship.isin(['Sibling', 'Cousin'])) & (dat['16p12.1 deletion']=='Noncarrier'), 'Group']='Noncarrier sibling/cousin'

dat=dat[(dat.Group.isin(['Proband', 'Carrier sibling/cousin', 'Noncarrier sibling/cousin'])) & (~dat[milestones].isnull().all(axis=1))]

# Reformat to long
df=pd.DataFrame(columns=['Sample', 'Relationship', 'Group', 'Milestone', 'Age (months)'])
for m in milestones:
	subdf=dat[['Sample', 'Relationship', 'Group', m]].copy()
	subdf.columns=['Sample', 'Relationship', 'Group', 'Age (months)']
	subdf['Milestone']=m
	subdf=subdf[~subdf['Age (months)'].isnull()]
	df=pd.concat([df, subdf], ignore_index=True)

# Do t-tests comparing each Group to each other and to the CDC guidelines
Groups=['Proband', 'Carrier sibling/cousin', 'Noncarrier sibling/cousin']
stat_lst=[]
for m in milestones:
	for g1 in Groups:
		g1_vals=df[(df.Milestone==m) & (df.Group==g1)]['Age (months)'].to_list()
		for g2 in Groups:
			if Groups.index(g1)>=Groups.index(g2):
				continue
			
			g2_vals=df[(df.Milestone==m) & (df.Group==g2)]['Age (months)'].to_list()
			
			res=stats.ttest_ind(g1_vals, g2_vals, alternative='greater')
			
			stat_lst.append([g1, g2, m, 'one tailed t-test', len(g1_vals), len(g2_vals), sum(g1_vals)/len(g1_vals), sum(g2_vals)/len(g2_vals), res.statistic, res.pvalue])

stat_df=pd.DataFrame(stat_lst, columns=['Group 1', 'Group 2', 'Milestone', 'Test', 'Group 1 sample size', 'Group 2 sample size', 'Group 1 mean', 'Group 2 mean', 'statistic', 'p value'])
stat_df['BH FDR']=np.nan
stat_df.loc[~stat_df['p value'].isnull(), 'BH FDR']=stats.false_discovery_control(stat_df[~stat_df['p value'].isnull()]['p value'].to_numpy(), method='bh')

# Save
stat_df.to_csv(OUTPUT_STATS, index=False)

# Add stars for plotting
stat_df['star']='n.s.'
stat_df.loc[stat_df['p value']<=0.05, 'star']='*'
stat_df.loc[stat_df['BH FDR']<=0.05, 'star']='**'

# Plot
fig, ax=plt.subplots(figsize=(10, 5))

sns.boxplot(data=df, y='Milestone', x='Age (months)', hue='Group', hue_order=['Noncarrier sibling/cousin', 'Carrier sibling/cousin', 'Proband'], fliersize=0, palette=['#98B6B1', '#FF6663', '#006C67'], zorder=-1, linewidth=0.75)

color_dict={'Noncarrier sibling/cousin':'#98B6B1', 'Carrier sibling/cousin':'#FF6663', 'Proband':'#006C67'}

for i, m in enumerate(milestones):
	# Plot significance as stars
	max_val=df[df.Milestone==m]['Age (months)'].max()
	
	# Add stars
	offset=1
	for g1 in Groups:
		for g2 in Groups:
			if Groups.index(g1)>=Groups.index(g2):
				continue
			
			start_off=0.05
			if g2=='Noncarrier sibling/cousin':
				start_off=-0.35
			end_off=-0.05
			if g1=='Proband':
				end_off=0.35
			
			text_off=0
			if g1=='Carrier sibling/cousin':
				text_off=-0.2
			if g2=='Carrier sibling/cousin':
				text_off=0.2
		
			txt=stat_df[(stat_df['Group 1']==g1) & (stat_df['Group 2']==g2) & (stat_df.Milestone==m)]['star'].to_list()[0]
		
			if txt!='n.s.':
				plt.plot([max_val+offset, max_val+offset], [i+start_off, i+end_off], color='k', linewidth=0.75)
				plt.text(max_val+offset+0.2, i+text_off, txt, fontsize=4, ha='left', va='center')
			
			offset+=1.2
	
plt.tight_layout()
plt.savefig(OUTPUT_FIG)