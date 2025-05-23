import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors

import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

# Create bar graphs for phenotypes in 16p12.1 deletion families from the DD cohort

# Input and Output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
OUTPUT_FIG="/path/to/output/figure.pdf"

# Subset necessary data
df=pd.read_csv(TABS1A)
df=df[df['Estonian Biobank Sample']!='X']
df=df[df.Relationship.isin(['Proband', 'Sibling', 'Cousin', 'Mother', 'Father'])]
df=df[df['16p12.1 deletion'].isin(['Carrier', 'Noncarrier'])]

# Assign samples into groups
df['Group']='Proband'
df.loc[(df.Relationship.isin(['Sibling', 'Cousin'])) & (df['16p12.1 deletion']=='Carrier'), 'Group']='Carrier sibling/cousin'
df.loc[(df.Relationship.isin(['Sibling', 'Cousin'])) & (df['16p12.1 deletion']=='Noncarrier'), 'Group']='Noncarrier sibling/cousin'
df.loc[(df.Relationship.isin(['Mother', 'Father'])) & (df['16p12.1 deletion']=='Carrier'), 'Group']='Carrier parent'
df.loc[(df.Relationship.isin(['Mother', 'Father'])) & (df['16p12.1 deletion']=='Noncarrier'), 'Group']='Noncarrier parent'

# Restrict to samples with phenotype information
adult_domains=['Depression (Questionnaire)', 'Anxiety (Questionnaire)', 'Sleep trouble (Questionnaire)', 'Psychosis (Questionnaire)', 'Addiction (Questionnaire)', 'Mood lability (Questionnaire)']
child_domains=['ID/DD (Child domain)', 'Behavioral features (Child domain)', 'Psychiatric features (Child domain)', 'Nervous System Abnormalities (Child domain)', 'Congenital Anomalies (Child domain)', 'Growth/Skeletal Defects (Child domain)']

df=df[(((df.Group.isin(['Carrier parent', 'Noncarrier parent'])) & (~df[adult_domains].isnull().all(axis=1)))) | (((df.Group.isin(['Proband', 'Carrier sibling/cousin', 'Noncarrier sibling/cousin'])) & (~df[child_domains].isnull().all(axis=1))))]
df=df[['Sample', 'Relationship', '16p12.1 deletion','Group']+adult_domains+child_domains]

# Create stacked barplots
pdf=PdfPages(OUTPUT_FIG)

# Child barplots
pro_colors = [mcolors.rgb2hex(mcolors.LinearSegmentedColormap.from_list('Proband', ['#ffffff', '#006C67'], 7)(i)) for i in range(1, 7)]
sc_colors = [mcolors.rgb2hex(mcolors.LinearSegmentedColormap.from_list('Carrier sibling/cousin', ['#ffffff', '#FF6663'], 7)(i)) for i in range(1, 7)]
snc_colors = [mcolors.rgb2hex(mcolors.LinearSegmentedColormap.from_list('Noncarrier sibling/cousin', ['#ffffff', '#3A95B3'], 7)(i)) for i in range(1, 7)]
	
legends=[]
for group in ['Noncarrier sibling/cousin', 'Carrier sibling/cousin', 'Proband']:
	last=[0]*len(child_domains)
	for score in range(0, 6):
		vals = []
		szs=[]
		for pheno in child_domains:
			num = df[(df[pheno]==score) & (df.Group==group)].shape[0]*100 / df[(~df[pheno].isnull()) & (df.Group==group)].shape[0]
			vals.append(num)
			sz=df[(~df[pheno].isnull()) & (df.Group==group)].shape[0]
			szs.append(sz)
		
		if group=='Proband':
			color=pro_colors[score]
			offset=-0.25
			offset2=0
		elif group=='Carrier sibling/cousin':
			color=sc_colors[score]
			offset=0
			offset2=2
		elif group=='Noncarrier sibling/cousin':
			color=snc_colors[score]
			offset=0.25
			offset2=0
		
		plt.bar([i+offset for i in range(len(child_domains))], vals, bottom=last, color=color, width=0.25)
		
		# Add sample sizes
		if score==0:
			for i in range(len(szs)):
				plt.text(x=i+offset, y=102+offset2, s='n='+str(szs[i]), ha='center', va='center', fontsize=6)
				
		newlast= [vals[i]+last[i] for i in range(len(last))]
		last=newlast
		
	plt.legend([i for i in range(0, 6)], bbox_to_anchor=(1, 1))
plt.xlabel("Phenotypic Domain")
plt.ylabel("Percentage of Group")
plt.title('Child Phenotypic Domains')
plt.xticks(ticks=[i for i in range(len(child_domains))], labels=[i.split(' (')[0] for i in child_domains], rotation=90)
plt.ylim(0, 110)

pdf.savefig()
plt.close()

# Adult barplots
pheno_colors = sns.color_palette("icefire", 15)
for par in ['Carrier parent', 'Noncarrier parent']:
    last=[0]*len(adult_domains)
    for score in range(0, 2):
        vals = []
        sizes=[]
        for pheno in adult_domains:
            num = df[(df[pheno]==score) & (df.Group==par)].shape[0]*100/df[(~df[pheno].isnull()) & (df.Group==par)].shape[0]
            vals.append(num)
            size=df[(~df[pheno].isnull()) & (df.Group==par)].shape[0]
            sizes.append(size)
        if par=='Carrier parent':
            color=pheno_colors[-1*(score+1)]
            plt.bar([i-0.2 for i in range(len(adult_domains))], vals, bottom=last, color=color, width=0.4)
        else:
            color=pheno_colors[score]
            plt.bar([i+0.2 for i in range(len(adult_domains))], vals, bottom=last, color=color, width=0.4)
        if score==0:
            for i in range(len(sizes)):
                if par=='Carrier parent':
                    plt.text(i-0.25, 105, 'n='+str(sizes[i]), ha='center')
                else:
                    plt.text(i+0.25, 105, 'n='+str(sizes[i]), ha='center')
        newlast= [vals[i]+last[i] for i in range(len(last))]
        last=newlast
plt.xlabel("Phenotypic Domain")
plt.ylabel("Percentage of Parents")
plt.legend(['Carrier_'+str(i) if i<5 else 'NonCarrier_'+str(i-5) for i in range(0, 10)], bbox_to_anchor=(1, 1))
plt.title('Parent Phenotypes')
plt.ylim(0, 110)
plt.xticks(ticks = [i for i in range(len(adult_domains))], labels=[i.split(' (')[0] for i in adult_domains], rotation=90)
pdf.savefig()
plt.close()

pdf.close()
