import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Wedge
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.colors as mcolors

matplotlib.rcParams['pdf.fonttype'] = 42

# Check for cases where a proband and their non-carrier sibling share a "dual diagnosis" variant

# Input and output files
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
PATHOGENIC_SNVS="/output/final/pathogenic/SNVs.csv" # Use the output of script 2_update_pathogenic_SNVs.py
CNVS="/path/to/CNV/calls.txt" # Use the output from script 1_Variant calling/DD cohort/3_CNV_calling_annotation/merge_all_cnvs/2_annotate_loeuf.py
SNVS="/path/to/SNV/variants.csv" # Use the output of script 1_Variant calling/DD cohort/2_SNV_annotation/coding_annotations/14_loeuf_scores.py

OUTPUT_PDF="/path/to/ouput/pedigrees.pdf" # These will be the pedigrees shown in Fig. S2C

# Load files
df=pd.read_csv(TABS1A)
fam_map=dict(zip(df.Sample.to_list(), df.Family.to_list()))

# Annotate families for pathogenic variants
psnv=pd.read_csv(PATHOGENIC_SNVS)
psnv['Family']=psnv.Sample.map(fam_map)
psnv['variant_id']=psnv.Chrom+'_'+psnv.Pos.astype(str)+'_'+psnv.Ref+'_'+psnv.Alt+'.'+psnv.Family

pcnv=pd.read_csv(CNVS, sep='\t')
pcnv=pcnv[pcnv.NEJM!='.'][['Sample', 'NEJM']]

pcnv['Family']=pcnv.Sample.map(fam_map)
pcnv['variant_id']=pcnv.NEJM+'.'+pcnv.Family
pcnv.drop_duplicates(inplace=True)

# Identify pathogenic SNVs in family members
snv=pd.read_csv(SNVS)
snv['Family']=snv.Sample.map(fam_map)
snv['variant_id']=snv.Chrom+'_'+snv.Pos.astype(str)+'_'+snv.Ref+'_'+snv.Alt+'.'+snv.Family

snv=snv[snv.variant_id.isin(psnv.variant_id.to_list())]

# Check for non-carrier siblings who have "dual diagnosis" hits
ncsibs=df[(df.Relationship=='Sibling') & (df['16p12.1 deletion']=='Noncarrier')].Sample.to_list()
sibsamps=snv[snv.Sample.isin(ncsibs)].Sample.to_list()+pcnv[pcnv.Sample.isin(ncsibs)].Sample.to_list()

# Check the relevant families
fams=df[df.Sample.isin(sibsamps)].Family.to_list()
df=df[df.Family.isin(fams)]

# Annotate each family member for whether they have pathogenic variants
df['Patho_var']=((df.Sample.isin(snv.Sample.to_list())) | (df.Sample.isin(pcnv.Sample.to_list())))

# Ensure the proband has a pathogenic variant
ppros=df[(df.Patho_var) & (df.Relationship=='Proband')]
df=df[df.Family.isin(ppros.Family.to_list())]

# For each family, create pedigrees showing the pathogenic variants, deletion, and phenotypes
child_domains=['ID/DD (Child domain)', 'Behavioral features (Child domain)', 'Psychiatric features (Child domain)', 'Nervous System Abnormalities (Child domain)', 'Congenital Anomalies (Child domain)', 'Growth/Skeletal Defects (Child domain)']
adult_domains=['Anxiety (Questionnaire)', 'Sleep trouble (Questionnaire)', 'Psychosis (Questionnaire)', 'Addiction (Questionnaire)', 'Mood lability (Questionnaire)']
df=df[['Sample', 'Family', 'Relationship', 'Sex', '16p12.1 deletion']+child_domains+adult_domains]

# Annotate the specific pathogenic variants in the families
vardf=pd.concat([pcnv[['Sample', 'NEJM', 'variant_id']], snv[['Sample', 'Gene_id', 'Gene_symbol', 'variant_id']]])
vardf=vardf[vardf.Sample.isin(df.Sample.to_list())]
vardf.fillna('', inplace=True)
vardf['vid']=vardf.variant_id+'.'+vardf.NEJM+vardf.Gene_symbol

vids=list(vardf.vid.unique())
for vid in vids:
	df[vid]=df.Sample.isin(vardf[vardf.vid==vid].Sample.to_list())

# Create pedigrees
def plot_shape(ax, sex, x, y, shape_width=0.1, fc='None', ec='k', lw=1):
	# location represents the center of the shape
	if sex=='M':
		ax.add_artist(Rectangle((x-(shape_width/2), y-(shape_width/2)), shape_width, shape_width, fc=fc, edgecolor=ec, lw=lw, zorder=2))
	elif sex=='F':
		ax.add_artist(Circle((x, y), shape_width/2, fc=fc, edgecolor=ec, lw=lw, zorder=2))

def x_positions(df, m_pos=0.6, f_pos=0.4, sib_width=0.1, au_width=0.17, shape_width=0.1, partner_width=0.1):
	pos_dict={}
	pos_dict[df[df.Relationship=='Mother'].Sample.to_list()[0]]=m_pos
	pos_dict[df[df.Relationship=='Father'].Sample.to_list()[0]]=f_pos
	
	# Get child positions
	nkids=df[df.Relationship.isin(['Sibling', 'Proband'])].shape[0]
	total_width=(shape_width*nkids)+(sib_width*(nkids-1))
	par_mid=(m_pos+f_pos)/2
	kids=df[df.Relationship=='Proband'].Sample.to_list()+df[df.Relationship=='Sibling'].Sample.to_list()
	start_loc=par_mid-(total_width/2)+(shape_width/2)
	for k in kids:
		pos_dict[k]=start_loc
		start_loc+=sib_width+shape_width
	
	# Get extended family positions
	for fs in ['M', 'F']:
		par_loc={'M':m_pos, 'F':f_pos}[fs]
		# Aunt/Uncle positions
		au=df[(df.Family_side==fs) & (df.Relationship.isin(['Aunt', 'Uncle']))].Sample.to_list()
		direction={'M':1, 'F':-1}[fs]
		start_loc=par_loc+(direction*au_width)
		for a in au:
			pos_dict[a]=start_loc
			start_loc+=(direction*au_width)
	
		# Grandparent positions
		for sex in ['M', 'F']:
			gp=df[(df.Family_side==fs) & (df.Sex==sex) & (df.Relationship.isin(['Grandmother', 'Grandfather']))].Sample.to_list()
			if len(gp)==0:
				continue
			gp=gp[0]
			xmid=par_loc+(direction*((len(au)*au_width)/2))
			gdirection={'M':-1, 'F':1}[sex]
			xloc=xmid+(gdirection*partner_width)
			pos_dict[gp]=xloc
	
	return df.Sample.map(pos_dict)

def draw_connecting(df, ax, buffer=0.1):
	# Draw lines connecting siblings
	kiddf=df[df.Relationship.isin(['Proband', 'Sibling'])]
	kid_min=kiddf.x.min()
	kid_max=kiddf.x.max()
	kid_y=kiddf.y.min()
	ax.plot([kid_min, kid_max], [kid_y+buffer, kid_y+buffer], color='k')
	
	# Connect kids to parents
	kid_mid=(kid_min+kid_max)/2
	ax.plot([kid_mid, kid_mid], [kid_y+buffer, kid_y+(buffer*2)], color='k')
	
	# Draw lines connecting aunts/uncles/parents
	for fs in ['M', 'F']:
		par={'M':'Mother', 'F':'Father'}[fs]
		fsdf=df[(df.Relationship==par) | ((df.Relationship.isin(['Aunt', 'Uncle'])) & (df.Family_side==fs))]
		pau_min=fsdf.x.min()
		pau_max=fsdf.x.max()
		pau_y=fsdf.y.min()
		ax.plot([pau_min, pau_max], [pau_y+buffer, pau_y+buffer], color='k')
	
		# Connnect aunts/uncles/parents to grandparents, if available
		gp=df[(df.Family_side==fs) & (df.Relationship.isin(['Grandmother', 'Grandfather']))]
		if gp.shape[0]==0:
			continue
		pau_mid=(pau_min+pau_max)/2
		ax.plot([pau_mid, pau_mid], [pau_y+buffer, pau_y+(buffer*2)], color='k')

def plot_phenos(ax, sex, x, y, phenos, shape_width=0.1):
	# Add in blocks of color to represent phenotypes
	pheno2=[i.replace(' (Questionnaire)', '').replace(' (Child domain)', '') for i in phenos]
	domain_colors={'ID/DD':'#6A8A82', 'Behavioral features':'#9385BF', 'Psychiatric features':'#F9B163', 'Nervous System Abnormalities':'#81B1D2', 'Congenital Anomalies':'#FFCE00', 'Growth/Skeletal Defects':'#EF7F71',
					'Anxiety':'#50723C','Sleep trouble': '#3D7EA6','Psychosis':'#AB6C82', 'Addiction':'#C34A36', 'Mood lability':'#D88C1C'}
	if sex=='M':
		# For males, draw rectangles with width=rect_width/(number of phenotypes)
		w=shape_width/len(pheno2)
		locx=x-(shape_width/2)
		for p in pheno2:
			ax.add_artist(Rectangle((locx, y-(shape_width/2)), w, 0.1, facecolor=domain_colors[p], edgecolor=None, zorder=1))
			locx+=w
	elif sex=='F':
		# For females, draw pie slices for phenotypes
		d_theta=360/len(pheno2)
		theta=90
		for p in pheno2:
			ax.add_artist(Wedge((x, y), shape_width/2, theta, theta+d_theta, facecolor=domain_colors[p], edgecolor=None, zorder=1))
			theta+=d_theta

df['Family_side']='M'

def draw_pedigree(df, family, pdf):
	subdf=df[df.Family==family].copy()

	# Define positions
	subdf['y']=subdf.Relationship.map({'Proband':0.15, 'Sibling':0.15, 'Cousin':0.15, 'Mother':0.35, 'Father':0.35, 'Aunt':0.35, 'Uncle':0.35, 'Grandmother':0.55, 'Grandfather':0.55})
	subdf['x']=x_positions(subdf)

	# Identify VIDs for the family
	fam_vids=[i for i in vids if subdf[i].any(axis=0)]
	clean_vids=[]
	fv_colors={}
	for i in range(len(fam_vids)):
		fv=fam_vids[i]
		# Clean up for mapping
		chrom=fv.split('_')[0]
		pos=fv.split('_')[1]
		ref=fv.split('_')[2]
		alt=fv.split('_')[3].split('.')[0]
		ng=fv.split('.')[2]
		
		fv_clean=f'{chrom}:{pos} {ref}>{alt} {ng}'
		subdf[fv_clean]=subdf[fv]
		clean_vids.append(fv_clean)
		
		fv_colors[fv_clean]=list(mcolors.TABLEAU_COLORS)[i]

	# For each member, plot shapes
	fig, ax = plt.subplots(figsize=(6,6), layout='constrained')
	shape_width=0.1
	partner_width=0.1
	for idx, row in subdf.iterrows():
		sex=row.Sex
		ec={'Carrier':'red', 'Noncarrier':'k'}[row['16p12.1 deletion']]
		lw={'Carrier':2, 'Noncarrier':1}[row['16p12.1 deletion']]
		plot_shape(ax, sex, row.x, row.y, ec=ec, lw=lw)
		
		# Add Sample name
		ax.text(row.x, row.y-0.075, row.Sample, ha='center')
		
		# Draw lines
		# If child, draw a line up
		if row.Relationship in ['Sibling', 'Proband', 'Cousin']:
			ax.plot([row.x, row.x], [row.y+shape_width/2, row.y+0.1], color='k')
		for i in range(2):
			fs=['M', 'F'][i]
			par=['Mother', 'Father'][i]
			if subdf[(subdf.Relationship.isin(['Grandmother', 'Grandfather', 'Aunt', 'Uncle'])) & (subdf.Family_side==fs)].shape[0]>0:
				if row.Relationship==par or (row.Relationship in ['Aunt', 'Uncle'] and row.Family_side==fs):
					ax.plot([row.x, row.x], [row.y+shape_width/2, row.y+0.1], color='k')
		# If parent, draw line to midpoint
		if row.Relationship in ['Mother', 'Father', 'Grandmother', 'Grandfather']:
			direction={'M':1, 'F':-1}[sex]
			ax.plot([row.x+(direction*shape_width/2), row.x+(direction*shape_width/2)+(direction*partner_width/2)], [row.y, row.y], color='k')
		
		# Annotate domain and questionnaire phenotypes
		# If phenotypes are all null, fill in grey
		if row.Relationship in ['Sibling', 'Proband', 'Cousin']:
			phenos=row[child_domains].index[row[child_domains]>0].to_list()
			if row[child_domains].isnull().all():
				plot_shape(ax, sex, row.x, row.y, fc='lightgrey', ec='None')
		else:
			phenos=row[adult_domains].index[row[adult_domains]>0].to_list()
			if row[adult_domains].isnull().all():
				plot_shape(ax, sex, row.x, row.y, fc='lightgrey', ec='None')
		
		if len(phenos)>0:
			plot_phenos(ax, sex, row.x, row.y, phenos)
		
		# Annotate pathogenic variants
		if row[clean_vids].sum()>0:
			nvids=row[clean_vids].sum()
			step=shape_width/(nvids+1)
			mult=1
			for fv in clean_vids:
				if row[fv]:
					ax.text(row.x-(shape_width/2)+(mult*step), row.y-0.08, '*', ha='center', va='top', color=fv_colors[fv], weight='bold')
					mult+=1
			

	# Draw connecting lines
	draw_connecting(subdf, ax)

	# Add in a custom legend
	domain_colors={'ID/DD':'#6A8A82', 'Behavioral features':'#9385BF', 'Psychiatric features':'#F9B163', 'Nervous System Abnormalities':'#81B1D2', 'Congenital Anomalies':'#FFCE00', 'Growth/Skeletal Defects':'#EF7F71',
					'Anxiety':'#50723C','Sleep trouble': '#3D7EA6','Psychosis':'#AB6C82', 'Addiction':'#C34A36', 'Mood lability':'#D88C1C'}
	y=0.9
	x=0.05
	for d in child_domains:
		plot_shape(ax, 'M', x, y, shape_width=0.02, fc=domain_colors[d.split(' (')[0]], ec='None')
		ax.text(x+0.02, y, d.split(' (')[0], ha='left', va='center')
		y+=-0.04
	y=0.9
	x=0.5
	for d in adult_domains:
		plot_shape(ax, 'M', x, y, shape_width=0.02, fc=domain_colors[d.split(' (')[0]], ec='None')
		ax.text(x+0.02, y, d.split(' (')[0], ha='left', va='center')
		y+=-0.04

	plot_shape(ax, 'M', x, y, shape_width=0.02, fc='lightgrey', ec='None')
	ax.text(x+0.02, y, 'No data', ha='left', va='center')

	y=0.65
	x=0.05
	for fv in clean_vids:
		ax.text(x, y, f'* {fv}', ha='left', va='center', color=fv_colors[fv])
		y+=-0.04

	y=0.95
	x=0.05
	plot_shape(ax, 'M', x, y, shape_width=0.02, fc='None', ec='red', lw=2)
	ax.text(x+0.02, y, '16p12.1 deletion carrier', ha='left', va='center')

	plt.xlim(0, 1)
	plt.ylim(0, 1)

	ax.axis('off')

	pdf.savefig()
	plt.close()

pdf=PdfPages(OUTPUT_PDF)
for fam in list(df.Family.unique()):
	draw_pedigree(df, fam, pdf)

pdf.close()
