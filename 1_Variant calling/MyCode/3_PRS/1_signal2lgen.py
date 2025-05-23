import pandas as pd

# Convert the BAF values in the singal file to MAP and LGEN files for plink

# Input and output files
FILE_LIST="/path/to/list/of/signal/files.csv" # File list includes information on sample ID, sex and array file name
OMNI_SNP="/path/to/OMNI/SNP/Position/file.txt"
GSA_SNP="/path/to/GSA/SNP/Position/file.txt"
OUTPUT_DIR="/path/to/output/directory"

# Get singal files
files=pd.read_csv(file_list)

# Load SNP positiongs
# OMNI SNP locations from BIM file provided by Matt Oetjens
omnisnp=pd.read_csv(OMNI_SNP, sep='\t',
					header=None, names=['Chr', 'Name', 'Pos_CM', 'Pos', 'A1', 'A2'])
# GSA SNP locations from BIM file provided by Matt Oetjens
gsasnp=pd.read_csv(GSA_SNP, sep='\t',
					header=None, names=['Chr', 'Name', 'Pos_CM', 'Pos', 'A1', 'A2'])

# Annotate with A+B alleles for each site
def anno_b(df, snpdf):
	out=pd.merge(df, snpdf[['Name', 'A1', 'A2']], on='Name', how='inner')
	return out

# Use BAF information to determine genotypes:
# 0-0.25: AA
# 0.25-0.75: AB
# 0.75-1: BB

for idx, row in files.iterrows():
	file=row['signal_filename']
	type=''
	if 'OMNI' in row['array']:
		type='OMNI'
	if 'GSA' in row['array']:
		type='GSA'
	sample=row['SEQN_ID']

	print(file, type, sample)

	df=pd.read_csv(file, sep='\t')
	# Rename columns
	cols=df.columns.to_list()
	cols=[]
	for c in df.columns.to_list():
		if 'Log R' in c:
			cols.append('LRR')
		elif 'B Allele' in c:
			cols.append('BAF')
		else:
			cols.append(c)
	df.columns=cols

	if type=='OMNI':
		df=anno_b(df, omnisnp)
	elif type=='GSA':
		df=anno_b(df, gsasnp)

	df['Allele_call_1']=df['A1']
	df.loc[df.BAF>0.75, 'Allele_call_1']=df.loc[df.BAF>=0.75, 'A2']
	df['Allele_call_2']=df['A1']
	df.loc[df.BAF>0.25, 'Allele_call_2']=df.loc[df.BAF>0.25, 'A2']

	# Save in LGEN format
	df['FID']=sample
	df['IID']=sample

	df[['FID', 'IID', 'Name', 'Allele_call_1', 'Allele_call_2']].to_csv(f'{OUTPUT_DIR}/{type}/{sample}.lgen', sep=' ', index=False, header=False)

# Also make FAM and MAP files

# FAM
files['FID']=files.SEQN_ID
files['sex']=files.PT_SEX.map({'Female':2, 'Male':1})
files['Mother']=0
files['Father']=0
files['Phenotype']=0
files[['FID', 'SEQN_ID', 'Father', 'Mother', 'sex', 'Phenotype']].to_csv(f'{OUTPUT_DIR}/signal_to_lgen.fam', sep=' ', index=False, header=False)

# MAP
# Make separate MAP for OMNI and GSA
omnisnp[['Chr', 'Name', 'Pos_CM', 'Pos']].to_csv(f'{OUTPUT_DIR}/singal_to_lgen_OMNI.map', sep=' ', index=False, header=False)
gsasnp[['Chr', 'Name', 'Pos_CM', 'Pos']].to_csv(f'{OUTPUT_DIR}/singal_to_lgen_GSA.map', sep=' ', index=False, header=False)
