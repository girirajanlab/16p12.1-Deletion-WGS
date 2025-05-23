import pandas as pd

# Explode CNV calls by genes

# Input and output files
INPUT='path/to/input/cnvs.csv' # Use the output of script 6_annotate_gencode.py
OUTPUT='path/to/output/file.csv'

# Load file
df=pd.read_csv(INPUT)

# Remove any CNVs that do not affect genes
df = df[~df.Genes.isna()]

# Explode by gene
for i, row in df.iterrows():
	sample = row['PatientID']
	genes = row['Genes']
	gene_ids = row['Gene_ids']
	dup_or_del = row['Type']
	
	genes = genes.split(';')
	gene_ids = gene_ids.split(';')
	for j in range(len(genes)):
		gene = genes[j]
		gene_id = gene_ids[j]
		app = [sample, gene, gene_id, dup_or_del]
		new_list.append(app)

df = pd.DataFrame(new_list, columns=['Sample', 'Gene', 'Gene_id', 'CNV_Type'])

# Standardize sample and CNV names
def std_name(s):
	if '.' in s:
		return s
	items = s.split('x')
	s = items[0] + '.x' + items[1]
	return s


def std_cnv_type(s):
	if s == 'dup':
		return 'Dup'
	if s == 'del':
		return 'Del'
	print(s)

df['Sample'] = df.Sample.apply(std_name)
df['CNV_Type'] = df.CNV_Type.apply(std_cnv_type)

# Drop any duplicate variants
df = df.drop_duplicates(['Sample', 'Gene_id', 'CNV_Type'])

# Save
df.to_csv(OUTPUT)

