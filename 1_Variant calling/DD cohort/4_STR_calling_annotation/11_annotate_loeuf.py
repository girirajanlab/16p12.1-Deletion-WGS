import pandas as pd

# Input and output files
INPUT='/path/to/input/file.csv' # Use the output of script 10_annotate_gencode_genes.py here
OUTPUT='/path/to/output/table.csv' # Note that the output of this script will be the final coding SNV annotations
loeuf_file='/path/to/parsed/gnomad/loeuf/annotations.csv'

# Load files
df = pd.read_csv(INPUT)
gnomad_df = pd.read_csv(loeuf_file, sep='\t')
gnomad_df = gnomad_df.set_index('gene_id')

# Remove the version number from gene IDs
df['Gene_id_'] = df['Gene_id'].apply(lambda s: s.split('.')[0])

# Annotate LOEUF scores
df['LOEUF'] = ''
for i, row in df.iterrows():
	gene_id = row['Gene_id_']
	if gene_id not in gnomad_df.index:
		continue
	df.at[i, 'LOEUF'] = gnomad_df.at[gene_id, 'oe_lof_upper']

# Save
df.to_csv(OUTPUT, index=False)
