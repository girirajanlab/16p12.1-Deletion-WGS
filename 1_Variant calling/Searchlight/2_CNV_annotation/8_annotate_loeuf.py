import pandas as pd

# Input and output files
input_file='/path/to/input/file.csv' # Use the output of script 7_explode_genes.py here
output_file='/path/to/output/table.csv' # Note that the output of this script will be the final CNV annotations
loeuf_file='/path/to/parsed/gnomad/loeuf/annotations.csv'

# Load files
df = pd.read_csv(input_file)
gnomad_df = pd.read_csv(loeuf_file, sep='\t')
gnomad_df = gnomad_df.set_index('gene_id')

# Annotate variants with LOEUF score
# If a variant is annotated with multiple genes/LOEUF scores, return the lowest
df['LOEUF'] = ''
for i, row in df.iterrows():
	gene_ids = row['Gene_id_'].split(';')
	gene_ids_in_gnomad = [s for s in gene_ids if s in gnomad_df.index]
	
	if len(gene_ids_in_gnomad) == 0:
		# none of the gene ids are in gnomad
		print(gene_ids)
		continue
		
	lowest_loeuf = 1000
	for gene_id in gene_ids_in_gnomad:
		loeuf = gnomad_df.at[gene_id, 'oe_lof_upper']
		if loeuf < lowest_loeuf:
			lowest_loeuf = loeuf
			
			
	df.at[i, 'LOEUF'] = lowest_loeuf

# Save
df.to_csv(output_file, index=False)
