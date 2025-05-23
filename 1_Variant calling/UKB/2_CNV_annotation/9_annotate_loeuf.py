import pandas as pd

# Annotate CNV genes with LOEUF score

# Input and ouput files
input_file="/path/to/input/cnvs.txt" # Use the output of script 8_explode_genes.py here
output_file="/path/to/output/file.txt" # Note that the output of this script will be the final CNV annotations

# Reference files
loeuf_file='/path/to/parsed/gnomad/loeuf/annotations.csv'

# Load files
calls = pd.read_csv(input_file, sep='\t')
gnomad_df = pd.read_csv(loeuf_file, sep='\t')

# Annotate variants with LOEUF score
# If a variant is annotated with multiple genes/LOEUF scores, return the lowest
def get_loeuf(row):
	# Get IDs
	ids = row['Gene_ID'].split(' ')

	# Restrict to only IDs in gnomAD
	ids2 = [i for i in ids if i in gnomad_df.gene_id.to_list()]

	# No IDs are in gnomAD list
	if len(ids2)==0:
		return '.'

	loeuf = min(gnomad_df[gnomad_df.gene_id.isin(ids2)]['oe_lof_upper'].to_numpy())
	return loeuf

calls['LOEUF'] = calls.apply(get_loeuf, axis = 1)

# Save
calls.to_csv(output_file, sep = '\t', index = False)
