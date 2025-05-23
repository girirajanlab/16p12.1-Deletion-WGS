import pandas as pd

# Separate calls into gene-level calls

# Input and ouput files
input_file="/path/to/input/cnvs.bed" # Use the output of script 7_annotate_gencode.py here
output_file="/path/to/output/file.txt"

# Load files
calls=pd.read_csv(input_file, sep='\t')

# Explode by gene
calls['gene_count']=calls.gene_ids.str.count(' ')+1
calls['symbol_count']=calls.gene_symbols.str.count(' ')+1
max_genes=max(calls.gene_count.unique())

gene_df = pd.DataFrame(columns = ['Sample', 'Gene_ID', 'Gene_Symbol', 'Type', 'Pathogenic', 'call_ID'])
for idx, row in calls.iterrows():
	# Get call information
	sample = row['Sample']
	gene_ids = row['gene_ids'].split(' ')
	gene_symbols = row['gene_symbols'].split(' ')
	type = row['Type']
	ID = row['variant_id']
	pathogenic = row['Pathogenic_Name']

	# Exclude 16p12.1 CNVs
	if pathogenic=='16p12.1':
		continue

	out_lines = []

	# NOTE: There are some cases where a gene symbol has multiple IDs!!
	# Remove version number from gene IDs
	gene_ids2 = [i.split('.')[0] for i in gene_ids]

	# We only want to report the unique gene symbols
	id_dict = {}
	for n, symbol in enumerate(gene_symbols):
		if symbol not in id_dict:
			id_dict[symbol] = [gene_ids2[n]]
		else:
			id_dict[symbol].append(gene_ids2[n])

	for key in id_dict.keys():
		out_lines.append([sample, ' '.join(id_dict[key]), key, type, pathogenic, ID])

	# Add to dataframe
	gene_df = pd.concat([gene_df, pd.DataFrame(out_lines, columns = ['Sample', 'Gene_ID', 'Gene_Symbol', 'Type', 'Pathogenic', 'call_ID'])], axis = 0, ignore_index = True)

gene_df.to_csv(output_file, sep = '\t', index = False)