import pandas as pd

# Separate calls into gene-level calls

# Input and output files
input_file='path/to/cnv/file.txt' # This is the output from script 11_inhertiance_annotation.py
output_file='path/to/output/file.txt'

# Reference files
gencode_path='/path/to/parsed/gencode/exon/file/gencode.v19.parsed.genes.csv'
# This file contains all genes from the Gencode v19 GTF, but parsed into a table format

# Load calls
calls = pd.read_csv(input_file, sep = '\t')

# Gene annotations
genes = pd.read_csv(gencode_path)

gene_df = pd.DataFrame(columns = ['Sample', 'Gene_ID', 'Gene_Symbol', 'Type', 'Pathogenic_Name', 'Inheritance'])
for idx, row in calls.iterrows():
	# Get call information
	sample = row['PatientID']
	gene_ids = row['gene_ids'].split(' ')
	gene_symbols = row['gene_names'].split(' ')
	type = row['Type']
	pathogenic = row['Pathogenic_Name']
	inh=row['inheritance']

	out_lines = []

	# Make sure gene ID list and gene symbol list are the same length
	if len(gene_ids) != len(gene_symbols):
		print('Gene IDs and symbols do not match!')
		print(' '.join([str(i) for i in row]))
		print(' ')

		# For cases where symbols and ids don't match, manually link them
		# All cases have more gene IDs than gene symbols, so look up the gene IDs and get the corresponding symbols
		gene_symbols = []
		for id in gene_ids:
			gene_symbols.append(genes[genes.gene_id==id]['gene_name'].to_string(index = False, header = False).strip())

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
			out_lines.append([sample, ' '.join(id_dict[key]), key, type, pathogenic, inh])

	else:
		# Iterate through gene lists

		# Remove version number from gene IDs
		gene_ids2 = [i.split('.')[0] for i in gene_ids]

		for i in range(len(gene_ids2)):
			out_lines.append([sample, gene_ids2[i], gene_symbols[i], type, pathogenic, inh])

	# Add to dataframe
	gene_df = pd.concat([gene_df, pd.DataFrame(out_lines, columns = ['Sample', 'Gene_ID', 'Gene_Symbol', 'Type', 'Pathogenic_Name', 'Inhertiance'])], axis = 0, ignore_index = True)

gene_df.to_csv(output_file, sep = '\t', index = False)
