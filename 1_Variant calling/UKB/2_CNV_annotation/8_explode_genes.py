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

calls['call_ID']=calls.variant_id
calls['Pathogenic']=calls.Pathogenic_Name

for i in range(max_genes):
	subdf=calls[calls.gene_count>=i+1].copy()
	subdf['Gene_ID']=subdf.gene_ids.str.split(' ', expand=True)[i]
	subdf['Gene_Symbol']=subdf.gene_symbols.str.split(' ', expand=True)[i]

	gene_df=pd.concat([gene_df, subdf[['Sample', 'Gene_ID', 'Gene_Symbol', 'Type', 'Pathogenic', 'call_ID']]], axis=0, ignore_index=True)

gene_df.to_csv(output_file, sep = '\t', index = False)
