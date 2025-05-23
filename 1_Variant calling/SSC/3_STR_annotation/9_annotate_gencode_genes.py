import pandas as pd

# Input and output files
INPUT_TABLE="/path/to/input/table.tsv" # Use the output of script 8_filter_calls.py
OUTPUT_TABLE="/path/to/output/table.csv"

# Reference files
GENCODE='/path/to/parsed/gencode/annotations.csv'

# Load files
df=pd.read_csv(INPUT_TABLE, sep='\t')

gencode_df = pd.read_csv(GENCODE)

gencode_chrom_array = gencode_df['Chrom'].to_numpy()
gencode_start_array = gencode_df['Start'].to_numpy()
gencode_end_array   = gencode_df['End'].to_numpy()

# Rename columns
columns = ['Chrom', 'Pos', 'Sample', 'zscore', 'longest_allele', 'both_alleles',
       'cohort_mean', 'cohort_mode', 'cohort_std', 'end', 'motif', 'period',
       'variant_id', 'Func.wgEncodeGencodeBasicV38',
       'Gene.wgEncodeGencodeBasicV38', 'GeneDetail.wgEncodeGencodeBasicV38',
       'ExonicFunc.wgEncodeGencodeBasicV38',
       'AAChange.wgEncodeGencodeBasicV38']
df.columns = columns

# Helper functions
def get_gene_attribute(chrom, pos, attribute='gene_name', offset=0):
	same_chrom = gencode_chrom_array == chrom
	after_start = gencode_start_array <= pos + offset
	before_end = gencode_end_array >= pos - offset
	
	overlap = same_chrom & after_start & before_end
	
	gene_attributes = gencode_df[overlap][attribute].to_list()
	
	# increase search space if none found
	if len(gene_attributes) == 0:
		print(chrom+'_'+str(pos))
		return get_gene_attribute(chrom, pos, attribute, offset=offset+10)
		
	return ';'.join(gene_attributes)

# Annotate calls with genes and attributes
df['Gene_symbol'] = df[['Chrom', 'Pos']].apply(lambda s: get_gene_attribute(s[0], s[1], 'gene_name'), axis=1)
df['Gene_id'] = df[['Chrom', 'Pos']].apply(lambda s: get_gene_attribute(s[0], s[1], 'gene_id'), axis=1)
df['Gene_biotype'] = df[['Chrom', 'Pos']].apply(lambda s: get_gene_attribute(s[0], s[1], 'gene_type'), axis=1)

# Remove any genes not annotated by ANNOVAR
for i, row in df.iterrows():
	annovar_genes = row['Gene.wgEncodeGencodeBasicV38'].split(';')
	gencode_genes = row['Gene_symbol'].split(';')
	gencode_ids = row['Gene_id'].split(';')
	gencode_biotypes = row['Gene_biotype'].split(';')
	
	# create dict of gene ids and biotypes
	symbol2id = {}
	symbol2biotype = {}
	for j in range(len(gencode_genes)):
		gencode_gene = gencode_genes[j]
		gene_id = gencode_ids[j]
		gene_biotype = gencode_biotypes[j]
		
		symbol2id[gencode_gene] = gene_id
		symbol2biotype[gencode_gene] = gene_biotype
	
	# remove gencode genes that aren't annotated by annovar
	genes_to_remove = []
	for gene in gencode_genes:
		if gene not in annovar_genes:
			genes_to_remove.append(gene)
	
	for gene in genes_to_remove:
		gencode_genes.remove(gene)
	
	# sort gencode and annovar genes
	gencode_genes = sorted(list(set(gencode_genes)))
	annovar_genes = sorted(list(set(annovar_genes)))
	
	# make sure gene ids and biotypes match order of genes
	gencode_ids = [symbol2id[s] for s in gencode_genes]
	biotype_ids = [symbol2biotype[s] for s in gencode_genes]
	
	df.at[i, 'Gene_symbol'] = ';'.join(gencode_genes)
	df.at[i, 'Gene.wgEncodeGencodeBasicV38'] = ';'.join(annovar_genes)
	df.at[i, 'Gene_id'] = ';'.join(gencode_ids)
	df.at[i, 'Gene_biotype'] = ';'.join(biotype_ids)

# Restrict to protein coding genes
df = df[df['Gene_biotype'].apply(lambda s: 'protein_coding' in s)]

# Save
df.to_csv(OUTPUT_TABLE, index=False)
