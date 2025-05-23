import pandas as pd

# Annotate gencode genes

# Input and output files
input_path='path/to/cnv/files' # These are the files from script 10_merge_files.sh
output_path='path/to/output/files'

# Reference files
gencode_path='/path/to/parsed/gencode/exon/file/gencode.v19.parsed.exons.csv'
# This file contains all exons from the Gencode v19 GTF, but parsed into a table format

gencode = pd.read_csv(gencode_path)
gencode = gencode[gencode.gene_type=='protein_coding']

# Get genes
def get_genes(s):
	chr = s.Chr
	start = s.Start
	end = s.End

	# Get exons from gencode file that are at least partially contained within the CNV length
	gene_df = gencode[(gencode.Chrom==chr) & (((gencode.Start>=start) & (gencode.Start<=end)) | ((gencode.End>=start) & (gencode.End<=end)))]

	if gene_df.shape[0]>0:
		# Get a list of gene IDs and gene names
		gene_info = gene_df[['gene_id', 'gene_name']].drop_duplicates()
		s['gene_ids'] = ' '.join(gene_info.gene_id)
		s['gene_names'] = ' '.join(gene_info.gene_name)

	else:
		s['gene_ids'] = '.'
		s['gene_names'] = '.'

	return s

for cnv_type in ['cnv', 'pathogenic_cnv']:
	df=pd.read_csv(f'{input_path}/filtered_{cnv_type}.bed', sep='\t')
	df=df.apply(get_genes, axis=1)
	df.to_csv(f'{output_path}/{cnv_type}_annotated_genes.bed', sep='\t', index=False)
