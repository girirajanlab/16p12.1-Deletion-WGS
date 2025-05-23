import pandas as pd

# Annotate gencode genes

# Input and output files
INPUT='path/to/input/cnvs.txt' # Use the output from script 9_filter_calls.py
OUTPUT='path/to/output/file.txt'

# Reference files
gencode_path='/path/to/parsed/gencode/exon/file/gencode.v19.parsed.exons.csv'
# This file contains all exons from the Gencode v19 GTF, but parsed into a table format

# Load files
df=pd.read_csv(INPUT, sep='\t')
gencode = pd.read_csv(gencode_path)

# Filter genes
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


df=df.apply(get_genes, axis=1)

# Remove any CNVs that do not overlap genes
df=df[df.gene_ids!='.']

df.to_csv(OUTPUT, sep='\t', index=False)
