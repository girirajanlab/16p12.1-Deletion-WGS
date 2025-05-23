import pandas as pd

# Annotate gencode genes

# Input and output files
INPUT='path/to/input/cnvs.txt' # Use the output from script 2_anno_segdup_centel.sh
OUTPUT='path/to/output/file.csv' # This is the final output for SSC CNV filtering

# Reference files
gencode_path='/path/to/parsed/gencode/exon/file/gencode.v19.parsed.exons.csv'
# This file contains all exons from the Gencode v19 GTF, but parsed into a table format

# Load files
df=pd.read_csv(INPUT, sep='\t')
gencode = pd.read_csv(gencode_path)

# Filter genes for only protein coding
gencode = gencode[gencode.gene_type=='protein_coding']

gene_chroms = gene_df.Chrom.to_numpy()
gene_starts = gene_df.Start.to_numpy()
gene_ends = gene_df.End.to_numpy()

# Get genes
df['Genes'] = ''
for i, row in df.iterrows():
	chrom = row['Chr']
	start = row['Start(hg19)']
	end = row['Stop(hg19)']
	
	min_end = np.minimum(gene_ends, end)
	max_start = np.maximum(gene_starts, start)
	intersect = (gene_chroms == chrom) & (min_end - max_start >= 0)
	intersect_df = gene_df.loc[intersect]
	genes = intersect_df['gene_name'].to_list()
	genes = list(set(genes))
	df.at[i, 'Genes'] = ';'.join(genes)
	
	genes = intersect_df['gene_id'].to_list()
	genes = list(set(genes))
	df.at[i, 'Gene_ids'] = ';'.join(genes)

# Save
df.to_csv(OUTPUT, index=False)
