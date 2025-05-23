import pandas as pd

# Annotate gencode genes

# Input and output files
INPUT_OMNI1='path/to/input/omni1/cnvs.xlsx'
INPUT_OMNI2='path/to/input/omni2/cnvs.xlsx'
OUTPUT='path/to/output/file.csv'

# Reference files
gencode_path='/path/to/parsed/gencode/exon/file/gencode.v19.parsed.exons.csv'
# This file contains all exons from the Gencode v19 GTF, but parsed into a table format

# Load files
omni1=pd.read_excel(INPUT_OMNI)
omni2=pd.read_excel(INPUT_OMN2)
gencode = pd.read_csv(gencode_path)

# Concat batches
df=omni1.append(omni2)
df=df.reset_index(drop=True)

# Filter genes
gencode = gencode[gencode.gene_type=='protein_coding']

# Annotate genes
gene_chroms = gene_df.Chrom.to_numpy()
gene_starts = gene_df.Start.to_numpy()
gene_ends = gene_df.End.to_numpy()

df['Genes'] = ''
for i, row in df.iterrows():
	chrom = row['Chromosome']
	start = row['Start']
	end = row['End']
	
	
	# Get exons from gencode file that are at least partially contained within the CNV length
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