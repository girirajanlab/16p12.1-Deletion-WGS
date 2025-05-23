import pandas as pd
import numpy as np

# Annotate calls with genes and filter to keep only calls that overlap genes

# Input and output files
input_file="/path/to/input/cnvs.bed" # Use the output of script 6_frequency_filter.py here
output_file="/path/to/output/file.bed"

# Reference files
gencode_genes_path='/path/to/parsed/gencode/genes/file/gencode.v39.parsed.genes.csv'
gencode_exons_path='/path/to/parsed/gencode/exon/file/gencode.v39.parsed.exons.csv'
# These files contains all genes/exons from the Gencode v39 GTF, but parsed into a table format

# Load files
calls=pd.read_csv(input_file)
genes=pd.read_csv(gencode_genes_path)

# Restrict annotations to protein-coding genes
genes = genes[genes.gene_type=='protein_coding']

# Helper functions
def gene_overlap(row):
	chr = row['Chr']
	start = row['Start']
	end = row['End']

	overlap_df = genes[(genes.Chrom==chr) & (((genes.Start >= start) & (genes.Start <= end)) |
			((genes.End >= start) & (genes.End <= end)))]

	if overlap_df.shape[0] == 0:
		return np.nan

	return list(overlap_df.gene_id.unique())

def exon_overlap(row):
	chr = row['Chr']
	start = row['Start']
	end = row['End']

	# Get exon overlaps
	overlap_df = exons[(exons.Chrom==chr) & (((exons.Start >= start) & (exons.Start <= end)) |
				((exons.End >= start) & (exons.End <= end)))]

	if overlap_df.shape[0]==0:
		return np.nan

	# Get overlap of gene ids from gene overlaps and gene ids from exon overlaps
	exon_genes = list(overlap_df.gene_id.unique())
	gene_genes = row['gene_overlap']
	gene_ids = [gene for gene in gene_genes if gene in exon_genes]

	if len(gene_ids) == 0:
		return np.nan

	return gene_ids

# Get overlapping genes
calls['gene_overlap'] = calls.apply(gene_overlap, axis = 1)
calls = calls[~calls.gene_overlap.isnull()]

# Ensure all annotated genes overlap exons
exons=pd.read_csv(gencode_exons_path)
exons = exons[exons.gene_type=='protein_coding']

# Restrict list to only those annotated from the gene set to save time/memory
anno_genes = []
for set in calls.gene_overlap.to_list():
	anno_genes += set
exons = exons[exons.gene_id.isin(anno_genes)]

calls['exon_overlap'] = calls.apply(exon_overlap, axis = 1)
calls = calls[~calls.exon_overlap.isnull()]

# Add gene symbols and IDs
calls['gene_symbols'] = calls['exon_overlap'].apply(lambda row: ' '.join(genes[genes.gene_id.isin(row)]['gene_name'].to_list()))
def remove_dot(ids):
	no_dot = [i.split('.')[0] for i in ids]
	return ' '.join(no_dot)
calls['gene_ids'] = calls['exon_overlap'].apply(remove_dot)

# Remove extra columns
cols = ['Chr', 'Start', 'End', 'Type', 'Zygosity', 'Length', 'NumSNP', 'Sample', 'StartSNP', 'EndSNP', 'Merge',
	'Control_num', 'Control_freq', 'MyCode_num', 'MyCode_freq', 'variant_id', 'Pathogenic_Name', 'gene_ids', 'gene_symbols']
calls=calls[cols]

# Save
calls.to_csv(output_file, sep='\t', index=False)
