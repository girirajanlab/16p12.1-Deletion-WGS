import pandas as pd

# Generate tables of genes preferrentially expressed in specific brain cell types

# Input and Output files
MEDIANS="/path/to/Huuman/M1/10x/medians.csv" # Cell expression medians downloaded from https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/medians.csv
METADATA="/path/to/Huuman/M1/10x/metadata.csv" # Metadata downloaded from https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv
GENE_ANNO="/path/to/gene/annotations.csv" # Use the output of script 2_Analysis preparation/Gene_Annotations/5_add_loeuf.py

SUBCLASS_OUT="/path/to/output/dataframe/subclass/preferrentially/expressed/genes.csv"
SUPERCLASS_OUT="/path/to/output/dataframe/superclass/preferrentially/expressed/genes.csv"
COMBINED_OUT="/paht/to/output/dataframe/combined/preferentially/expressed/genes.csv"

# Load data
cell_df = pd.read_csv(METADATA)
df = pd.read_csv(MEDIANS)

# Parse cell class labels
cell_df = cell_df[['cluster_label', 'subclass_label']]
cell_df['superclass_label'] = cell_df["cluster_label"].str.split(" ", 1).str[0]
cell_df = cell_df.drop_duplicates()
cell_df = cell_df.set_index('cluster_label', drop=True)

# Add class labels to data
cell_type_labels = list(df.columns[1:])
df = df.set_index('feature', drop=True)
df.columns = cell_type_labels

df=pd.concat([df, cell_df[["subclass_label"]].transpose()])
df=pd.concat([df, cell_df[["superclass_label"]].transpose()])

superclass_list=list(set(cell_df['superclass_label'].values.tolist()))
subclass_list=list(set(cell_df['subclass_label'].values.tolist()))

# Calculate summed read counts for each class type
subclass_df=pd.DataFrame()
for subclass in subclass_list:
	sub_df = df.loc[:, df.loc['subclass_label'] == subclass]
	sub_df = sub_df.drop(['superclass_label','subclass_label'])
	gene_sum = sub_df.sum(axis=1)
	subclass_df[subclass]=gene_sum

superclass_df=pd.DataFrame()
for superclass in superclass_list:
	subdf = df.loc[:, df.loc['superclass_label'] == superclass]
	subdf = subdf.drop(['superclass_label','subclass_label'])
	gene_sum = subdf.sum(axis=1)
	superclass_df[superclass]=gene_sum

# Calculate medians and standard deviations of genes
sub_gene_medians = subclass_df.median(axis=1)
sub_gene_std_dev = subclass_df.std(axis=1)

sup_gene_medians = superclass_df.median(axis=1)
sup_gene_std_dev = superclass_df.std(axis=1)

# Create a frame annotated with cell tye and whether a gene is preferrentially expressed in that cell type and save

# Subclass
new_df = pd.DataFrame(index=subclass_df.index)
for cell_type in subclass_df.columns:
	new_df[cell_type] = subclass_df[cell_type] > (sub_gene_medians + 2*sub_gene_std_dev)
	new_df[cell_type] = new_df[cell_type].astype(int)
new_df.to_csv(SUBCLASS_OUT, index=True)

# Superclass
new_df = pd.DataFrame(index=subclass_df.index)
for cell_type in subclass_df.columns:
	new_df[cell_type] = superclass_df[cell_type] > (sup_gene_medians + 2*sup_gene_std_dev)
	new_df[cell_type] = new_df[cell_type].astype(int)
new_df.to_csv(SUPERCLASS_OUT, index=True)

# Merge data from both sub- and super-classes into single file
sub=pd.read_csv(SUBCLASS_OUT, index_col=0)
super=pd.read_csv(SUPERCLASS_OUT, index_col=0)

df=pd.merge(sub, super, right_index=True, left_index=True, suffixes=['_subclass', '_superclass'])

genedf=pd.read_csv(GENE_ANNO)

df=pd.merge(df, genedf[['gene_symbol', 'gene_id']], left_index=True, right_on='gene_symbol', how='inner')

# Save
df.to_csv(COMBINED_OUT, index=False)