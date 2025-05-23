import pandas as pd

# Append ANNOVAR annotations to calls

# Input and output files
ANNOVAR_TABLE="/path/to/annovar/annotated/table.hg38_multianno.txt" # Use the output from script 5_annovar.sh
INPUT_TABLE="/path/to/input/table.tsv" # Use the output from script 4_anno_locus_info.py
OUTPUT_TABLE="/path/to/output/table.tsv"

# Load files
anno=pd.read_csv(ANNOVAR_TABLE, sep='\t')
anno['variant_id'] = anno['Chr'] + '_' + anno['Start'].astype(str) + '_' + anno['Alt']
anno = anno.set_index('variant_id')
anno = anno.drop_duplicates().copy()

df=pd.read_csv(INPUT_TABLE, sep='\t')
df['variant_id'] = df['chrom'] + '_' + df['pos'].astype(str) + '_' + df['alt_allele']

# Apply ANNOVAR annotations
anno_cols = anno.columns[5:]
anno = anno.to_dict()

cols = ['chrom', 'pos', 'sample', 'zscore', 'longest_allele', 'both_alleles',
       'cohort_mean', 'cohort_mode', 'cohort_std', 'end', 'motif',
       'period', 'variant_id']
# Drop ref allele and alt allele
df = df[cols].copy()

for col in anno_cols:
	df[col] = '.'

total = df.shape[0]
for i, row in df.iterrows():
	variant_id = row['variant_id']
	for col in anno_cols:
		annotation = anno[col][variant_id]
		df.at[i, col] = annotation

# Save
df.to_csv(OUTPUT_TABLE, sep='\t', index=False)
