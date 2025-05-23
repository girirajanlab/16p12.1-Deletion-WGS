import pandas as pd

# Apply ANNOVAR annotations to all variants and filter

# Input and output files
ANNOVAR_INPUT="/path/to/annovar/annotated/variants.txt" # Use the output from script 8_annovar.sh here
EXP_INPUT="/path/to/str/expansion/table.tsv" # Use tht output from script 7_2SD_STR.py here
OUTPUT_FILE="/path/to/output/table.tsv"

# Load files
anno=pd.read_csv(ANNOVAR_INPUT, sep='\t')
anno['variant_id'] = anno['Chr'] + '_' + anno['Start'].astype(str) + '_' + anno['Alt']
anno = anno.set_index('variant_id')
anno = anno.drop_duplicates().copy()

df=pd.read_csv(EXP_INPUT, sep='\t')
df['variant_id'] = df['chrom'] + '_' + df['pos'].astype(str) + '_' + df['alt_allele']

# Create a dictionary of annotations
anno_cols = anno.columns[5:]
anno = anno.to_dict()

# Drop ref and alt alleles
cols = ['chrom', 'pos', 'end', 'sample', 'zscore',
       'longest_allele', 'cohort_mode', 'motif', 'motif_period', 'variant_id']
df = df[cols].copy()

# Apply annotations
for col in anno_cols:
	df[col] = '.'

for i, row in df.iterrows():
	variant_id = row['variant_id']
	for col in anno_cols:
		annotation = anno[col][variant_id]
		df.at[i, col] = annotation

# Filter for exonic STRs
df = df[df['Func.wgEncodeGencodeBasicV19'] == 'exonic']

# Save
df.to_csv(OUTPUT_FILE, sep='\t', index=False)

