import pandas as pd

# Apply an intracohort frequency filter
# Variants should be present in <= 10 individuals in the cohort

# Input and output files
input_file='/path/to/input/file.csv' # Use the output of script 7_vcf2table.sh here
output_file='/path/to/output/table.csv' # Note that the output of this script will be the final enhancer SNV annotations

# Load in inputs
df = pd.read_csv(input_file, sep = '\t',
		names = ['Sample', 'Chrom', 'Pos', 'Ref', 'Alt', 'Qual',
				'gg_AF', 'gg_AF_nfe', 'gg_AF_asj', 'gg_AF_eas', 'gg_AF_amr', 'gg_AF_afr',
				'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB'])

# Filter intracohort count
df['variant_id'] = df['Chrom'] + '_' + df['Pos'].astype(str) + '_' + df['Ref'] + '_' + df['Alt']

variant_counts = df['variant_id'].value_counts().to_dict()
df['cohort_count'] = df['variant_id'].map(variant_counts)

df = df[df['cohort_count'] <= 10]

# Save
df.to_csv(output_file, index=False, sep = '\t')
