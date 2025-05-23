import pandas as pd

# Apply an intracohort frequency filter
# Variants should be present in <= 10 individuals in the cohort

# Input and output files
input_file='/path/to/input/file.csv' # Use the output of script 10_filter_variant_types.py here
output_file='/path/to/output/table.csv'

df = pd.read_csv(input_file)

# Filter for intracohort frequency
df['variant_id'] = df['Chrom'] + '_' + df['Pos'].astype(str) + '_' + df['Ref'] + '_' + df['Alt']
variant_counts = df['variant_id'].value_counts()
df['cohort_count'] = df['variant_id'].map(variant_counts)

df = df[df['cohort_count'] <= 10]

df.to_csv(output_file, index=False)
