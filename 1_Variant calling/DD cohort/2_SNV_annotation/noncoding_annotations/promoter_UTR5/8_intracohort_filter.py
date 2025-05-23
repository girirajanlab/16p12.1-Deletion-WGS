import pandas as pd

# Apply an intracohort frequency filter
# Variants should be present in <= 10 individuals in the cohort

# Input and output files
input_file='/path/to/input/file.txt' # Use the output of script 7_format_table.sh here
output_promoter_file='/path/to/output/promoter_table.txt' # Note that the output of this script will be the final promoter and 5' UTR SNV annotations
output_UTR5_file='/path/to/output/UTR5_table.txt'

# Load input table
df = pd.read_csv(input_file, sep = '\t')

# Apply intracohort filter
df['variant_id'] = df['Chrom'] + '_' + df['Pos'].astype(str) + '_' + df['Ref'] + '_' + df['Alt']
variant_counts = df['variant_id'].value_counts()
df['cohort_count'] = df['variant_id'].map(variant_counts)

df = df[df['cohort_count'] <= 10]

# Removing noncoding RNA  5' UTR variants
df = df[df['Func.wgEncodeGencodeBasicV19']!='ncRNA_UTR5']

# Separate promoter and 5' UTR variants and save
df[df['Func.wgEncodeGencodeBasicV19'] == 'upstream'].to_csv(output_promoter_table, sep = '\t', index = False)
df[df['Func.wgEncodeGencodeBasicV19'] == 'UTR5'].to_csv(output_UTR5_table, sep = '\t', index = False)
