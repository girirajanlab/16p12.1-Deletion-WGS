import pandas as pd

# Filter CNV calls for size. Calls should be
# 1. Less than 50kb
# 2. Greater than 99bp

# Input and output files
input_file='/path/to/input_table.txt' # Use the file generated from 7_adjacent_filter.py
output_file='/path/to/output_table.txt'
sample='sample_id'

df = pd.read_csv(input_file, sep = '\t')

# Get the length
df['Size'] = df.End - df.Pos

# Add column for sample name - will be needed later
df['Sample'] = sample

# Remove CNVs >50kb
df = df[df.Size < 50000]
# ALso remove CNVs < 100 bp
df = df[df.Size >= 100]

# Save
df.to_csv(output_file, sep = '\t', index = False)
