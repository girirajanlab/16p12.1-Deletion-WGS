import pandas as pd
import os

# Remove calls that are:
# 1. < 50kb, or
# 2. Contain < 5 SNPs

# Input and output files
input_file="/path/to/input/CNV/calls.txt" # Use the output from script 3_merge_calls.py
output_path="/path/to/output/directory"

# Combine all calls into a single file
df = pd.read_csv(input_file, sep='\t')

# Size filter
df = df[df.Length >= 50000]
df = df[df.NumSNP >= 5]

# Separate deletions and duplications for downstream filtering
dels = df[df.Type=='DEL']
dups = df[df.Type=='DUP']

# Save
dels.to_csv(f'{output_dir}/dels.bed', sep = '\t', index = False)
dups.to_csv(f'{output_dir}/dups.bed', sep = '\t', index = False)

