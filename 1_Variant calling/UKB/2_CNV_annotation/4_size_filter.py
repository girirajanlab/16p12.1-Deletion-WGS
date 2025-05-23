import pandas as pd
import os

# Remove calls that are:
# 1. < 50kb, or
# 2. Contain < 5 SNPs

# Input and output files
input_path="/path/to/input/CNV/calls" # Use the output directory from script 3_merge_calls.py
output_path="/path/to/output/directory"

# Combine all calls into a single file
out_df = pd.DataFrame(columns = ['Chr', 'Start', 'End', 'Type', 'Zygosity', 'Length', 'NumSNP', 'Sample', 'StartSNP', 'EndSNP', 'Merge'])
for (root,dirs,files) in os.walk(input_path, topdown=True):
	for f in files:
		df=pd.read_csv(f, sep='\t')
		df = df[df.Length >= 50000]
		df = df[df.NumSNP >= 5]
		out_df=pd.concat([out_df, df], axis=0)

# Separate deletions and duplications for downstream filtering
dels = out_df[out_df.Type=='DEL']
dups = out_df[out_df.Type=='DUP']

# Save
dels.to_csv(f'{output_dir}/dels.bed', sep = '\t', index = False)
dups.to_csv(f'{output_dir}/dups.bed', sep = '\t', index = False)

