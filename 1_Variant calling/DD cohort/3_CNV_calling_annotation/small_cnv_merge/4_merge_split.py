import pandas as pd

# Merge files from all samples
# Then split into separate files for deletions and duplications

# Input and output files
input_files='/path/to/list/of/input/files.csv' # Use files generated from script 3_merge_calls.py
output_path='/path/to/output/files/'

# Merge all calls into a single file
filenames=open(input_files, 'r').readlines()
df=pd.DataFrame()
for f in filenames:
	fdf=pd.read_csv(f, sep='\t')
	df=pd.concat([df, fdf])

# Split CNVs into Dels and Dups
del_df = df[df.Type=='DEL']
dup_df = df[df.Type=='DUP']

# Save
del_df.to_csv(f'{output_path}/dels.bed', sep = '\t', index = False)
dup_df.to_csv(f'{output_path}/dups.bed', sep = '\t', index = False)
