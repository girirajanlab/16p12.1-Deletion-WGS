import pandas as pd

# Split CNVs into two groups by size:
# 1. 50kb or greater
# 2. Less than 50kb and greater than 100bp
# Also merge calls across samples

# Input and output files
input_files='/path/to/list/of/input/filenames.txt' # Use the files generated from 2_adjacent_filter.py
output_path='/path/to/output/files'

output_large=output_path+'/large_cnvnator.bed'
output_small=output_path+'/small_cnvnator.bed'

# Iterate through filenames to create individual tables with calls from all samples
filenames=open(input_files, 'r').readlines()

# Initialize dataframes
ldf=pd.DataFrame()
sdf=pd.DataFrame()

for f in filenames:
	# Read in merged variants
	df = pd.read_csv(input_file, sep = '\t')

	# Get the length
	df['Size'] = df.End - df.Pos

	# Add column for sample name - will be needed later
	df['Sample'] = sample

	#Separate >=50kb CNVs and <50kb CNVs
	large_df = df[df.Size >= 50000]
	small_df = df[df.Size < 50000]
	# Also remove CNVs < 100bp
	small_df = small_df[small_df.Size >= 100]

	ldf=pd.concat([ldf, large_df])
	sdf=pd.concat([sdf, small_df])

# Save
ldf.to_csv(output_large, sep = '\t', index = False)
sdf.to_csv(output_small, sep = '\t', index = False)
