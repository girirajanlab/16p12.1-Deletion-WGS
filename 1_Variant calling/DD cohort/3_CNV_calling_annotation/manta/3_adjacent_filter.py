import pandas as pd
import sys

# Merge adjacent calls per individual
# Merge calls if they:
# 1. Overlap (any %)
# 2. Have a gap of less than 50kb AND <20% of the combined CNV length

# Input and output files
input_file='/path/to/input_table.txt' # Use the file generated from 2_get_cnvs.py
intermediate_output_dir='/path/to/directory/for/intermediate/files'
output_file='/path/to/output_table.txt'

# Read file as dataframe
vcf = pd.read_csv(input_file, sep = '\t')
# Make a new column for END for easier filtering
vcf['End'] = vcf.INFO.str.split('END=', expand = True)[1].str.split(';', expand = True)[0]
vcf['End'] = pd.to_numeric(vcf['End'])

# Functions
# Format IDs
def get_ids(id, rows):
	ids = [id.split(';')] + rows.ID.str.split(';').tolist()
	ids = [item for items in ids for item in items if items]
	ids.sort()
	ids2 = ';'.join(list(set(ids)))
	return(ids2)

# Function to read VCF line by line and do merging
def merge_overlap_gap(df):
	to_skip = []
	df2 = pd.DataFrame(columns = ['Chr', 'Pos', 'End', 'Type', 'ID'])
	for idx, row in df.iterrows():

		# Skip line, if needed
		if idx in to_skip:
			continue

		chr = row['Chr']
		start = row['Pos']
		end = row['End']
		type = row['Type']
		id = row['ID']


		# Check for overlap
		# 1. Complete overlap
		# One CNV is completely surrounded by another
		#   <---->
		# <--------->
		# Current CNV is completely covered - skip and only annotate the larger CNV
		co1_rows = df[(df.Type==type) & (df.Chr==chr) & (df.Pos<=start) & (df.End>=end) & (df.ID!=id)]
		if co1_rows.shape[0] > 0:
			continue
		# Another CNV is completely covered - skip other CNV
		co2_rows = df[(df.Type==type) & (df.Chr==chr) & (start<=df.Pos) & (end>=df.End) & (df.ID!=id)]
		if co2_rows.shape[0] > 0:
			to_skip = to_skip + co2_rows.index.to_list()

		# 2. Incomplete overlap
		#  <---->
		#     <---->
		# Current CNV is first
		io1_rows = df[(df.Type==type) & (df.Chr==chr) & (start<=df.Pos) & (end>=df.Pos) & (end<=df.End) & (df.ID!=id)]
		if io1_rows.shape[0] > 0:
			# Remove rows with a smaller index than the current row - these would already have been merged
			io1_rows = io1_rows[io1_rows.index > idx]
			if io1_rows.shape[0] > 0:
				# Write merged call to a new df
				ids = get_ids(id, io1_rows)
				df2.loc[len(df2.index)] = [chr, min(io1_rows.Pos.to_list()+[start]), max(io1_rows.End.to_list()+[start]), type, ids]
				continue
		# Current CNV is second
		io2_rows = df[(df.Type==type) & (df.Chr==chr) & (df.Pos<=start) & (df.End>=start) & (df.End<=end) & (df.ID!=id)]
		if io2_rows.shape[0] > 0:
			# Remove rows with a smaller index than the current row - these would already have been merged
			io2_rows = io2_rows[io2_rows.index > idx]
			if io2_rows.shape[0] > 0:
				# Write merged call to a new df
				ids = get_ids(id, io2_rows)
				df2.loc[len(df2.index)] = [chr, min(io2_rows.Pos.to_list()+[start]), max(io2_rows.End.to_list()+[start]), type, ids]
				continue

		# Check for small gaps
		# <---> <--->
		# Current CNV is first
		sg1_rows = df[(df.Type==type) & (df.Chr==chr) & (df.Pos>end) & (df.Pos<end+50000) & ((df.Pos-end)<((df.End-start)*0.2)) & (df.ID!=id)]
		if sg1_rows.shape[0] > 0:
			# Remove rows with a smaller index than the current row - these would already have been merged
			sg1_rows = sg1_rows[sg1_rows.index > idx]
			if sg1_rows.shape[0] > 0:
				# Write merged call to new df
				ids = get_ids(id, sg1_rows)
				df2.loc[len(df2.index)] = [chr, min(sg1_rows.Pos.to_list()+[start]), max(sg1_rows.End.to_list()+[start]), type, ids]
				continue
		# Current CNV is second
		sg2_rows = df[(df.Type==type) & (df.Chr==chr) & (start>df.End) & (start<df.End+50000) & ((start-df.End)<((end-df.Pos)*0.2)) & (df.ID!=id)]
		if sg2_rows.shape[0] > 0:
			# Remove rows with a smaller index than the current row - these would already have been merged
			sg2_rows = sg2_rows[sg2_rows.index > idx]
			if sg2_rows.shape[0] > 0:
				# Write merged call to new df
				ids = get_ids(id, sg2_rows)
				df2.loc[len(df2.index)] = [chr, min(sg2_rows.Pos.to_list()+[start]), max(sg2_rows.End.to_list()+[start]), type, ids]
				continue

		# If it doesn't merge, just write the current CNV
		df2.loc[len(df2.index)] = [chr, start, end, type, id]

	return(df2)

# Run once
merge1 = merge_overlap_gap(vcf)

# Run iteratively until calls stop merging
df1 = vcf
df2 = merge1
counter = 1
while not df1.equals(df2):
	df2.to_csv(f'{intermediate_output_dir}/{str(counter)}.txt', index = False, sep = '\t')
	new_df = merge_overlap_gap(df2)
	counter+=1
	df1 = df2
	df2 = new_df

df2.to_csv(output_file, index = False, sep = '\t')
