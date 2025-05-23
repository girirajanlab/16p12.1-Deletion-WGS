import pandas as pd

# Input and output files
input_path='path/to/input/files'
output_path='path/to/output/files'

# Remove CNVs with an inracohort frequency > 10
for cn in ['dels', 'dups']:
	df=pd.read_csv(f'{input_path}/intracohort_{cn}.bed', sep = '\t',
					names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'Intracohort_count'])
	
	# Intracohort filter
	df = df[df.Intracohort_count <= 10]

	# Save
	df.to_csv(f'{output_path}/frequency_filter_{cn}.bed', sep = '\t', index = False)
