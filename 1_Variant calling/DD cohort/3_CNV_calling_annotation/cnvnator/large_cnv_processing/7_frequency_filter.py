import pandas as pd

# Input and output files
input_path='path/to/input/files'
output_path='path/to/output/files'

# Remove CNVs with
# A microarray control frequency > 0.001 (if not a known pathogenic CNV)
# or, an inracohort frequency > 10
for cnv_type in ['cnv', 'pathogenic_cnv']:
	for cn in ['dels', 'dups']:
		df=pd.read_csv(f'{input_path}/{cnv_type}_control_intracohort_{cn}.bed', sep = '\t',
						names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'Microarray_count', 'Intracohort_count'])

		# Control population filter
		# 18572 unique individuals in file
		df['microarray_freq'] = df.Microarray_count/18572
		if cnv_type=='cnv':
			df = df[df.microarray_freq <= 0.001]
		
		# Intracohort filter
		df = df[df.Intracohort_count <= 10]

		# Save
		df.to_csv(f'{output_path}/{cnv_type}_frequency_filter_{cn}.bed', sep = '\t', index = False)

