import pandas as pd

# Filter non-pathogenic CNV calls for those with a frequency in gnomADSV <= 0.001

# Input and output files
lookup_path='path/to/lookup/files' # These are the files from script 7_gnomadSV_annotation.sh
cnv_path='path/to/cnv/files' # These are the files from script 5_intracohort_filter.py
output_path='path/to/output/files'

def get_AF(s, lookup):
	# There are some CNVs that have 50% reciprocal overlap with more than one SV in gnomAD
	# Base filter on the highest AF, but report both
	chr = s.Chr
	start = s.Start
	end = s.End
	# Find lines in the lookup table with the same location
	overlap_lines = lookup[(lookup.Chr==chr) & (lookup.Start==start) & (lookup.End==end)]
	# Get allele frequencies
	afs = list(overlap_lines.AF.unique())
	afs.sort()
	s['gnomADSV_AF'] = ' '.join(afs)

	if afs==['.']:
		s['gnomADSV_AFfilter'] = True
	elif max([float(i) for i in afs if i!='.']) <= 0.001:
		s['gnomADSV_AFfilter'] = True
	else:
		s['gnomADSV_AFfilter'] = False

	return s

for cn in ['dels', 'dups']:
	lookup=pd.read_csv(f'{lookup_path}/gnomadSV_anno_{cn}.bed', sep='\t')
	df=pd.read_csv(f'{cnv_path}/frequency_filter_{cn}.bed', sep='\t')

	# Get AF and apply filter (if not pathogenic CNVs)
	df = df.apply(lambda x: get_AF(x, lookup), axis = 1)
	df = df[df.gnomADSV_AFfilter]
	
	# Save
	df.to_csv(f'{output_path}/gnomadSV_filter_{cn}.bed', index = False, sep = '\t')
