import pandas as pd
import sys

# Filter CNVs to remove those with breakpoints that overlap known STR regions

# Input and output files
input_file=sys.argv[1]
output_path=sys.argv[2]

output_file=f'{output_path}/str_all_filter.bed'

# Reference files
# This is a reference file of STRs identified by GangSTR from the Gymrek lab
# It can be downloaded here: https://s3.amazonaws.com/gangstr/hg19/genomewide/hg19_ver13_1.bed.gz
# And is described here: https://github.com/gymreklab/GangSTR
str_file='/path/to/hg19/strs/hg19_ver13_1.bed'

# Read in CNV calls
cnvs = pd.read_csv(input_file, sep = '\t')

# Read in STR regions
strs = pd.read_csv(str_file, sep = '\t', names = ['Chr', 'Start', 'End', 'motif_length', 'motif'])

# Remove any CNV with a breakpoint in an STR region
def str_remove(chr, start, end):
	str_chr = strs[strs.Chr == chr]
	# Get rows in STR file that overlap a breakpoint
	str_breakpoint = str_chr[((str_chr.Start <= start) & (str_chr.End>=start)) | ((str_chr.Start <= end) & (str_chr.End >= end))]
	# If any STRs overlap breakpoints, remove it from output
	if str_breakpoint.shape[0] == 0:
		return True
	else:
		return False
cnvs = cnvs[cnvs.apply(lambda row: str_remove(row['Chr'], row['Start'], row['End']), axis = 1)]

# Save
cnvs.to_csv(output_file, sep = '\t', index = False, header = False)
