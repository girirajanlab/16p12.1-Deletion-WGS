import pandas as pd
import subprocess
import datetime

# Annotate inheritance of CNVs

# Step to avoid truncation of long columns - relevant for getting filenames
pd.options.display.max_colwidth = 999

# Input and output files
INPUT_CALLS="/path/to/input/calls.bed" # Use the output of script 10_annotate_gencode.py
OVERLAP="/path/to/overlapped/calls.bed" # Use the output of script 11_inhertiance_lookup.sh
FAM="/path/to/sample/FAM/file.fam"
SAMPLE_LIST="/path/to/list/of/CNV/samples.txt"

OUTPUT_FILE="/path/to/output.bed"

# Use the lookup table to find samples that have the same CNV (defined as 50% reciprocal overlap)
calls = pd.read_csv(INPUT_CALLS, sep = '\t')
lookup = pd.read_csv(OVERLAP, sep = '\t',
	names = ['Chr', 'Start', 'End', 'Type', 'Sample', 'variant_id', 'Chr2', 'Start2', 'End2', 'Type2', 'Sample2', 'variant_id2'])

fam = pd.read_csv(FAM, sep = '\t', names = ['Family', 'Sample', 'Father', 'Mother', 'Sex', 'Carrier_Status'])
penncnv = pd.read_csv('../../2022_02_14/PennCNV_calls_combined_final.csv')

# Find 50% reciprocal overlap calls in the final output and return the samples they belong to
def overlap_samps(vid):
	overlaps = lookup[(lookup.variant_id==vid) & (lookup.variant_id2!=vid)]
	samps = list(overlaps.Sample2.unique())

	return samps

# Function to check for reciprocal overlap
def recip_overlap(start1, end1, start2, end2, length1, length2, threshold = 0.5):
	overlap = (min(end1, end2) - max(start1, start2))
	if overlap >= (threshold*max(length1, length2)):
		return True
	return False

def inheritance(row, father = False, mother = False):
	iid = row['PatientID']

	# Get family information from FAM file
	mother_code = fam[fam.Sample==iid]['Mother'].to_string(index = False, header = False).strip()
	father_code = fam[fam.Sample==iid]['Father'].to_string(index = False, header = False).strip()

	# Check if we have CNV calls for the parents
	if mother_code == '0' and father_code == '0':
		return '.'
	if mother_code not in penncnv.PatientID.to_list() and father_code not in penncnv.PatientID.to_list():
		return '.'

	# Get carrier status of parents
	if mother_code == '0':
		mother_return = '0'
	elif fam[fam.Sample==mother_code]['Carrier_Status'].to_string(index = False).strip()=='unknown':
		mother_return = 'M'
	elif int(fam[fam.Sample==mother_code]['Carrier_Status'])==1:
		mother_return = 'MC'
	elif int(fam[fam.Sample==mother_code]['Carrier_Status'])==0:
		mother_return = 'MNC'
	else:
		mother_return = 'M'

	if father_code == '0':
		father_return = '0'
	elif fam[fam.Sample==father_code]['Carrier_Status'].to_string(index = False).strip()=='unknown':
		father_return = 'F'
	elif int(fam[fam.Sample==father_code]['Carrier_Status'])==1:
		father_return = 'FC'
	elif int(fam[fam.Sample==father_code]['Carrier_Status'])==0:
		father_return = 'FNC'
	else:
		father_return = 'F'

	# Check if parents have the same variant
	other_samps = overlap_samps(row['variant_id'])

	if mother_code in other_samps:
		mother = True
	if father_code in other_samps:
		father = True

	# Return inheritance
	if mother and father:
		return 'both'
	elif mother:
		return mother_return
	elif father:
		return father_return
	else:
		return '.'

calls['inheritance'] = calls.apply(inheritance, axis = 1)

calls.to_csv(OUTPUT_FILE, sep = '\t', index = False)
