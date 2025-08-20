import pandas as pd
import subprocess
import datetime

# Annotate inheritance of CNVs

# Step to avoid truncation of long columns - relevant for getting filenames
pd.options.display.max_colwidth = 999

# Input and output files
INPUT_CALLS="/path/to/input/calls.bed" # Use the output of script 12_finalize_calls.py
OVERLAP="/path/to/overlapped/calls.bed" # Use the output of script 13_inhertiance_lookup.sh
FAM="/path/to/sample/FAM/file.fam"
SAMPLE_LIST="/path/to/list/of/CNV/samples.txt"

PENNCNV_CALLS="/path/to/PennCNV/calls.csv" # Use the output of script 3_CNV_calling_annotation\PennCNV\10_annotate_gencode.py
CNVNATOR_RAW_LOC="/path/to/CNVnator/calls" # Paths to unfiltered calls
CNVNATOR_MERGE_LOC="/path/to/CNVnator/calls" # Paths to calls after merging adjacent calls

OUTPUT_FILE="/path/to/output.bed"

# Use the lookup table to find samples that have the same CNV (defined as 50% reciprocal overlap)
calls = pd.read_csv(INPUT_CALLS, sep = '\t')
lookup = pd.read_csv(OVERLAP, sep = '\t',
	names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'Microarray_count', 'Intracohort_count', 'microarray_freq', 'gnomADSV_AF', 'gnomADSV_AFfilter',
		'gene_ids', 'gene_names', 'NEJM_Name', 'variant_id',
		'Chr2', 'Start2', 'End2', 'Type2', 'Name2', 'Length2', 'Sample2', 'Microarray_count2', 'Intracohort_count2', 'microarray_freq2', 'gnomADSV_AF2', 'gnomADSV_AFfilter2',
		'gene_ids2', 'gene_names2', 'NEJM_Name2', 'variant_id2'])

fam = pd.read_csv(FAM, sep = '\t', names = ['Family', 'Sample', 'Father', 'Mother', 'Sex', 'Carrier_Status'])
cnv_samps = pd.read_csv(SAMPLE_LIST)
penncnv = pd.read_csv(PENNCNV_CALLS)

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

# Use bcftools to search through the merged output from the callers to search for calls that didn't pass QC
def rawfile_lookup(file, chr, start, end, length, svtype):
	# Due to 50% reciprocal overlap requirement, any overlapping CNV can't start more than LENGTH away from the START of the query CNV
	# It also can't start any more than half way through the query CNV
	search_start = start-length
	search_end = start+(length/2)

	# Reformat type to match VCFs
	type = svtype.replace('<', "").replace('>', "")

	command = "bcftools view -r %s:%i-%i -t %s:%i-%i -i 'INFO/SVTYPE=\"%s\"' %s -H" % (chr, search_start, search_end, chr, search_start, search_end, type, file)
	match_lines = subprocess.run(command, capture_output = True, shell = True).stdout.decode()

	if len(match_lines) == 0:
		return False

	match_list = [i.split('\t') for i in match_lines.split('\n') if i!='']
	# Convert to df
	match_df = pd.DataFrame(match_list, columns = ['Chr','Pos', 'ID', 'Ref', 'Alt', 'Qual', 'FILTER', 'INFO', 'FORMAT', 'FORMAT_VALUES'])
	# Separate out end and length
	match_df['End'] = match_df.INFO.str.split('END=', expand = True)[1].str.split(';', expand = True)[0]
	match_df['Pos'] = pd.to_numeric(match_df['Pos'])
	match_df['End'] = pd.to_numeric(match_df['End'])

	match_df['Length'] = match_df.End - match_df.Pos

	# Get only calls that meet 50% reciprocal overlap
	match_df = match_df[match_df.apply(lambda row: recip_overlap(start, end, row['Pos'], row['End'], length, row['Length']), axis = 1)]

	if match_df.shape[0] > 0:
		return True

	return False

# Look through merged files for overlap
def mergefile_lookup(file, chr, start, end, length, svtype):
	# Read in parent file as dataframe
	df = pd.read_csv(file, sep = '\t')
	# Unify column names
	df.columns = ['Chr', 'Pos', 'End', 'Type', 'ID']

	# Filter line
	search_start = start-length
	search_end = start+(length/2)
	df = df[(df.Chr==chr) & (df.Type==svtype) & (df.Pos>=search_start) & (df.Pos<=search_end)]

	# Add length
	df['Length'] = df.End - df.Pos

	# Get 50% reciprocal overlap calls
	df = df[df.apply(lambda row: recip_overlap(start, end, row['Pos'], row['End'], length, row['Length']), axis = 1)]

	if df.shape[0] > 0:
		return True

	return False

# Look through PennCNV output
def penncnv_lookup(sample, chr, start, end, length, svtype):
	# Get lines for parent from PennCNV calls
	df = penncnv[penncnv.Sample==sample]
	if df.shape[0]==0:
		return False
	# Get correct chromosome and type
	# Change type because PennCNV calls have type differently
	new_type = svtype.replace('<', "").replace('>', "").lower()
	df = df[(df.Chromosome==chr) & (df.Type==new_type)]
	if df.shape[0]==0:
		return False
	# Get correct location
	search_start = start-length
	search_end = start+(length/2)
	df = df[(df.Start >= search_start) & (df.Start <= search_end)]
	if df.shape[0]==0:
		return False
	# Check for 50% reciprocal overlap
	df = df[df.apply(lambda row: recip_overlap(start, end, row['Start'], row['End'], length, row['Length']), axis = 1)]
	if df.shape[0]>0:
		return True
	return False

def inheritance(row, father = False, mother = False):
	iid = row['Sample']

	# Get family information from FAM file
	mother_code = fam[fam.Sample==iid]['Mother'].to_string(index = False, header = False).strip()
	father_code = fam[fam.Sample==iid]['Father'].to_string(index = False, header = False).strip()

	# Check if we have CNV calls for the parents
	if mother_code == '0' and father_code == '0':
		return '.'
	if mother_code not in cnv_samps.Sample.to_list() and father_code not in cnv_samps.Sample.to_list():
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
		return '.'
	elif mother:
		return mother_return
	elif father:
		return father_return

	# If CNV wasn't found in the final output for either person, look into the output files for each CNV caller for each person
	chr = row['Chr']
	start = row['Start']
	end = row['End']
	length = row['Length']
	type = row['Type']

	if mother_code != '0' and mother_code in cnv_samps.Sample.to_list():
		# Filtered outputs
		mother_file = f"{CNVNATOR_RAW_LOC}/{mother_code}.filtered.vcf.gz"
		mother = rawfile_lookup(mother_file, chr, start, end, length, type)

		if not mother:
			# Merged outputs
			mother_merge = f"{CNVNATOR_MERGE_LOC}/{mother_code}.filtered.vcf.gz"
			mother = mergefile_lookup(mother_merge, chr, start, end, length, type)

		if not mother:
			# PennCNV outputs
			if mother_code in penncnv.Sample.to_list():
				mother = penncnv_lookup(mother_code, chr, start, end, length, type)

	if father_code != '0' and father_code in cnv_samps.Sample.to_list():
		# Filtered outputs
		father_file = f"{CNVNATOR_RAW_LOC}/{father_code}.filtered.vcf.gz"
		father = rawfile_lookup(father_file, chr, start, end, length, type)

		if not father:
			# Merged outputs
			father_merge = f"{CNVNATOR_MERGE_LOC}/{father_code}.filtered.vcf.gz"
			father = mergefile_lookup(father_merge, chr, start, end, length, type)

		if not father:
			# PennCNV outputs
			if father_code in penncnv.Sample.to_list():
				father = penncnv_lookup(father_code, chr, start, end, length, type)

	if mother and father:
		return '.'
	elif mother:
		return mother_return
	elif father:
		return father_return
	elif mother_code=='0' or father_code=='0':
		return '.'
	elif mother_code not in cnv_samps.Sample.to_list() or father_code not in cnv_samps.Sample.to_list():
		return '.'
	else:
		return '.'

calls['inheritance'] = calls.apply(inheritance, axis = 1)

calls.to_csv(OUTPUT_FILE, sep = '\t', index = False)
