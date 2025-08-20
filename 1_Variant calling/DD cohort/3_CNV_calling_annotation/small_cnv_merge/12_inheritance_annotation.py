import pandas as pd
import subprocess
import datetime

# Step to avoid truncation of long columns - relevant for getting filenames
pd.options.display.max_colwidth = 999

# Use the lookup table to find samples that have the same CNV (defined as 50% reciprocal overlap)

# Input and output files
GENE_ANNO="/path/to/gene/annotated/calls.bed" $ Use the output of script 10_annotate_gencode.py
OVERLAP="/path/to/CNV/overlap.bed" # Use the output of script 11_inheritance_lookup.sh
FAM="/path/to/sample/FAM/file.fam"
SAMPLE_LIST="/path/to/list/of/CNV/samples.txt"

CNVNATOR="/path/to/list/of/CNVnator/files.txt"
MANTA="/path/to/list/of/manta/files.txt"
DELLY="/path/to/list/of/delly/files.txt"
LUMPY="/path/to/list/of/lumpy/files.txt"

CNVNATOR_RAW_LOC="/path/to/CNVnator/calls" # Paths to unfiltered calls
MANTA_RAW_LOC="/path/to/manta/calls"
DELLY_RAW_LOC="/path/to/delly/calls"
LUMPY_RAW_LOC="/path/to/lumpy/calls"

CNVNATOR_MERGE_LOC="/path/to/CNVnator/calls" # Paths to calls after merging adjacent calls
MANTA_MERGE_LOC="/path/to/manta/calls"
DELLY_MERGE_LOC="/path/to/delly/calls"
LUMPY_MERGE_LOC="/path/to/lumpy/calls"

OUTPUT_FILE="/path/to/output.bed"

# Load files
calls = pd.read_csv(GENE_ANNO, sep = '\t')
lookup = pd.read_csv(OVERLAP, sep = '\t',
	names = ['Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Sample', 'Intracohort_count', 'gnomADSV_AF', 'gnomADSV_AFfilter', 'gene_ids', 'gene_names', 'variant_id',
		'Chr2', 'Start2', 'End2', 'Type2', 'Name2', 'Length2', 'Sample2', 'Intracohort_count2', 'gnomADSV_AF2', 'gnomADSV_AFfilter2', 'gene_ids2', 'gene_names2', 'variant_id2'])
fam = pd.read_csv(FAM, sep = '\t', names = ['Family', 'Sample', 'Father', 'Mother', 'Sex', 'Carrier_Status'])

cnv_samps = pd.read_csv(SAMPLE_LIST, names = ['Sample'])

# Caller output file locations
cnvnator = pd.read_csv(CNVNATOR)
manta =  pd.read_csv(MANTA)
delly = pd.read_csv(DELLY)
lumpy = pd.read_csv(LUMPY)

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

	command = "bcftools view -r %s:%i-%i -t %s:%i-%i -i 'INFO/SVTYPE=\"%s\"' %s -H" % (chr, search_start, search_end, chr, search_start, search_end, svtype, file)
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

# Look through merged files from each caller for overlap
def mergefile_lookup(file,chr, start, end, length, svtype):
	# Read in parent file as dataframe
	df = pd.read_csv(file, sep = '\t')
	# Unify column names
	df.columns = ['Chr', 'Pos', 'End', 'Type', 'ID']
	# Unify SV type name
	type_dict = {'<DEL>':'DEL', '<DUP>':'DUP', 'DEL':'DEL', 'DUP':'DUP'}
	df['Type'] = df.Type.map(type_dict)

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
		return 'both'
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

	caller_list_files = [cnvnator, manta, delly, lumpy]
	caller_locs = [CNVNATOR_RAW_LOC, MANTA_RAW_LOC, DELLY_RAW_LOC, LUMPY_RAW_LOC]
	merge_locs = [CNVNATOR_MERGE_LOC, MANTA_MERGE_LOC, DELLY_MERGE_LOC, LUMPY_MERGE_LOC]
	if mother_code != '0':
		for i, caller in enumerate(caller_list_files):
			if mother_code in caller.Sample.to_list():
				loc = caller_locs[i]
				# Filtered outputs
				mother_file = f"{caller_loc}/{mother_code}.filtered.vcf.gz"
				mother = rawfile_lookup(mother_file, chr, start, end, length, type)

				if mother:
					break

				# Merged outputs
				merge_loc=merge_locs[i]
				mother_merge = f"{merge_loc}/{mother_code}_final.bed"
				mother = mergefile_lookup(mother_merge, chr, start, end, length, type)

				if mother:
					break


	if father_code != '0':
		for i, caller in enumerate(caller_list_files):
			if father_code in caller.Sample.to_list():
				loc = caller_locs[i]
				# Filtered outputs
				father_file = f"{caller_loc}/{father_code}.filtered.vcf.gz"
				father = rawfile_lookup(father_file, chr, start, end, length, type)

				if mother:
					break

				# Merged outputs
				merge_loc=merge_locs[i]
				father_merge = f"{merge_loc}/{father_code}_final.bed"
				father = mergefile_lookup(father_merge, chr, start, end, length, type)

				if mother:
					break

	if mother and father:
		return 'both'
	elif mother:
		return mother_return
	elif father:
		return father_return
	elif mother_code=='0' or father_code=='0':
		return '.'
	elif mother_code not in cnv_samps.Sample.to_list() or father_code not in cnv_samps.Sample.to_list():
		return '.'
	else:
		return 'de novo'

calls['inheritance'] = calls.apply(inheritance, axis = 1)

calls.to_csv(OUTPUT_FILE, sep = '\t', index = False)
