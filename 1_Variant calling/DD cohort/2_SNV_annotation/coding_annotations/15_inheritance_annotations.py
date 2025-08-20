import pandas as pd
import subprocess
import datetime

# Step to avoid truncation of long columns - relevant for getting filenames
pd.options.display.max_colwidth = 999

# Input and output files
variant_file="/path/to/input/variants.csv" # Use the output of script 14_loeuf_scores.py here
fam_file="/path/to/sample/FAM/file.fam"
vcf_files="/path/to/vcf/locations/file.csv" # File should contain the sample name and path to the output gVCF from 1_GATK/3_split_norm.sh
output_table="/path/to/output/table/with/inheritance.txt"

# Load files
var_file = pd.read_csv(variant_file)
ped = pd.read_csv(fam_file, sep = '\t', names = ['Family', 'Sample', 'Father', 'Mother', 'Sex', 'Carrier_Status'])
vcfs = pd.read_csv(vcf_files)

def find_in_vcf(vcf, chr, pos, alt):
	# Use bcftools to get relevant line in parent file
	command = "bcftools view -r %s:%s -H %s" % (chr, pos, vcf)
	rel_line = subprocess.run(command, capture_output = True, shell = True).stdout.decode()

	# If no matching line is found - should not happen as these are gVCFs, but just in case
	if len(rel_line) == 0:
		print('No matching line found!')
		print(chr, pos, vcf)
		return [-1, -1]

	line = rel_line.split('\t')
	format = line[8].split(':')
	info = line[9].split(':')
	# Look for allele depth information
	# NOTE: AD will only be reported in cases where parent has non-reference allele
	# Check if 'AD' is reported
	if 'AD' in format:
		ad_loc = format.index('AD')
		ad = info[ad_loc]
		ad = [int(i) for i in ad.split(',')]
		# Be careful of multi-allelic sites!
		# Check the length of AD
		# If it is more than 2, only return the read depth information for the reference allele and specific alternate allele
		# Alleles are always in the order Ref, Alt1, Alt2, Alt3, . . . , <NON_REF>
		if len(ad) > 2:
			alts = line[4].split(',')
			# Check if the alternate allele is in the list
			# If not, return the AD for <NON_REF> (should be 0)
			if alt in alts:
				alt_loc = alts.index(alt)+1
			else:
				print(alts)
				print(alt)
				alt_loc = alts.index('<NON_REF>')+1
			ad = [ad[0], ad[alt_loc]]
		return ad

	# If AD is not in parent file, they only have reference allele
	# Confirm this is true by checking genotype info
	geno_loc = format.index('GT')
	gt = info[geno_loc]

	if gt != '0/0':
		print('Parent as non-reference allele but is missing AD information!')
		return [-1, -1]
	# If they are homozygous for the reference, return an alt allele depth of 0 and the read depth as the depth of the reference
	dp_loc = format.index('DP')
	dp = info[dp_loc]
	return [int(dp), 0]


def add_inherit(row, father = False, mother = False, unknown = False):
	if int(row.name)%500 == 0:
		print('Started line ' + str(row.name) + " at: " + str(datetime.datetime.now()))

	iid = row['Sample']

	# Get parents
	father_code = ped[ped.Sample==iid]['Father'].to_string(index = False, header = False).strip()
	mother_code = ped[ped.Sample==iid]['Mother'].to_string(index = False, header = False).strip()

	# Check if we have WGS data available for the parents
	if mother_code not in vcfs.Sample.to_list() and father_code not in vcfs.Sample.to_list():
		return '.'
	if mother_code == '0' and father_code == '0':
		return '.'

	if 'NOSAMPLE' in mother_code:
		mother_code='0'
	if 'NOSAMPLE' in father_code:
		father_code='0'

	# Get carrier status of parents
	if mother_code == '0':
		mother_return = '0'
	elif ped[ped.Sample==mother_code]['Carrier_Status'].to_string(index = False, header = False).strip()=='unknown':
		mother_return = 'M'
	elif int(ped[ped.Sample==mother_code]['Carrier_Status'])==1:
		mother_return = 'MC'
	elif int(ped[ped.Sample==mother_code]['Carrier_Status'])==0:
		mother_return = 'MNC'
	else:
		mother_return = 'M'

	if father_code == '0':
		father_return = '0'
	elif ped[ped.Sample==father_code]['Carrier_Status'].to_string(index = False, header = False).strip()=='unknown':
		father_return = 'F'
	elif int(ped[ped.Sample==father_code]['Carrier_Status'])==1:
		father_return = 'FC'
	elif int(ped[ped.Sample==father_code]['Carrier_Status'])==0:
		father_return = 'FNC'
	else:
		father_return = 'F'

	# Check if parents have the same variant
	vid = row['variant_id']
	other_samps = var_file[var_file.variant_id == vid]['Sample'].to_list()

	if mother_code in other_samps:
		mother = True
	if father_code in other_samps:
		father = True

	if mother and father:
		return 'both'
	elif mother:
		return mother_return
	elif father:
		return father_return

	# If neither parent has variant, see if variant was present in original VCF
	chr = row['Chrom']
	pos = row['Pos']
	alt = row['Alt']

	if mother_code != '0' and mother_code in vcfs.Sample.to_list():
		mother_vcf = vcfs[vcfs.Sample == mother_code]['filename'].to_string(index = False, header = False).strip()

		ad = find_in_vcf(mother_vcf, chr, pos, alt)
		if ad:
			# To be called as de novo, read depth in parents needs to be at least 15
			# Additionally, the number of alt reads in the parents should be 0
			# To be called as inherited, at least 20% of the reads in the parent need to be the alternate allele
			alt_read = int(ad[1])
			ref_read = int(ad[0])

			if alt_read < 0 or ref_read < 0:
				print(mother_vcf, chr, pos, alt)
				print(error)
				unknown = True
			if alt_read==0 and ref_read==0:
				unknown = True
			elif float(alt_read)/(alt_read+ref_read) >= 0.2:
				mother = True
			elif alt_read > 0:
				unknown = True
			if sum([ref_read, alt_read]) < 15:
				unknown = True

	if father_code != '0' and father_code in vcfs.Sample.to_list():
		father_vcf = vcfs[vcfs.Sample == father_code]['filename'].to_string(index = False, header = False).strip()

		ad = find_in_vcf(father_vcf, chr, pos, alt)
		if ad:
			alt_read = int(ad[1])
			ref_read = int(ad[0])

			if alt_read < 0 or ref_read < 0:
				print(father_vcf, chr, pos, alt)
				print(error)
				unknown = True
			if alt_read==0 and ref_read==0:
				unknown = True
			elif float(alt_read)/(alt_read+ref_read) >= 0.2:
				father = True
			elif alt_read > 0:
				unknown = True
			if sum([ref_read, alt_read]) < 15:
				unknown = True


	if mother and father:
		return 'both'
	elif mother:
		return mother_return
	elif father:
		return father_return

	# To be called de novo, make sure we have WGS available for both parents
	if mother_code=='0' or father_code=='0':
		unknown = True

	if mother_code not in vcfs.Sample.to_list() or father_code not in vcfs.Sample.to_list():
		unknown = True

	if unknown:
		return '.'
	else:
		return 'de novo'

var_file['inheritance'] = var_file.apply(lambda row: add_inherit(row), axis = 1)

var_file.to_csv(output_table, sep = '\t', index = False)
