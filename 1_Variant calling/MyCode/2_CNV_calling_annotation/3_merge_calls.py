import pandas as pd

# Merge calls within a sample if
# 1. 2 calls overlap
# 2. The gap is < 50kb and <20% of the combined CNV length

# Input and output files
QC_file="/path/to/qc/file.txt" # Use the output of script 2_check_QC.py
SNP_Pos="/path/to/SNP/Position/file.txt"
pcnv_output="/path/to/PennCNV/output/directory" # Use the OUTPUT_DIR from script 1_run_PennCNV.sh
output_file="/path/to/output/file.txt"

# Load files
filter = pd.read_csv(QC_file, sep = '\t')
passed = filter[(filter.Pass_auto == True) & (filter.Pass_chrX==True)]['Sample'].to_list()

snp_df = pd.read_csv(SNP_Pos, sep = '\t')
snp_df=snp_df[(~snp_df.Chr.isnull()) & (~snp_df.Pos.isnull())]
snp_df=snp_df[snp_df.Chr!='0']
snp_df.sort_values(by=['Chr', 'Pos'], inplace=True)
snp_df.reset_index(inplace=True, drop=True)

columns=['Location', 'NumSNP', 'Length', 'State', '`echo', 'Name\t.Log R Ratio\t.B Allele Freq', ';', '(gunzip', '-c', 'File',  '|', 'tail', '-n',  '+2)`', 'StartSNP', 'EndSNP']

calls=pd.DataFrame()
for BATCH in ['OMNI', 'GSA']:
	auto=pd.read_csv(f'{pcnv_output}/QC/goodcnv/autosome_{BATCH}.goodcnv', delim_whitespace=True, header=None, names = columns)
	chrX=pd.read_csv(f'{pcnv_output}/QC/chrX/autosome_{BATCH}.goodcnv', delim_whitespace=True, header=None, names = columns)

	calls=pd.concat([calls, auto])
	calls=pd.concat([calls, chrX])
calls=calls[['Location', 'NumSNP', 'Length', 'State', 'File', 'StartSNP', 'EndSNP']]

# Helper functions
def to_merge(start, end, start2, end2, threshold = 0.2):
	gap = max(start, int(start2)) - min(end, int(end2))
	# Overlap
	if gap <= 0:
		return True
	# Is gap < 50kb?
	if gap > 50000:
		return False
	# Is gap < 20% of the combined length?
	total_len = max(end, int(end2)) - min(start, int(start2))
	if gap <= (threshold*total_len):
		return True
	return False

def test_merge(row):
	chr = row['Chr']
	start = int(row['Start'])
	end = int(row['End'])
	type = row['Type']
	length = int(row['Length'])
	sample = row['Sample']
	id = row.name

	# Find other calls that have a 50% reciprocal overlap with the current call
	same_chrom = calls[(calls.Chr==chr) & (calls.Type==type) & (calls.Sample==sample)]
	# Remove current call
	not_curr = same_chrom[same_chrom.index != id]
	# If there are no other calls that meet these criteria
	if not_curr.shape[0]==0:
		return False
	# Find calls with overlap or small gaps
	merge = not_curr[not_curr.apply(lambda row: to_merge(start, end, row['Start'], row['End']), axis = 1)]
	if merge.shape[0]==0:
		return False

	return merge.index.to_list()

def combine_ids(row):
	in_ids = row['direct_merge']
	out_ids = in_ids

	# Iteratively add in extended overlaps until the list stops changing
	while True:
		for id in in_ids:
			add_ids = calls[calls.index==id]['direct_merge']
			for a_id in add_ids:
				for b_id in a_id:
					if b_id not in out_ids:
						out_ids.append(b_id)
		if in_ids==out_ids:
			break

		in_ids = out_ids
	out_ids.sort()
	return(out_ids)

def merge_indirect(indirect):
	to_merge = calls[calls.index.isin(indirect)]

	# Get information for the final call
	chrom = ' '.join(list(to_merge.Chr.unique())) # Should only give one chromosome
	start = min(to_merge.Start)
	end = max(to_merge.End)
	type = ' '.join(list(to_merge.Type.unique())) # Should also only be one type
	zygosity = ' '.join(list(to_merge.Zygosity.unique())) # May be multiple values
	length = int(end) - int(start)
	sample = ' '.join(list(to_merge.Sample.unique())) # Should only be one sample

	# To get StartSNP, EndSNP, and NumSNP, consult SNP lists
	startsnp_idx = min(snp_df[snp_df.Name.isin(to_merge.StartSNP.to_list())].index)
	startsnp = snp_df.iloc[startsnp_idx]['Name']
	endsnp_idx = max(snp_df[snp_df.Name.isin(to_merge.EndSNP.to_list())].index)
	endsnp = snp_df.iloc[endsnp_idx]['Name']
	numsnp = endsnp_idx - startsnp_idx

	# Return merged call as series
	s = pd.Series([chrom, start, end, type, zygosity, length, numsnp, sample, startsnp, endsnp],
			index = ['Chr', 'Start', 'End', 'Type', 'Zygosity', 'Length', 'NumSNP', 'Sample', 'StartSNP', 'EndSNP'])
	return s

# Remove samples that failed QC
calls = calls[calls.File.isin(passed)]

# Add new fields for easy merging and reformatting
calls['Chr'] = calls.Location.str.split(':', expand = True)[0]
calls['Start'] = calls.Location.str.split(':', expand = True)[1].str.split('-', expand = True)[0]
calls['End'] = calls.Location.str.split(':', expand = True)[1].str.split('-', expand = True)[1]
calls['Length'] = calls.Length.str.split('=', expand = True)[1].str.replace(',', '')
calls['StartSNP'] = calls.StartSNP.str.split('=', expand = True)[1]
calls['EndSNP'] = calls.EndSNP.str.split('=', expand = True)[1]
calls['NumSNP'] = calls.NumSNP.str.split('=', expand = True)[1].str.replace(',', '')

# Translations of the CN and state
# NOTE: "state3,cn=2" is only for the X chromosome, it indicates a chrX duplication in a male
# CN=2 would not be reported in an autosomal region
type_dict = {'state1,cn=0': 'DEL', 'state2,cn=1':'DEL', 'state5,cn=3':'DUP', 'state6,cn=4':'DUP', 'state3,cn=2':'DUP'}
calls['Type'] = calls.State.map(type_dict)
zygosity_dict = {'state1,cn=0':'homozygous', 'state2,cn=1':'heterozygous', 'state5,cn=3':'heterozygous', 'state6,cn=4':'homozygous', 'state3,cn=2':'hemizygous'}
calls['Zygosity'] = calls.State.map(zygosity_dict)
# Convert relevant columns to numeric
for col in ['Start', 'End', 'Length', 'NumSNP']:
	calls[col]=calls[col].astype(int)

# Merge calls, as needed
calls['direct_merge'] = calls.apply(test_merge, axis = 1)

# Separate out calls that do not merge
no_merge = calls[calls.direct_merge==False]
# Convert to final format
no_merge_out = no_merge[['Chr', 'Start', 'End', 'Type', 'Zygosity', 'Length', 'NumSNP', 'Sample', 'StartSNP', 'EndSNP']]
no_merge_out['Merge'] = False

# Merge calls that do need to be merged
calls = calls[calls.direct_merge != False]
print(calls)
# Make a tiling path for merging
calls['indirect_merge'] = calls.apply(combine_ids, axis = 1)
# Write the merge of the indirect overlap lists
# Only perform merge on unique sets of IDs
unique_sets = []
out_calls = []
for i in calls['indirect_merge']:
        if i not in unique_sets:
                unique_sets.append(i)
                out_calls.append(merge_indirect(i))
merge_out = pd.DataFrame(out_calls)
merge_out['Merge'] = True

# Write all calls to file
final_out = pd.concat([no_merge_out, merge_out], axis = 0)
final_out.to_csv(output_file, sep = '\t', index = False)

