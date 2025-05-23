import gzip

# Convert VCFs into tables
# Run this step on all chromosome VCFs
INPUT_VCF="/path/to/input/file.vcf.gz" # Use the merged VCFs ouput from script 1_run_mergeSTR_statSTR.sh
OUTPUT_TABLE="/path/to/output/file.tsv"

# Helper functions
def get_gt_info_field(gt_info, field_i):
	if gt_info == '.':
		return '.'
	gt_info = gt_info.split(':')
	return gt_info[field_i]

# Parse VCF
fin=gzip.open(INPUT_VCF)
fout=open(OUTPUT_TABLE, 'w')

for line in fin:
	line = line.decode()
	
	# write header
	if line.startswith('#CHROM'):
		sline = line.strip().split('\t')
		samples = sline[9:]
		outline = 'chrom\tpos\t'
		outline = outline + '\t'.join(samples) + '\n'
		fout.write(outline)
		continue
	
	if line.startswith('#'):
		continue
	
	sline = line.strip().split('\t')

	chrom = sline[0]
	pos = sline[1]
	gt_fields = sline[8]
	gt_infos = sline[9:]
		
	# find which gt_field is REPCN
	repcn_i = -1 # dummy
	gt_fields = gt_fields.split(':')
	for i in range(len(gt_fields)):
		field = gt_fields[i]
		if field == 'REPCN':
			repcn_i = i
			break

	# check that repcn_i was found
	if repcn_i == -1:
		print('REPCN not found {} {}\n'.format(chrom, pos))
		continue

	# get rep counts
	repcns = [get_gt_info_field(s, repcn_i) for s in gt_infos]

	# create outline
	outline = chrom + '\t' + str(pos) + '\t'
	for repcn in repcns:
		outline = outline + repcn + '\t'
		
	# remove trailing tab and add new line
	outline = outline[:-1] + '\n'
	fout.write(outline)

fout.close()
fin.close()
