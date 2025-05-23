import numpy as np

# Identify samples with STR lengths greater than 2 SD of the cohort mean

# Input and output files
INPUT_TABLE="/path/to/input/table.tsv" # Use the table generated from script 2_vcf2table.py
OUTPUT_TABLE="/path/to/output/table.tsv"

# Load files
fin=open(INPUT_TABLE, 'r')
fout=open(OUTPUT_TABLE, 'w')

# Iterate through input table and identify STR expansions
for line in fin:
	# header
	if line.startswith('chrom'):
		# get sample info
		samples = line.strip().split('\t')[2:]
		# write out header
		outline = 'chrom\tpos\tsample\tzscore\tlongest_allele\t'
		outline = outline + 'both_alleles\tcohort_mean\t'
		outline = outline + 'cohort_mode\tcohort_std\n'
		fout.write(outline)
		continue

	sline = line.strip().split('\t')
	chrom = sline[0]
	pos   = sline[1]
	reps  = sline[2:]

	reps_all = []
	num_missing = 0
	for alleles in reps:
		if alleles == '.':
			num_missing = num_missing + 1
			continue
		alleles = alleles.split(',')
		for allele in alleles:
			reps_all.append(int(allele))

	# if more than 20% of cohort is missing calls skip
	if num_missing > (len(samples) * .2):
		continue

	# if all alleles are the same go to next one
	# there will be no outliers here
	if len(list(set(reps_all))) <= 1:
		continue

	# get mean and standard deviation
	mean = np.mean(reps_all)
	sd   = np.std(reps_all)
	mode = max(set(reps_all), key=reps_all.count)

	# go through each sample again and if > 2sd write out
	for i in range(len(samples)):
		sample  = samples[i]
		alleles = reps[i]
		if alleles == '.':
			continue
		alleles = alleles.split(',')
		alleles = [int(s) for s in alleles]
		longest_allele = max(alleles)
		zscore = (longest_allele - mean) / sd
		if zscore > 2:
			pass
			soutline = [chrom, pos, sample, zscore, 
						longest_allele, reps[i], mean,
						mode, sd]
			soutline = [str(s) for s in soutline]
			outline = '\t'.join(soutline)
			outline = outline + '\n'
			fout.write(outline)
fout.close()
fin.close()
