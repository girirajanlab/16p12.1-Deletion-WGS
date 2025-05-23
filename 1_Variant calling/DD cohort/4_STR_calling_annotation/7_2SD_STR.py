import pandas as pd
import gzip

# Identify samples and loci with STR lengths > 2SD than the cohort mean

# Input and output files
VCF="/path/to/input.vcf.gz" # Use the cohort VCF output from script 6_run_statSTR.sh
STATS="/paht/to/statSTR/statistics.tab" # Use the statistics output from script 6_run_statSTR.sh
OUT_TAB="/path/to/str/output_table.tsv"
OUT_EXP="/path/to/expansion/output_table.tsv"

# Reference files
GANGSTR_REF="/path/to/hg19/gangstr/reference/hg19_ver13.1.bed"  # This file can be downloaded here: https://s3.amazonaws.com/gangstr/hg19/genomewide/hg19_ver13_1.bed.gz

# Annotate BED file
bed=pd.read_csv(GANGSTR_REF, sep='\t', header=None)
bed['variant_id'] = bed[0] + '_' + bed[1].astype(str)
bed = bed.set_index('variant_id')

# Parse VCF
fin = gzip.open(VCF)
fout=open(OUT_TAB, 'w')
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
	ref_allele = sline[3]
	alt_alleles = sline[4]
	gt_infos = sline[9:]
	
	# get period of str
	variant_id = chrom + '_' + str(pos)
	period = bed.loc[variant_id, 3]
	
	# get lengths of all of the alleles
	if alt_alleles == '.':
		alleles = [ref_allele]
	else:
		alleles = [ref_allele] + alt_alleles.split(',')
	alleles_lengths = [len(s)/period for s in alleles]
	
	# strip extra info from gts
	gt_infos = [s.split(':')[0] for s in gt_infos]
	
	# process gts and create outline
	outline = chrom + '\t' + str(pos) + '\t'
	for gt_info in gt_infos:
		if gt_info == '.':
			outline = outline + '.\t'
			continue
		gt = gt_info.split('/')
		lengths = [alleles_lengths[int(s)] for s in gt]
		lengths = [str(int(s)) for s in lengths]
		outline = outline + ','.join(lengths) + '\t'
	
	# remove trailing tab and add new line
	outline = outline[:-1] + '\n'
	fout.write(outline)
fout.close()
fin.close()

# Identify STR expansions
vcf=pd.read_csv(OUT_TAB, sep='\t', comment='#')
vcf = vcf.set_index(['chrom', 'pos'])

stats=pd.read_csv(STATS, sep='\t')
stats = stats.set_index(['chrom', 'start'], drop=False)

# Parse reference file
bed = bed.set_index([0,1], drop=False)
bed.columns = ['chrom', 'start', 'end', 'period', 'motif']
bed['ref_length'] = bed.end - bed.start + 1
def construct_ref_allele(motif, ref_length):
    return motif* int(ref_length/len(motif))
bed['ref_allele'] = [construct_ref_allele(s[0], s[1]) for s in zip(bed['motif'], bed['ref_length'])]

fexp=open(OUT_EXP, 'w')

# Write header
outline = '{}\t{}\t{}\t{}\t{}\t'.format('chrom', 'pos', 'end', 'ref_allele', 'alt_allele')
outline = outline + '{}\t{}\t{}\t{}\t{}\t{}\n'.format('sample', 'zscore', 'longest_allele', 'cohort_mode', 'motif', 'motif_period')
fexp.write(outline)

# Iterate through calls to identify expansions
samples = list(vcf.columns)
for i, vcf_row in vcf.iterrows():
    chrom = i[0]
    pos   = i[1]

    # get mean and std dev at locus
    loc_alleles = []
    vcf_row = vcf_row.to_dict()
    for samp in samples:
        alleles = vcf_row[samp]
        if alleles == '.':
            continue
        alleles = alleles.split(',')
        for a in alleles:
            a = int(a)
            loc_alleles.append(a)
            
    standard_dev = np.std(loc_alleles)
    mean = np.mean(loc_alleles)

    # get other stats to write out
    mode = int(stats.loc[i, 'mode'])
    end = int(stats.loc[i, 'end']) - 1
    ref_allele = bed.loc[i, 'ref_allele']
    period = bed.loc[i, 'period']
    motif = bed.loc[i, 'motif']
    
    for samp in samples:
        alleles = vcf_row[samp]
        if alleles == '.':
            continue
        alleles = alleles.split(',')
        alleles = [int(s) for s in alleles]
        longest_allele = max(alleles)
        if (longest_allele - mean) > (2 * standard_dev):
            zscore = (longest_allele - mean) / standard_dev
            alt_allele = construct_ref_allele(motif, longest_allele * period)
            outline = '{}\t{}\t{}\t{}\t{}\t'.format(chrom, pos, end, ref_allele, alt_allele)
            outline = outline + '{}\t{}\t{}\t{}\t{}\t{}\n'.format(samp, zscore, longest_allele, mode, motif, period)
            fexp.write(outline)

fexp.close()

