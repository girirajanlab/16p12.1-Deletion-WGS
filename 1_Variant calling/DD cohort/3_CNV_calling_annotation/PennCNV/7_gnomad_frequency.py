import pandas as pd
import numpy as np

# Annotate calls with their frequency in gnomAD SV v2.1 based on 50% reciprocal overlap

# Input and output files
INPUT="/path/to/input/file.csv" # Use the output from script 6_intracohort_freq.py
OUTPUT="/paht/to/output/file.csv"

# Reference files
gnomad_sv="/path/to/gnomad_sv_v2.1/sv/sites.bed.gz" # This file can be downloaded here: https://gnomad.broadinstitute.org/downloads#v2-structural-variants

# Load in files
svs=pd.read_csv(INPUT)
gnomad=pd.read_csv(gnomad_sv, sep='\t')

# Filter gnomAD SV calls
# create gnomad column for chromosome with chr
gnomad['chrom'] = gnomad['#chrom'].apply(lambda s: 'chr{}'.format(s))

# filter for only dups and dels
gnomad = gnomad[gnomad.svtype.isin(['DEL', 'DUP'])].copy()

# filter for pass only
gnomad = gnomad[gnomad.FILTER == 'PASS'].copy()

# set gnomad frequency to type float
gnomad['AF'] = gnomad.AF.astype(float)

# Loop over every chromosome to identify CNVs with 50% reciprocal overlaps in gnomAD SV
chromosomes = list(svs.Chromosome.unique())
svtypes = ['dup', 'del']
svs['gnomad_freq_overall'] = '.'

for chrom in chromosomes:
    print(chrom)
    for svtype in svtypes:
        # get all structural variants in chromosome and svtype
        sub_svs = svs[(svs.Chromosome == chrom) & (svs.Type == svtype)].copy()
        sub_gnomad = gnomad[(gnomad['chrom'] == chrom) & (gnomad.svtype == svtype.upper())].copy()
        
        # reset index in sub_gnomad
        sub_gnomad = sub_gnomad.reset_index(drop=True)
        
        # create numpy structures for end, start, and svlength
        # (numpy is faster than pandas)
        chrom_ends = sub_gnomad['end'].to_numpy()
        chrom_starts = sub_gnomad['start'].to_numpy()
        chrom_svlengths = chrom_ends - chrom_starts

        # iterate over all structural variants in cohort
        # and compare against other structural variants gnomad
        for i, row in sub_svs.iterrows():
            start = row['Start']
            end = row['End']
            length = row['Length']

            min_end   = np.minimum(chrom_ends, end)
            max_start = np.maximum(chrom_starts, start)
            max_length = np.maximum(chrom_svlengths, length)

            overlap = (min_end - max_start) > .5 * max_length

            # in some cases the structural variant overlaps with more than one gnomad SV
            # in those cases get the sum of all SV that it overlaps with
            allele_frequency = sub_gnomad.loc[overlap]['AF'].sum()
            
            # save to table
            svs.at[i, 'gnomad_freq_overall'] = allele_frequency

# Save
svs.to_csv(OUTPUT, index=False)
