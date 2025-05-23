import pandas as pd
import numpy as np

# Calculate the intracohort frequency based on 50% reciprocal overlap of CNVs

# Input and output files
INPUT="/path/to/input/file.txt" # Use the calls from script 5_anno_segdup_centel.sh here
OUTPUT="/path/to/output/file.csv"

# Load in input calls
svs=pd.read_csv(INPUT, sep='\t')

# Loop over every CNV to identify calls with a 50% reciprocal overlap
chromosomes = list(svs.Chromosome.unique())
svtypes = ['dup', 'del']

for chrom in chromosomes:
    for svtype in svtypes:
        # get all structural variants in chromosome of the svtype
        sub_svs = svs[(svs.Chromosome == chrom) & (svs.Type == svtype)].copy()

        # create numpy structures for end, start, and svlength
        chrom_ends = sub_svs['End'].to_numpy()
        chrom_starts = sub_svs['Start'].to_numpy()
        chrom_svlengths = sub_svs['Length'].to_numpy()

        # iterate over all structural variants in chrom
        # and compare against other structural variants in chrom
        for i, row in sub_svs.iterrows():
            start = row['Start']
            end = row['End']
            length = row['Length']

            min_end   = np.minimum(chrom_ends, end)
            max_start = np.maximum(chrom_starts, start)
            max_length = np.maximum(chrom_svlengths, length)

            overlap = (min_end - max_start) > .5 * max_length

            intra_cohort_count = overlap.sum()
            svs.at[i, 'intracohort_count'] = intra_cohort_count
            
svs['intracohort_count'] = svs['intracohort_count'].astype(int)

# Save
svs.to_csv(OUTPUT, index=False)