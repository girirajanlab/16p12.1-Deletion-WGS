import pandas as pd

# Annotate STR calls with locus information from reference file

# Input and output files
INPUT_TABLE="/path/to/input/table.tsv" # Use the output from script 3_2SD_STR.py
OUTPUT_TABLE="/path/to/output/table.tsv"
ANNOVAR_INPUT="/path/to/output/table/for/annovar.tsv" # This file will be used in script 5_annovar.sh

# Reference files
GANGSTR_REF="/path/to/hg19/gangstr/reference/hg38_ver16.bed.gz"  # This file can be downloaded here: https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver16.bed.gz

# Load files
bed=pd.read_csv(GANGSTR_REF, sep='\t', header=None)
bed = bed.set_index([0,1], drop=False)
bed.columns = ['chrom', 'start', 'end', 'period', 'motif']
bed['ref_length'] = bed.end - bed.start + 1

df = pd.read_csv(INPUT_TABLE, sep='\t')
df = df.set_index(['chrom', 'pos'], drop=False)

# Helper functions
def construct_allele(motif, ref_length):
    return motif* int(ref_length/len(motif))

# Annotate calls with files from reference file
bed_at_df = bed.loc[df.index]

df['end'] = bed_at_df['end']
df['motif'] = bed_at_df['motif']
df['ref_length'] = bed_at_df['ref_length']
df['period'] = bed_at_df['period']
df['alt_length'] = df['longest_allele'] * df['period']

# Identify ref and alt alleles
df['ref_allele'] = [construct_allele(s[0], s[1]) for s in zip(df['motif'], df['ref_length'])]
df['alt_allele'] = [construct_allele(s[0], s[1]) for s in zip(df['motif'], df['alt_length'])]

# Save annotations
df.to_csv(OUTPUT_TABLE, sep='\t', index=False)

# Filter for columns needed for annotation and save
cols = ['chrom', 'pos', 'end', 'ref_allele', 'alt_allele']
df = df[cols]
df = df.drop_duplicates()
df.to_csv(ANNOVAR_INPUT, sep='\t', index=False)
