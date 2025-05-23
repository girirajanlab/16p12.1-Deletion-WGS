import pandas as pd

# Filter variants to include
# 1. Exonic and splicing variants
# 2. Missense and splice variants with CADD >= 25 or LOF variants

# Input and output files
input_table='/path/to/input/table.tsv' # Use the all-sample output table of script 4_filter_pop_frequency.sh
output_file='/path/to/output/table.csv'

# Read table
df = pd.read_csv(input_table, sep = '\t', header=None,
		names = ['Sample', 'Chrom', 'Pos', 'Ref', 'Alt', 'Qual',
                'Func.wgEncodeGencodeBasicV38', 'Gene.wgEncodeGencodeBasicV38', 'GeneDetail.wgEncodeGencodeBasicV38', 'ExonicFunc.wgEncodeGencodeBasicV38', 'AAChange.wgEncodeGencodeBasicV38',
                'ge_AF', 'ge_AF_nfe', 'ge_AF_asj', 'ge_AF_eas', 'ge_AF_amr', 'ge_AF_afr',
                'gg_AF', 'gg_AF_nfe', 'gg_AF_asj', 'gg_AF_eas', 'gg_AF_amr', 'gg_AF_afr',
                'CADD_PHRED', 'CADD_RawScore',
                'GT', 'DP', 'AD', 'GQ', 'PL'])

# Replace any missing CADD scores with NAN
df.loc[df.CADD_PHRED=='.', 'CADD_PHRED']=np.nan

# Fix mistakes in parsing unicode hex
df['Func.wgEncodeGencodeBasicV38'] = df['Func.wgEncodeGencodeBasicV38'].str.replace('\\x3b', ';')
df['Gene.wgEncodeGencodeBasicV38'] = df['Gene.wgEncodeGencodeBasicV38'].str.replace('\\x3b', ';')
df['ExonicFunc.wgEncodeGencodeBasicV38'] = df['ExonicFunc.wgEncodeGencodeBasicV38'].str.replace('\\x3b', ';')

# Keep only exonic and splicing
df=df[(df['Func.wgEncodeGencodeBasicV38'].str.contains('exonic')) | (df['Func.wgEncodeGencodeBasicV38'].str.contains('splicing'))]
df=df[(df['Func.wgEncodeGencodeBasicV38']!='ncRNA_exonic') & (df['Func.wgEncodeGencodeBasicV38']!='ncRNA_splicing')]

# Filter variant types
missense_variant_classes = ['nonsynonymous_SNV', 'nonframeshift_deletion', 'nonframeshift_insertion', 'nonframeshift_substitution']
lof_variant_classes = ['stopgain', 'stoploss', 'frameshift_deletion', 'frameshift_insertion', 'frameshift_substitution']

df['Mut_type']='exonic'
df.loc[df['Func.wgEncodeGencodeBasicV38'].str.contains('splicing'), 'Mut_type']='splice'
df.loc[df['ExonicFunc.wgEncodeGencodeBasicV38'].isin(missense_variant_classes), 'Mut_type']='missense'
df.loc[df['ExonicFunc.wgEncodeGencodeBasicV38'].isin(lof_variant_classes), 'Mut_type']='lof'

df=df[df['Mut_type'].isin(['splice', 'missense', 'lof'])]

# Apply a CADD score filter only to missense and splice
df.CADD_PHRED=df.CADD_PHRED.astype(float)
df=df[(df.CADD_PHRED>=25) | (df.Mut_type=='lof')]

# Apply an intracohort filter 
# Variants should be present in <= 10 individuals in the cohort
df['variant_id']=df.Chrom.astype(str)+'_'+df.Pos.astype(str)+'_'+df.Ref.astype(str)+'_'+df.Alt.astype(str)
counts=df.variant_id.value_counts().to_dict()
df['intracohort_frequency']=df.variant_id.map(counts)
df=df[df.intracohort_frequency<=10]

# in some cases a variant has Func="exonic;splicing" region or Func="ncRNA_exonic;splicing"
# in which case, only keep the gene whose function we're more interested in
# exonic > splicing > ncRNA_exonic
for i, row in df.iterrows():
	locations = row['Func.wgEncodeGencodeBasicV38'].split(';')
	genes = row['Gene.wgEncodeGencodeBasicV38'].split(';')

	if len(locations) == 1:
		continue

	if 'exonic' in locations:
		index = locations.index('exonic')
		df.at[i, 'Func.wgEncodeGencodeBasicV38'] = locations[index]
		df.at[i, 'Gene.wgEncodeGencodeBasicV38'] = genes[index]

	elif 'splicing' in locations:
		index = locations.index('splicing')
		df.at[i, 'Func.wgEncodeGencodeBasicV38'] = locations[index]
		df.at[i, 'Gene.wgEncodeGencodeBasicV38'] = genes[index]

# Save to file
df.to_csv(output_file, index=False)
