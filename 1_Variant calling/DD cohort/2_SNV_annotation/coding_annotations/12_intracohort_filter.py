import pandas as pd

# Apply an intracohort frequency filter
# Variants should be present in <= 10 individuals in the cohort

# Input and output files
input_file='/path/to/input/file.csv' # Use the output of script 11_filter_variant_types.py here
output_file='/path/to/output/table.csv'

df = pd.read_csv(input_file)

# Reorder columns
columns = ['Sample', 'Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Mut_type', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
		'Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19', 'GeneDetail.wgEncodeGencodeBasicV19', 'ExonicFunc.wgEncodeGencodeBasicV19', 'AAChange.wgEncodeGencodeBasicV19',
		'ge_AF', 'ge_AF_nfe', 'ge_AF_asj', 'ge_AF_eas', 'ge_AF_amr', 'ge_AF_afr',
		'gg_AF', 'gg_AF_nfe', 'gg_AF_asj', 'gg_AF_eas', 'gg_AF_amr', 'gg_AF_afr',
		'CADD_PHRED', 'CADD_RawScore', 'ClinVar_CLNDN', 'ClinVar_CLNDISDB', 'ClinVar_CLNREVSTAT', 'ClinVar_CLNSIG', 'ClinVar_ALLELEID',
		'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB']
df = df[columns]

# Filter for intracohort frequency
df['variant_id'] = df['Chrom'] + '_' + df['Pos'].astype(str) + '_' + df['Ref'] + '_' + df['Alt']
variant_counts = df['variant_id'].value_counts()
df['cohort_count'] = df['variant_id'].map(variant_counts)

df = df[df['cohort_count'] <= 10]

df.to_csv(output_file, index=False)
