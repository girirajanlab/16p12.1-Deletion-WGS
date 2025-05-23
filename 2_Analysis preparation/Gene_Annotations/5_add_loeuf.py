import pandas as pd

# Input and output files
LOEUF='/path/to/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz' # Downloaded from https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
GENE_ANNO="/paht/to/output/gene/annotations.csv"

# Add LOEUF scores
loeuf = pd.read_csv(LOEUF, sep = '\t', copmression='gzip')
loeuf = loeuf[['gene_id', 'oe_lof_upper']]
loeuf.index = loeuf.gene_id

df = pd.read_csv('intermediate_tables/4_wang_epilepsy_annotations.csv')
df['LOEUF'] = df.gene_id.map(loeuf.oe_lof_upper.to_dict())

df['Constrained_LOEUF'] = 0
df.loc[df.LOEUF <= 0.35, 'Constrained_LOEUF'] = 1

df['SFARI_gene_score']=df['gene-score']
df['Geisinger_DBD_Tier']=df.Tier

df = df[['gene_id', 'gene_symbol', 'ASD_risk_genes_TADA_FDR0.3',
       'ASD_coexpression_networks_Willsey2013', 'BrainExpressed_Kang2011',
       'Constrained_LOEUF', 'PSD_Genes2Cognition',
       'Developmental_delay_DDD', 'CHD8_targets_Cotney2015_Sugathan2014',
       'FMRP_targets_Darnell2011', 'Geisinger_DBD_Tier', 'DD_G2P',
       'SFARI_gene_score', 'SZDB_schizophrenia', 'Wang_Epilepsy', 'LOEUF']]

print(df)

for col in ['ASD_risk_genes_TADA_FDR0.3',
       'ASD_coexpression_networks_Willsey2013', 'BrainExpressed_Kang2011',
       'Constrained_LOEUF', 'PSD_Genes2Cognition',
       'Developmental_delay_DDD', 'CHD8_targets_Cotney2015_Sugathan2014',
       'FMRP_targets_Darnell2011', 'Geisinger_DBD_Tier', 'DD_G2P',
       'SFARI_gene_score', 'SZDB_schizophrenia', 'Wang_Epilepsy']:
	print(df[col].value_counts())

# Save to file
df.to_csv(GENE_ANNO, index = False)

# Also save a file of just gene names
df.drop_duplicates(subset='gene_id', keep='first', inplace=True)
df[['gene_id']].to_csv('All_genes.csv', index=False, header=False)