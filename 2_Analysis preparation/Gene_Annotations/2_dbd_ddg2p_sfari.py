import pandas as pd

df = pd.read_csv('intermediate_tables/1_werling_annotations.csv')

# Geisinger DBD Database
# See Data_Files README for more information
# LOF Variant table download includes "Tier"s for genes
# Tier meanings:
# Tier 1: Genes with three or more de novo pathogenic loss-of-function variants (High confidence)
# Tier 2: Genes with two de novo pathogenic loss-of-function variant (High confidence)
# Tier 3: Genes with one de novo pathogenic loss-of-function variant (Emerging candidate)
# Tier 4: Genes with no de novo pathogenic loss-of-function variant (Emerging candidate)
# Tier AR: Genes with autosomal recessive inheritance (High confidence)
dbd = pd.read_csv('Data_Files/Full-LoF-Table-Data.csv')

df = pd.merge(df, dbd[['Gene', 'Tier']], left_on = 'gene_symbol', right_on='Gene', how = 'left')
print(df.Tier.value_counts())
# Save memory
del dbd

# Gene2Phenotype Developmental Delay set
# See Data_Files README for more information
ddg2p = pd.read_csv('Data_Files/DDG2P.csv.gz')
# We only care if a gene is present in this list
df['DD_G2P']=0
df.loc[df.gene_symbol.isin(ddg2p['gene symbol'].to_list()), 'DD_G2P']=1
print(df.DD_G2P.value_counts())
# Save memory
del ddg2p

# SFARI Gene database
# See Data_Files README for more information
sfari = pd.read_csv('Data_Files/SFARI-Gene_genes_01-11-2022release_03-14-2022export.csv')
df = pd.merge(df, sfari[['ensembl-id', 'gene-score']], left_on = 'gene_id', right_on = 'ensembl-id', how = 'left')
print(df['gene-score'].value_counts())
# Save memory
del sfari

# Save these annotations to a file
df.to_csv('intermediate_tables/2_dbd_ddg2p_sfari_annotations.csv', index = False)