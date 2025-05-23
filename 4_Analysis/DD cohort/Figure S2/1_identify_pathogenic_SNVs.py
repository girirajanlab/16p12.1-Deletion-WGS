import pandas as pd
import numpy as np

# Identify pathogenic SNVs

# Definitions of pathogenicity:
# 1. We will only consider LOF variants for the gene lists (missense is still ok for ClinVar)
# 2. We will filter ClinVar for the criteria provided, followed by manual curation for relevant developmental phenotypes and removal of autosomal recessive variants

# Input and output files
SNVS="/path/to/SNV/variants.csv" # Use the output of script 1_Variant calling/DD cohort/2_SNV_annotation/coding_annotations/14_loeuf_scores.py
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
GENE_ANNO="/path/to/gene/annotations.csv" # Use the output of script 2_Analysis preparation\Gene_Annotations\5_add_loeuf.py

OUTPUT_TAB="/path/to/preliminary/pathogenic/variants.csv" # Preliminary pathogenic ClinVar variants will be manually filtered for relevant developmental phenotypes and autosomal recessive variants will be removed

# Load data
snvs=pd.read_csv(SNVS)

df=pd.read_csv(TABS1A)
df=df[df['Estonian Biobank Sample']!='X']
df=df[df.Relationship=='Proband']

# Remove any unneeded participants (for example, Estonian samples, parents, siblings)
snvs=snvs[snvs.Sample.isin(df.Sample.to_list())]

# Pathogenic SNVS are:
patho_snvs = snvs.copy()
merge_cols = patho_snvs.columns.to_list()
patho_snvs['CLINVAR_Pathogenic'] = '.'
# Any SNV labelled as "Pathogenic" or "likely pathogenic" in ClinVar
# AND criteria is "criteria_provided,_multiple_submitters,_no_conflicts" or "reviewed_by_expert_panel"
patho_snvs.loc[(patho_snvs.ClinVar_CLNSIG.isin(['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'])) & (patho_snvs.ClinVar_CLNREVSTAT.isin(['criteria_provided,_multiple_submitters,_no_conflicts', 'reviewed_by_expert_panel'])),
 'CLINVAR_Pathogenic'] = 'X'

# In a Tier S SFARI gene
gene_anno = pd.read_csv(GENE_ANNO)
print(gene_anno.SFARI_gene_score.value_counts())
sfari_genes = gene_anno[gene_anno.SFARI_gene_score=='S'][['gene_id', 'SFARI_gene_score']]
sfari_snv = pd.merge(snvs, sfari_genes, left_on = 'Gene_id_', right_on = 'gene_id', how = 'inner')
sfari_snv['SFARI_Tier_S'] = 'X'
patho_snvs = pd.merge(patho_snvs, sfari_snv, on = merge_cols, how = 'left')

# In Tier 1 genes from the Geisinger DBD database
print(gene_anno.Geisinger_DBD_Tier.value_counts())
dbd_genes = gene_anno[gene_anno.Geisinger_DBD_Tier.isin(['1', 1])][['gene_id', 'Geisinger_DBD_Tier']]
dbd_snv = pd.merge(snvs, dbd_genes, left_on = 'Gene_id_', right_on = 'gene_id', how = 'inner')
dbd_snv['DBD_Tier_1'] = 'X'
patho_snvs = pd.merge(patho_snvs, dbd_snv, on = merge_cols, how = 'left')

# Fill NA with .
patho_snvs.fillna('.', inplace = True)
# Remove extra columns
print(patho_snvs.columns)
patho_snvs = patho_snvs[['Sample', 'Chrom', 'Pos', 'Ref', 'Alt', 'Mut_type',
	   'ClinVar_CLNDN', 'ClinVar_CLNDISDB', 'ClinVar_CLNREVSTAT', 'ClinVar_CLNSIG', 'ClinVar_ALLELEID',
       'Gene_symbol', 'Gene_id_', 'LOEUF', 'CLINVAR_Pathogenic',
       'SFARI_gene_score', 'SFARI_Tier_S', 'Geisinger_DBD_Tier', 'DBD_Tier_1']]

# Add final annotation
patho_snvs['Pathogenic'] = (patho_snvs[['CLINVAR_Pathogenic', 'SFARI_Tier_S', 'DBD_Tier_1']]=='X').any(axis=1)
print(patho_snvs.Pathogenic.value_counts())

# If variant is not LOF and not in ClinVar, remove pathogenic annotation
# Now that ClinVar are annotated, remove any non-ClinVar variant that is LOF
patho_snvs.loc[~((patho_snvs.CLINVAR_Pathogenic=='X') | (patho_snvs.Mut_type=='lof')), ['SFARI_Tier_S', 'DBD_Tier_1', 'Pathogenic']] = False

print(patho_snvs.Pathogenic.value_counts())

# Manually curate SNVs marked pathogenic through OMIM OR ClinVar annotations to check for neurodevelopmental phenotypes
patho_snvs=patho_snvs[patho_snvs.Pathogenic]
patho_snvs['CLINVAR_Manual_review']=np.nan
patho_snvs.loc[patho_snvs['CLINVAR_Pathogenic']=='X', 'CLINVAR_Manual_review']='X'

patho_snvs.drop_duplicates(inplace=True)
patho_snvs.sort_values(by=['Sample', 'Chrom'], inplace=True)

patho_snvs.to_csv(OUTPUT_TAB, index = False)