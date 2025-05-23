import pandas as pd

from upsetplot import plot
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42

# Update the pathogenic SNVs by incorporating the manual review

# Input and output files
MANUAL_REVIEW="/path/to/manually/reviewed/SNVs.csv" # Preliminary pathogenic SNVs were manually reviewed to remove SNVs with autosomal recessive MOI and those with unrelated phenotypes
TABS1A="/path/to/Table_S1A.csv" # Use the output of script 3_Data preparation/DD_cohort/1_make_table_s1a.py
CNVS="/path/to/CNV/calls.txt" # Use the output from script 1_Variant calling/DD cohort/3_CNV_calling_annotation/merge_all_cnvs/2_annotate_loeuf.py

PATHOGENIC_SNVS="/output/final/pathogenic/SNVs.csv"
UPSET_FIG="/path/to/upset/figure.pdf" # Figure presented in Fig S2A

# Load data
df = pd.read_csv(MANUAL_REVIEW)

df.loc[df.CLINVAR_AD=='N', 'CLINVAR_Manual_review']='N'

df['Pathogenic'] = '.'
df.loc[((df[['SFARI_Tier_S', 'DBD_Tier_1']]!='.').any(axis=1)) | (df.CLINVAR_Manual_review=='Y'), 'Pathogenic'] = 'X'

print(df.Pathogenic.value_counts())
df=df[df.Pathogenic=='X']

df.loc[df.CLINVAR_Manual_review=='N', 'CLINVAR_Pathogenic']='.'
df.loc[df.CLINVAR_AD=='N', 'CLINVAR_Pathogenic']='.'

# Reorder columns
df = df[['Sample', 'Chrom', 'Pos', 'Ref', 'Alt', 'Mut_type',
       'ClinVar_CLNDN', 'ClinVar_CLNDISDB', 'ClinVar_CLNREVSTAT', 'ClinVar_CLNSIG', 'ClinVar_ALLELEID',
       'Gene_symbol', 'Gene_id_', 'LOEUF', 'CLINVAR_Pathogenic',
       'SFARI_gene_score', 'SFARI_Tier_S', 'Geisinger_DBD_Tier', 'DBD_Tier_1', 'Pathogenic']]

# Save to file
df.to_csv(PATHOGENIC_SNVS, index = False)

# Make an upset plot show casing the types of mutations 16p12.1 deletion probands have
# Create an upset plot showing the probands with mutations in each category of Pathogenic variant
# Categories are:
# 1. LOF mutations in Tier S or 1 SFARI genes
# 2. LOF mutations in Tier 1 or 2 DBD genes
# 3. "Pathogenic" or "Likely pathogenic" variant in ClinVar with relevant phenotypes/mode of action (from a manual screen)

# Get all probands
cohort_info=pd.read_csv(TABS1A)
cohort_info=cohort_info[cohort_info['Estonian Biobank Sample']!='X']
cohort_info=cohort_info[cohort_info.Relationship=='Proband']
cohort_info=cohort_info[cohort_info.WGS=='X']

df.replace('.', 0, inplace=True)
df.replace('X', 1, inplace=True)

# Add in Pathogenic CNV data
cnv=pd.read_csv(CNVS, sep='\t')
cnv=cnv[cnv.NEJM!='.'][['Sample', 'NEJM']]
cnv['Pathogenic_CNV']=1

df=pd.merge(df, cnv[['Sample', 'Pathogenic_CNV']], on='Sample', how='outer')
df.fillna(0, inplace=True)

# Condense data on the proband level
proband_df=df[['Sample', 'CLINVAR_Pathogenic', 'SFARI_Tier_S', 'DBD_Tier_1', 'Pathogenic_CNV']].groupby('Sample').sum()
for col in ['CLINVAR_Pathogenic', 'SFARI_Tier_S', 'DBD_Tier_1', 'Pathogenic_CNV']:
    proband_df.loc[proband_df[col]>0, col]=True
    proband_df.loc[proband_df[col]==0, col]=False
proband_df['None']=True
proband_df.loc[(proband_df.CLINVAR_Pathogenic) | (proband_df.SFARI_Tier_S) | (proband_df.DBD_Tier_1) | (proband_df.Pathogenic_CNV), 'None']=False

# Update column names for plot
proband_df=proband_df[['CLINVAR_Pathogenic', 'SFARI_Tier_S', 'DBD_Tier_1', 'Pathogenic_CNV', 'None']]
proband_df.columns=['ClinVar', 'Tier S SFARI', 'Tier 1 DBD', 'Pathogenic CNV', 'None']

# Make upset plot
plot_df=proband_df.groupby(['ClinVar', 'Tier S SFARI', 'Tier 1 DBD', 'Pathogenic CNV', 'None']).size()
plot(plot_df)
plt.savefig(UPSET_FIG)
