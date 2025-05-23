#!bin/python
import pandas as pd

# SZDB is a database of schizophrenia genes from various sources
# We are using genes identified from GWAS, CNVs, or Exome sequencing

# GWAS
gwas_anno = pd.DataFrame()
# GWAS genes come from 2 sources - compile
clozuk = pd.read_csv('Data_Files/SZDB/GWAS-Genes/CLOZUK_GENE.txt', sep = '\t')
pgc2 = pd.read_csv('Data_Files/SZDB/GWAS-Genes/PGC2_GENE.txt', sep = '\t')
gwas_anno['gene_symbol'] = clozuk.Gene.to_list() + pgc2.Gene.to_list()
gwas_anno.drop_duplicates(keep = 'first', inplace = True)
gwas_anno['GWAS'] = 1
print(gwas_anno)

# CNV
cnv_anno = pd.DataFrame()
# There is something wrong with the header of this file - so skip it
cnv = pd.read_csv('Data_Files/SZDB/cnvgene.txt', sep = '\t', comment = '#', skiprows = [1], header = None,
    names = ['gene', 'genechr', 'genestart', 'geneend', 'gene_type', 'cytoband', 'cnvchr', 'cnvstart', 'cnvend', 'Putative_mechanism', 'CNV_test', 'Direction', 'case_number', 'control_number', 'Pregional', 'OR'])
cnv_anno['gene_symbol'] = cnv.gene
cnv_anno['CNV'] = 1
print(cnv_anno)

df = pd.merge(gwas_anno, cnv_anno, on = 'gene_symbol', how = 'outer')
print(df)

# Exome
exome_anno = pd.DataFrame()
exome = pd.read_csv('Data_Files/SZDB/Exome.txt', sep = '\t')
# Replace NA with 0
exome.fillna(0, inplace = True)
# Check deleteriousness annotations
print(exome["PoluPhen-2"].value_counts())
print(exome.SIFT.value_counts())
# Exome calls include variants which may not be deleterious
# Remove variants which are considered benign or unknown by both PolyPhen2 and SIFT
exome = exome[~((exome["PoluPhen-2"].isin(['Benign', 'benign', '-', 'unknown', '.'])) & (exome.SIFT.isin(['-', 'Tolerated', '.'])))]
print(exome["PoluPhen-2"].value_counts())
print(exome.SIFT.value_counts())

exome_anno['gene_symbol'] = exome.Gene.unique()
gene_count = exome.Gene.value_counts().to_dict()
exome_anno['Exome'] = exome_anno.gene_symbol.map(gene_count)
print(exome_anno)
print(exome_anno.Exome.value_counts())
# Restrict to only recurrent genes
exome_anno = exome_anno[exome_anno.Exome>=2]
print(exome_anno)

df = pd.merge(df, exome_anno, on = 'gene_symbol', how = 'outer')

# Look at overlaps
df.fillna('.', inplace = True)
df['Total'] = (df[['GWAS', 'CNV', 'Exome']]!='.').sum(1)
print(df)
print(df.Total.value_counts())

# We will use the union of all three sets as our data frame

# Annotate gene list
gene_list = pd.read_csv('intermediate_tables/2_dbd_ddg2p_sfari_annotations.csv')
gene_list['SZDB_schizophrenia']=0
gene_list.loc[gene_list.gene_symbol.isin(df.gene_symbol.to_list()), 'SZDB_schizophrenia']=1
print(gene_list.SZDB_schizophrenia.value_counts())
gene_list.to_csv('intermediate_tables/3_szdb_annotations.csv', index = False)
