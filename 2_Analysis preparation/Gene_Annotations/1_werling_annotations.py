import pandas as pd

# Input files
GENCODE="/path/to/parsed/gencode/annotations.csv"

# Organize post-hoc gene annotations
# First, get a list of all the protein-coding genes from the GENCODE V19 file
gencode = pd.read_csv(GENCODE)
gencode = gencode[gencode.gene_type=='protein_coding']

# Take these gene names and make a new DataFrame
df = gencode[['gene_name', 'gene_id']].copy()
df.drop_duplicates(keep = 'first', inplace = True)
colnames = df.columns.to_list()
colnames[0] = 'gene_symbol'
df.columns = colnames

# Remove the version number from the ENSEMBL ID
df['gene_id'] = df.gene_id.str.split('.', expand=True)[0]

# Werling annotations
# Annotations are taken from Supp Table 6 from Werling Nature Genetics 2018
# Paper link: https://pubmed.ncbi.nlm.nih.gov/29700473/
werling = pd.read_excel('Data_Files/Werling_NatGen_2018_Supp/suppTables/SuppTable06_GeneLists.xlsx')
colnames = werling.columns.to_list()
# Remove extra columns
colnames = colnames[0:9]
werling = werling[colnames]

# Get Werling annotations
df = pd.merge(df, werling, left_on = 'gene_symbol', right_on='Genes_background', how = 'left')

# Save to file
df[['gene_id', 'gene_symbol']+colnames].to_csv('intermediate_tables/1_werling_annotations.csv', index = False)