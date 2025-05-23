import pandas as pd

# Input and output files
CNVNATOR="/path/to/CNVnator/outputs.txt" # Use the output from cnvnator/large_cnv_processing/13_explode_genes.py here
PENNCNV="/path/to/PennCNV/outputs.txt" # Use the output from PennCNV/11_explode_genes.py here
SMALL_CNV="/path/to/merged/small/cnv/calls.txt" # Use the output from small_cnv_merge/11_explode_genes.py here
OUTPUT="/path/to/output/file.txt"

# Combine all calls together
df=pd.DataFrame()
for i in range(3):
	pipeline=[CNVNATOR, PENNCNV, SMALL_CNV]
	pdf=pd.read_csv(pipeline[i], sep='\t')
	df=pd.concat([df, pdf])

# Replace any NaN with .
df.loc[df.Pathogenic_Name.isnull(), 'Pathogenic_Name']='.'

# Drop duplicate genes
df.drop_duplicates(subset=['Sample', 'Type', 'Gene_ID', 'Gene_Symbol'], keep = 'first')

# Save
df.to_csv(OUTPUT, sep='\t')
