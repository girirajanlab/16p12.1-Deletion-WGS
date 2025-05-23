import pandas as pd

# Merge and finalize CNVnator large CNV calls

# Input and output files
input_path='path/to/cnv/files' # These are the files from script 11_annotate_gencode.py
output_file='path/to/output/file.txt'

# Merge calls into a single callset
combined=pd.DataFrame()
for cnv_type in ['cnv', 'pathogenic_cnv']:
	ctdf=pd.read_csv(f'{input_path}/{cnv_type}_annotated_genes.bed', sep='\t')
	combined=pd.concat([combined, ctdf], ignore_index = True)

# Remove any duplicate calls
combined['variant_id'] = combined['Chr']+'_'+combined['Start'].astype(str)+'_'+combined['End'].astype(str)+'_'+combined['Sample']
combined.drop_duplicates(subset='variant_id', keep = 'last', inplace = True)

# Remove chrY calls
combined = combined[combined.Chr!='chrY']

# Update CN annotations
combined['Type']=combined['Type'].map({'<DEL>':'DEL', '<DUP>':'DUP'})

# Remove any CNVs that don't overlap any exons
combined=combined[combined.gene_ids!='.']

# Reorder columns and save to file
combined = combined[['Sample', 'Chr', 'Start', 'End', 'Type', 'Name', 'Length', 'Intracohort_count', 'Microarray_count', 'microarray_freq', 'gnomADSV_AF',
                	'NEJM_Name', 'gene_ids', 'gene_names', 'variant_id']]
combined.to_csv(output_file, sep = '\t', index = False)
