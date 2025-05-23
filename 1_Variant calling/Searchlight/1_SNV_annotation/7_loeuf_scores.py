import pandas as pd

# Input and output files
input_file='/path/to/input/file.csv' # Use the output of script 6_annotate_gencode_genes.py here
output_file='/path/to/output/table.csv' # Note that the output of this script will be the final coding SNV annotations
loeuf_file='/path/to/parsed/gnomad/loeuf/annotations.csv'

# Load files
df = pd.read_csv(input_file)
gnomad_df = pd.read_csv(loeuf_file, sep='\t')
gnomad_df = gnomad_df.set_index('gene_id')

# Annotate variants with LOEUF score
# If a variant is annotated with multiple genes/LOEUF scores, return the lowest
def get_loeuf(gene_id):
        gene_ids = gene_id.split(';')
        gene_ids_in_gnomad = [s.split('.')[0] for s in gene_ids if s.split('.')[0] in gnomad_df.gene_id.to_list()]

        if len(gene_ids_in_gnomad) == 0:
                return np.nan

        loeuf = min(gnomad_df[gnomad_df.gene_id.isin(gene_ids_in_gnomad)]['oe_lof_upper'].to_numpy())
        return loeuf

snvs['LOEUF']=snvs.Gene_id.apply(get_loeuf)

# Save
df.to_csv(output_file, index=False)
