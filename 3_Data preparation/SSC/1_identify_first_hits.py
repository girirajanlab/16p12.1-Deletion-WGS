import pandas as pd

# Input and output files
CNVS="/path/to/CNV/calls.txt" # Use the output from script 1_Variant calling/SSC/2_CNV_annotation/3_annotate_gencode_genes.py
SNVS="/path/to/SNV/variants.csv" # Use the output of script 1_Variant calling/SSC/1_SNV_annotation/8_loeuf_scores.py

GENE_ANNOTATIONS="/path/to/gene/annotation/file.csv" # Use the output of 

SSC_SAMPLES="/path/to/SSC/samples/Individuals_by_Distribution_v15.csv" # This file can be accessed from SFARI base. Researchers can get access to SFARI resources following steps here: https://www.sfari.org/resource/sfari-base/

COHORT_OUTPUT="/path/to/output/cohort/file.csv" # This file will contain sample IDs annotated with the primary variants they have
SNV_OUTPUT="/path/to/output/SNV/file.csv" # This file will contain SNVs annotated with first hits
CNV_OUTPUT="/path/to/output/CNV/file.csv" # This file will contain CNVs annotated with first hits

# Identify probands and first hits
# CNVs
cnvs=pd.read_csv(CNVS)
cnvs=cnvs[cnvs.patientID.str.contains('.p1')]
cnvs=cnvs[~cnvs.Genes.isnull()]

# Identify the primar variant for each proband
cnvs['Length']=cnvs['Stop(hg19)']-cnvs['Start(hg19)']
cnvs['LargeRare']=cnvs.Length>=500000

# If a proband has multiple large, rare CNVs, use the largest as the primary variant
cnvs['pvid']=cnvs.patientID+'_'+cnvs.Chr+':'+cnvs['Start(hg19)'].astype(str)+'-'+cnvs['Stop(hg19)'].astype(str)+cnvs['Del/Dup']
prodf=cnvs[cnvs.LargeRare][['patientID', 'pvid', 'Length', 'Del/Dup']].sort_values(by=['patientID', 'Del/Dup', 'Length'], ascending=[True, True, False])
prodf.drop_duplicates(subset=['patientID', 'Del/Dup'], inplace=True, keep='first')

cnvs['primary']=cnvs.pvid.isin(prodf.pvid.to_list())

# SNVs
snvs=pd.read_csv(SNVS)
snvs=snvs[snvs.Sample.str.contains('.p1')]

# Annotate DBD Tier 1 genes
gene_anno=pd.read_csv(GENE_ANNOTATIONS)
snvs['DBD_Tier1']=snvs['Gene_id_'].isin(gene_anno[gene_anno.Geisinger_DBD_Tier=='1'].gene_id.to_list())

# Identify the primary variant
# If a proband has multiple SNVS in DBD Tier 1 genes, use the one with the lowest LOEUF score as the primary variant
snvs['pvid']=snvs.Sample+'_'+snvs.Chrom+':'+snvs.Pos.astype(str)+'_'+snvs.Ref+'_'+snvs.Alt
prosnv=snvs[snvs.DBD_Tier1][['Sample', 'LOEUF', 'pvid']].sort_values(by=['Sample', 'LOEUF'])
prosnv.drop_duplicates(subset=['Sample'], inplace=True, keep='first')

snvs['primary']=snvs.pvid.isin(prosnv.pvid.to_list())

# Create a frame with which probands have each kind of primary variant
ssc_pros=pd.read_csv(SSC_SAMPLES)
ssc_pros=ssc_pros[ssc_pros['SSC ID'].str.contains('.p1')]
primary_df=pd.DataFrame({'Sample':list(set(snvs.Sample.to_list()+cnvs.patientID.to_list()+ssc_pros['SSC ID'].to_list()))})

primary_df['DBD Tier 1 SNVs']=primary_df.Sample.isin(snvs[snvs.primary].Sample.to_list())
primary_df['Large rare deletions']=primary_df.Sample.isin(cnvs[(cnvs.primary) & (cnvs['Del/Dup']=='Del')].patientID.to_list())
primary_df['Large rare duplications']=primary_df.Sample.isin(cnvs[(cnvs.primary) & (cnvs['Del/Dup']=='Dup')].patientID.to_list())

primary_df['No primary variant']=~(primary_df[['DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications']].any(axis=1))

# Save
primary_df.sort_values(by='Sample', inplace=True)
primary_df.to_csv(COHORT_OUTPUT, index=False)

# Save annotated variant files
snvs.to_csv(SNV_OUTPUT, index=False)

# Split CNVs by gene
cnvs=cnvs[['familyID', 'patientID', 'Chr', 'Start(hg19)', 'Stop(hg19)', 'Length',
			'Del/Dup', 'Inheritance', 'Genes', 'LargeRare', 'primary']]
cnvs['Genes']=cnvs.Genes.str.split(';')
cnvs=cnvs.explode('Genes')

# If a proband has a CNV that overlaps the primary (on the other allele), count both as the primary
cnvs['pvid']=cnvs.patientID+'.'+cnvs.Genes+'.'+cnvs['Del/Dup']
primary_cnvs=cnvs[cnvs.primary].pvid.to_list()
cnvs.loc[cnvs.pvid.isin(primary_cnvs), 'primary']=True

# Annotate LOEUF scores
gene_anno.index=gene_anno.gene_symbol.to_list()
cnvs['LOEUF']=cnvs.Genes.map(gene_anno.LOEUF.to_dict())

# Annotate gene id
cnvs['Gene_id_']=cnvs.Genes.map(gene_anno.gene_id.to_dict())

cnvs.to_csv(CNV_OUTPUT, index=False)

