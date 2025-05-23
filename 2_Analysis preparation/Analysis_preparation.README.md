This README describes scripts related to processing of data (outside of participant genetic or phenotypic data) relevant for analyses

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- Brain_Celltype_Expression/
	- _1_parse_Allen_data.py_: Create dataframes of genes preferrentially expressed in specific brain cell types from Allen Brain Institute data
- Gene_Annotations/
	- _0_download_annotations.sh_: Download gene annotations for identifying disease-associated genes
	- _1_werling_annotations.py_: Parse disease gene annotations from [Werling et al Nature Genetics 2018](https://pubmed.ncbi.nlm.nih.gov/29700473/)
	- _2_dbd_ddg2p_sfari.py_: Parse disease gene annotations from the DBD Genes Database, Gene2Phenotype (G2P) database, and SFARI Gene
	- _3_szdb_anno.py_: Parse schizophrenia annotations from SZDB using GWAS, CNVs, or Exome sequencing results
	- _4_epilepsy_anno.py_: Parse epilepsy-associated gene annotations from [Wang et al. Seizure 2017](https://pubmed.ncbi.nlm.nih.gov/28007376/)
	- _5_add_loeuf.py_: Add gnomAD LOEUF annotations to genes
- NPC_RNAseq/
	- _1_kallisto.sh_: Trim reads and quantify transcript abundance with kallisto
	- _2_WGCNA.R_: Perform WGCNA on RNA-seq data
- STRING_Network/
	- _1_parse_string.py_: Annotate genes and updated interaction scores in STRING DB
	- _2_string_connectivity.py_: Quantify connectivity of genes in annotated STRING network

# Files
This lists all of the additional non-script files in the directory and subdirectories and a description
- Gene_Annotations/
	- Data_Files/
		- _Wang_Seizure_2017_epilepsy_gene_list.xlsx_: A list of epilesy associated genes from [Wang et al. Seizure 2017](https://pubmed.ncbi.nlm.nih.gov/28007376/)