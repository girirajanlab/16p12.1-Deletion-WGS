#!/bin/bash

mkdir -p Data_Files
cd Data_Files

# Werling NatGen2018 annotations
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0107-y/MediaObjects/41588_2018_107_MOESM3_ESM.zip
# Change filename for easier use
mv 41588_2018_107_MOESM3_ESM.zip Werling_NatGen_2018_Supp.zip
# Unzip folder
mkdir Werling_NatGen_2018_Supp
unzip -d Werling_NatGen_2018_Supp Werling_NatGen_2018_Supp.zip

# Geisinger DBD Database
wget https://dbd.geisingeradmi.org/downloads/Full-LoF-Table-Data.csv

# Genes2Phenotype DD
wget https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz

# SFARI Gene
wget --no-check-certificate https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes
# Rename file
mv download-csv.php@api-endpoint=genes SFARI-Gene_genes_01-11-2022release_03-14-2022export.csv

# SZDB2.0
# We are using an intersection of the SZDB2.0 CNV genes, GWAS genes, and Exome sequencing
mkdir SZDB
cd SZDB
wget http://szdb.org/download/cnvgene.txt
wget http://szdb.org/download/GWAS-Genes.zip
wget http://szdb.org/download/Exome.txt
unzip GWAS-Genes.zip
# These gene lists need to be filtered and parsed - this will be done as they are added to the final table

# Epilepsy gene list
# This list is a manually created compilation of all genes from all tables in Wang Seizure 2017

# Make another folder for intermediate files
mkddir -p ../intermediate_tables
