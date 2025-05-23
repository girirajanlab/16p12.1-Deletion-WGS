This README describes the scripts used for CNV calling and annotation in the MyCode cohort

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- _1_run_PennCNV.sh_: Call CNVs from microarray using PennCNV
- _2_check_QC.py_: Compile CNV calling QC
- _3_merge_calls.py_: Merge adjacent CNVs per sample
- _4_size_filter.py_: Remove CNVs less than 50kb or that encompass less than 5 SNPs
- _5_frequency_pathogenicity_anno.sh_: Annotate calls with their intracohort frquency, frequency in a control cohort, and overlap with known dosage-sensitive regions
- _6_frequency_filter.py_: Filter calls to only those present at less than 0.1% frequeny in a control cohort
- _7_annotate_gencode.py_: Annotate CNVs with genes they affect
- _8_explode_genes.py_: Explode CNV calls by gene
- _9_annotate_loeuf.py_: Annotate CNV genes with LOEUF score

