This README describes the scripts used for CNV calling and annotation in the Simons Searchlight cohort

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- _1_run_PennCNV.sh_: Run PennCNV to call CNVs from microarray data
- _2_parse_penncnv.py_: Merge individual and family CNV calls and parse PennCNV output
- _3_merge.py_: Merge adjacent calls if they overlap or have a gap less than 50kb and less than 20% of the combined CNV length and filter calls less than 50kb in length or that encompass less than 5 SNPs
- _4_control_frequency.py_: Annotate calls for their frequency in a control cohort
- _5_anno_segdup_centel.sh_: Annotate calls for their overlaps with segmental duplications and centromeres and telomeres
- _6_annotate_gencode.py_: Annotate CNVs with genes they affect
- _7_explode_genes.py_: Explode CNV calls by gene
- _8_annotate_loeuf.py_: Annotate CNV genes with LOEUF score

