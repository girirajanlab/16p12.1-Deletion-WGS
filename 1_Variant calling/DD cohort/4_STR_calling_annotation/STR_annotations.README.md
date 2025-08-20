This README describes the scripts used for STR calling and annotation in the DD cohort

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- _1_run_GangSTR_dumpSTR.sh_: Use GangSTR to call STRs and dumpSTR to filter STR calls for a family
- _2_run_mergeSTR.sh_: Merge STR calls across families and apply population level filters
- _3_exclude_loci_bed.sh_: Identify loci that failed population-level filters to exclude from downstream analysis
- _4_filter_calls.sh_: Remove loci that failed population filters, identify non-reference calls, and convert to BED format 
- _5_GangSTR_nonref.sh_: Call STRs by family for non-reference loci
- _6_run_statSTR.sh_: Use statSTR to calculate STR statistics for the cohort
- _7_2SD_STR.py_: Identify samples with STR lengths > 2 SD from the cohort mean
- _8_annovar.sh_: Annotate genes affected by variants with ANNOVAR
- _9_filter_annovar.py_: Apply ANNOVAR annotations to calls and restrict calls to exonic
- _10_annotate_gencode_genes.py_: Annotate variants with GENCODE genes
- _11_annotate_loeuf.py_: Annotate variants with LOEUF scores
- _12_inheritance_annotation.py_: Annotate variants with parent-of-origin

