This README describes the scripts used for STR annotation in SSC

# Scripts
This lists all of the scripts in the directory and a description
- _1_run_mergeSTR_statSTR.sh_: Merge and calculate statistics on STR calls
- _2_vcf2table.py_: Convert VCFs into tables
- _3_2SD_STR.py_: Identify samples with STR lengths > 2 SD from the cohort mean
- _4_anno_locus_info.py_: Add STR locus information from reference files
- _5_annovar.sh_: Annotate genes affected by variants with ANNOVAR
- _6_append_annotations.py_: Apply ANNOVAR annotations to calls
- _7_combine_chromosomes.sh_: Combine calls from all chromosomes
- _8_filter_calls.py_: Filter calls to include only exonic variants
- _9_annotate_gencode_genes.py_: Annotate variants with GENCODE genes
- _10_annotate_loeuf.py_: Annotate variants with LOEUF scores
- _11_identify_samples.py_: Identify samples genotyped for all chromosomes
