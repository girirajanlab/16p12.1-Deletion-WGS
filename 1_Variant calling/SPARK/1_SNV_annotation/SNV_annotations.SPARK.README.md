This README describes the scripts used for SNV annotation in the SPARK cohort

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- _1_annotate_variants.sh_: Annotate variants with functional consequences and deleteriousness metrics
- _2_somalier_ancestry.sh_: Calculate ancestry for samples using somalier
- _3_combine_somalier_reported.py_: Merge somalier-calculated ancestries with reported ancestry
- _4_filter_pop_frequency.sh_: Filter variants for those with ancestry-specific allele frequencies < 0.001 and save as table
- _5_filter_variants.py_: Filter variants for missense or splice variants with CADD>=25 or LOF variants and variants occur in <= 10 individuals in the cohort
- _6_annotate_gencode_genes.py_: Annotate variants with GENCODE genes
- _7_loeuf_scores.py_: Annotate variants with LOEUF scores
- _8_burden_table.py_: Create a table of rare variant counts

# Files
This lists all of the non-script annotation files in the directory and subdirectories and a description.
- toml_files/
	- _gnomad_annotation_file.toml_: A TOML file describing annotation of gnomAD allele frequencies for VCFanno
	- _cadd_annotation_file.toml_: A TOML file describing annotation of CADD scores for VCFanno