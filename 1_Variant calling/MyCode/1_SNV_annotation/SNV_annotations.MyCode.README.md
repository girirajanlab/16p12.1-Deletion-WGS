This README describes the scripts used for SNV annotation in the MyCode cohort

# Scripts
This lists all of the scripts in the directory and subdirectories and a description
- _1_filter_variants.sh_: Perform QC filtering on variants
- _2_annovar.sh_: Annotate variants with functional consequences using ANNOVAR
- _3_filter_exonic.sh_: Filter variants for those with exonic and splicing functions
- _4_annotate_gnomad.sh_: Annotate variants for gnomAD allele frequencies using VCFanno
- _5_somalier_ancestry.sh_: Calculate ancestry for samples using somalier
- _6_combine_somalier_reported.py_: Merge somalier-calculated ancestries with repotred ancestries
- _7_filter_pop_frequency.sh_: Filter variants for those with ancestry-specific allele frequencies < 0.001
- _8_annotate_cadd.sh_: Annotate variants with CADD scores using VCFanno
- _9_vcf2table.sh_: Merge all VCFs across samples into a single table
- _10_filter_variant_types.py_: Filter variants for missense or splice variants with CADD>=25 or LOF variants
- _11_intracohort_filter.py_: Apply an intracohort filter such that each variant occurs in <=10 individuals in the cohort
- _12_annotate_gencode_genes.py_: Annotate variants with GENCODE genes
- _13_loeuf_scores.py_: Annotate variants with LOEUF scores

# Files
This lists all of the non-script annotation files in the directory and subdirectories and the original scripts each was derived from and a description.
- toml_files/
	- _gnomad_annotation_file.toml_: A TOML file describing annotation of gnomAD allele frequencies for VCFanno
	- _cadd_annotation_file.toml_: A TOML file describing annotation of CADD scores for VCFanno
